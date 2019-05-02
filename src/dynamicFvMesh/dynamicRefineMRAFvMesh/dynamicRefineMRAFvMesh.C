/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dynamicRefineMRAFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceInterpolate.H"
#include "volFields.H"
#include "polyTopoChange.H"
#include "surfaceFields.H"
#include "syncTools.H"
#include "pointFields.H"
#include "sigFpe.H"
#include "cellSet.H"
#include "scalar.H"

#include "searchableSurface.H"
#include "searchableBox.H"
#include "searchableSphere.H"
#include "searchableCylinder.H"

#include "topoSetSource.H"
#include "boxToCell.H"
#include "sphereToCell.H"
#include "cylinderToCell.H"

#include "setFieldsDynamic/setFieldsDynamic.H"

// #include "WaveletDecomp.h"
// #include "Wavelet.h"
#include <cmath>
#include <sys/time.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicRefineMRAFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicRefineMRAFvMesh, IOobject);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::dynamicRefineMRAFvMesh::topParentID(label p)
{
    label nextP = meshCutter_.history().splitCells()[p].parent_;
    if( nextP < 0 )
    {
        return p;
    }
    else
    {
        return topParentID(nextP);
    }
}

// the PackedBoolList::count method would probably be faster
// since we are only checking for 'true' anyhow
Foam::label Foam::dynamicRefineMRAFvMesh::count
(
    const PackedBoolList& l,
    const unsigned int val
)
{
    label n = 0;
    forAll(l, i)
    {
        if (l.get(i) == val)
        {
            n++;
        }

        // debug also serves to get-around Clang compiler trying to optimsie
        // out this forAll loop under O3 optimisation
        if (debug)
        {
            //Info<< "n=" << n << endl;
        }
    }

    return n;
}


void Foam::dynamicRefineMRAFvMesh::readDict()
{
    dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict(typeName + "Coeffs")
    );

    List<Pair<word> > fluxVelocities = List<Pair<word> >
    (
        refineDict.lookup("correctFluxes")
    );
    // Rework into hashtable.
    correctFluxes_.resize(fluxVelocities.size());
    forAll(fluxVelocities, i)
    {
        correctFluxes_.insert(fluxVelocities[i][0], fluxVelocities[i][1]);
    }

    dumpLevel_ = Switch(refineDict.lookup("dumpLevel"));
    dumpID_ = Switch(refineDict.lookup("dumpID"));
}


Foam::autoPtr<Foam::mapPolyMesh>
Foam::dynamicRefineMRAFvMesh::decomposePolyhedra
(
)
{
    // Mesh changing engine.
    polyTopoChange meshMod(*this);

    meshCutter_.decomposePolyhedra(meshMod, protectedCell_);

    // Create mesh (with inflation), return map from old to new mesh.
    //autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);

    Info<< "Decomposed from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " cells to " << globalData().nTotalCells() << " cells." << endl;

    if (debug)
    {
        // Check map.
        for (label faceI = 0; faceI < nInternalFaces(); faceI++)
        {
            label oldFaceI = map().faceMap()[faceI];

            if (oldFaceI >= nInternalFaces())
            {
                FatalErrorInFunction
                    << "New internal face:" << faceI
                    << " fc:" << faceCentres()[faceI]
                    << " originates from boundary oldFace:" << oldFaceI
                    << abort(FatalError);
            }
        }
    }

//    // Remove the stored tet base points
//    tetBasePtIsPtr_.clear();
//    // Remove the cell tree
//    cellTreePtr_.clear();

    // Update fields
    updateMesh(map);


    // Move mesh
    /*
    pointField newPoints;
    if (map().hasMotionPoints())
    {
        newPoints = map().preMotionPoints();
    }
    else
    {
        newPoints = points();
    }
    movePoints(newPoints);
    */

    // Correct the flux for modified/added faces. All the faces which only
    // have been renumbered will already have been handled by the mapping.
    {
        const labelList& faceMap = map().faceMap();
        const labelList& reverseFaceMap = map().reverseFaceMap();

        // Storage for any master faces. These will be the original faces
        // on the coarse cell that get split into four (or rather the
        // master face gets modified and three faces get added from the master)
        labelHashSet masterFaces;

        forAll(faceMap, faceI)
        {
            label oldFaceI = faceMap[faceI];

            if (oldFaceI >= 0)
            {
                label masterFaceI = reverseFaceMap[oldFaceI];

                if (masterFaceI < 0)
                {
                    FatalErrorInFunction
                        << "Problem: should not have removed faces"
                        << " when refining."
                        << nl << "face:" << faceI << abort(FatalError);
                }
                else if (masterFaceI != faceI)
                {
                    masterFaces.insert(masterFaceI);
                }
            }
        }
        if (debug)
        {
            Pout<< "Found " << masterFaces.size() << " split faces " << endl;
        }

        HashTable<surfaceScalarField*> fluxes
        (
            lookupClass<surfaceScalarField>()
        );
        forAllIter(HashTable<surfaceScalarField*>, fluxes, iter)
        {
            if (!correctFluxes_.found(iter.key()))
            {
                WarningInFunction
                    << "Cannot find surfaceScalarField " << iter.key()
                    << " in user-provided flux mapping table "
                    << correctFluxes_ << endl
                    << "    The flux mapping table is used to recreate the"
                    << " flux on newly created faces." << endl
                    << "    Either add the entry if it is a flux or use ("
                    << iter.key() << " none) to suppress this warning."
                    << endl;
                continue;
            }

            const word& UName = correctFluxes_[iter.key()];

            if (UName == "none")
            {
                continue;
            }

            if (UName == "NaN")
            {
                Pout<< "Setting surfaceScalarField " << iter.key()
                    << " to NaN" << endl;

                surfaceScalarField& phi = *iter();

                sigFpe::fillNan(phi.internalField());

                continue;
            }

            if (debug)
            {
                Pout<< "Mapping flux " << iter.key()
                    << " using interpolated flux " << UName
                    << endl;
            }

            surfaceScalarField& phi = *iter();
            const surfaceScalarField phiU
            (
                fvc::interpolate
                (
                    lookupObject<volVectorField>(UName)
                )
              & Sf()
            );

            // Recalculate new internal faces.
            for (label faceI = 0; faceI < nInternalFaces(); faceI++)
            {
                label oldFaceI = faceMap[faceI];

                if (oldFaceI == -1)
                {
                    // Inflated/appended
                    phi[faceI] = phiU[faceI];
                }
                else if (reverseFaceMap[oldFaceI] != faceI)
                {
                    // face-from-masterface
                    phi[faceI] = phiU[faceI];
                }
            }

            // Recalculate new boundary faces.
            surfaceScalarField::GeometricBoundaryField& bphi =
                phi.boundaryField();
            forAll(bphi, patchI)
            {
                fvsPatchScalarField& patchPhi = bphi[patchI];
                const fvsPatchScalarField& patchPhiU =
                    phiU.boundaryField()[patchI];

                label faceI = patchPhi.patch().start();

                forAll(patchPhi, i)
                {
                    label oldFaceI = faceMap[faceI];

                    if (oldFaceI == -1)
                    {
                        // Inflated/appended
                        patchPhi[i] = patchPhiU[i];
                    }
                    else if (reverseFaceMap[oldFaceI] != faceI)
                    {
                        // face-from-masterface
                        patchPhi[i] = patchPhiU[i];
                    }

                    faceI++;
                }
            }

            // Update master faces
            forAllConstIter(labelHashSet, masterFaces, iter)
            {
                label faceI = iter.key();

                if (isInternalFace(faceI))
                {
                    phi[faceI] = phiU[faceI];
                }
                else
                {
                    label patchI = boundaryMesh().whichPatch(faceI);
                    label i = faceI - boundaryMesh()[patchI].start();

                    const fvsPatchScalarField& patchPhiU =
                        phiU.boundaryField()[patchI];

                    fvsPatchScalarField& patchPhi = bphi[patchI];

                    patchPhi[i] = patchPhiU[i];
                }
            }
        }
    }

    // Update numbering of cells/vertices.
    meshCutter_.updateMesh(map);

    // Update numbering of protectedCell_
    if (protectedCell_.size())
    {
        PackedBoolList newProtectedCell(nCells());

        forAll(newProtectedCell, cellI)
        {
            label oldCellI = map().cellMap()[cellI];
            newProtectedCell.set(cellI, protectedCell_.get(oldCellI));
        }
        protectedCell_.transfer(newProtectedCell);
    }

    // Debug: Check refinement levels (across faces only)
    meshCutter_.checkRefinementLevels(-1, labelList(0));

    // check for zero volume cells
    label nZeroVolCells = 0;

    const volScalarField::DimensionedInternalField& vols = V();

    forAll(vols, i)
    {
        if(vols.field()[i] == 0)
        {
            nZeroVolCells++;

            Info << "Zero volume cell: " << i << nl;
        }
    }

    if(nZeroVolCells > 0)
    {
        FatalErrorInFunction
            << "Found " << nZeroVolCells << " zero volume cells."
            << abort(FatalError);
    }

    if
    (
        returnReduce(map().nOldCells(), sumOp<label>())
            < globalData().nTotalCells()
    )
    {
        topoChanging(true);

        // Reset moving flag (if any). If not using inflation we'll not move,
        // if are using inflation any follow on movePoints will set it.
        moving(false);
    }
    else
    {
        topoChanging(false);
    }

    return map;
}


// Refines cells, maps fields and recalculates (an approximate) flux
Foam::autoPtr<Foam::mapPolyMesh>
Foam::dynamicRefineMRAFvMesh::refine
(
    const labelList &cellsToRefine
)
{
    // Mesh changing engine.
    polyTopoChange meshMod(*this);

//    struct timeval timeStart, timeStop, diff;
//    label sec, usec;
//    gettimeofday(&timeStart, NULL);

    // Play refinement commands into mesh changer.
    meshCutter_.setRefinement(cellsToRefine, meshMod);

//    gettimeofday(&timeStop, NULL);
//    timersub(&timeStop, &timeStart, &diff);
//    sec = diff.tv_sec;
//    usec = diff.tv_usec;
//    Info << "refTime setRefinement time:        " << sec << "." << usec << nl;


//    gettimeofday(&timeStart, NULL);

    // Create mesh (with inflation), return map from old to new mesh.
    //autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);

    Info<< "Refined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << globalData().nTotalCells() << " cells." << endl;

    if (debug)
    {
        // Check map.
        for (label faceI = 0; faceI < nInternalFaces(); faceI++)
        {
            label oldFaceI = map().faceMap()[faceI];

            if (oldFaceI >= nInternalFaces())
            {
                FatalErrorInFunction
                    << "New internal face:" << faceI
                    << " fc:" << faceCentres()[faceI]
                    << " originates from boundary oldFace:" << oldFaceI
                    << abort(FatalError);
            }
        }
    }

//    // Remove the stored tet base points
//    tetBasePtIsPtr_.clear();
//    // Remove the cell tree
//    cellTreePtr_.clear();

    // Update fields
    updateMesh(map);


    // Move mesh
    /*
    pointField newPoints;
    if (map().hasMotionPoints())
    {
        newPoints = map().preMotionPoints();
    }
    else
    {
        newPoints = points();
    }
    movePoints(newPoints);
    */

    // Correct the flux for modified/added faces. All the faces which only
    // have been renumbered will already have been handled by the mapping.
    {
        const labelList& faceMap = map().faceMap();
        const labelList& reverseFaceMap = map().reverseFaceMap();

        // Storage for any master faces. These will be the original faces
        // on the coarse cell that get split into four (or rather the
        // master face gets modified and three faces get added from the master)
        labelHashSet masterFaces
            (4*cellsToRefine.size());

        forAll(faceMap, faceI)
        {
            label oldFaceI = faceMap[faceI];

            if (oldFaceI >= 0)
            {
                label masterFaceI = reverseFaceMap[oldFaceI];

                if (masterFaceI < 0)
                {
                    FatalErrorInFunction
                        << "Problem: should not have removed faces"
                        << " when refining."
                        << nl << "face:" << faceI << abort(FatalError);
                }
                else if (masterFaceI != faceI)
                {
                    masterFaces.insert(masterFaceI);
                }
            }
        }
        if (debug)
        {
            Pout<< "Found " << masterFaces.size() << " split faces " << endl;
        }

        HashTable<surfaceScalarField*> fluxes
        (
            lookupClass<surfaceScalarField>()
        );
        forAllIter(HashTable<surfaceScalarField*>, fluxes, iter)
        {
            if (!correctFluxes_.found(iter.key()))
            {
                WarningInFunction
                    << "Cannot find surfaceScalarField " << iter.key()
                    << " in user-provided flux mapping table "
                    << correctFluxes_ << endl
                    << "    The flux mapping table is used to recreate the"
                    << " flux on newly created faces." << endl
                    << "    Either add the entry if it is a flux or use ("
                    << iter.key() << " none) to suppress this warning."
                    << endl;
                continue;
            }

            const word& UName = correctFluxes_[iter.key()];

            if (UName == "none")
            {
                continue;
            }

            if (UName == "NaN")
            {
                Pout<< "Setting surfaceScalarField " << iter.key()
                    << " to NaN" << endl;

                surfaceScalarField& phi = *iter();

                sigFpe::fillNan(phi.internalField());

                continue;
            }

            if (debug)
            {
                Pout<< "Mapping flux " << iter.key()
                    << " using interpolated flux " << UName
                    << endl;
            }

            surfaceScalarField& phi = *iter();
            const surfaceScalarField phiU
            (
                fvc::interpolate
                (
                    lookupObject<volVectorField>(UName)
                )
              & Sf()
            );

            // Recalculate new internal faces.
            for (label faceI = 0; faceI < nInternalFaces(); faceI++)
            {
                label oldFaceI = faceMap[faceI];

                if (oldFaceI == -1)
                {
                    // Inflated/appended
                    phi[faceI] = phiU[faceI];
                }
                else if (reverseFaceMap[oldFaceI] != faceI)
                {
                    // face-from-masterface
                    phi[faceI] = phiU[faceI];
                }
            }

            // Recalculate new boundary faces.
            surfaceScalarField::GeometricBoundaryField& bphi =
                phi.boundaryField();
            forAll(bphi, patchI)
            {
                fvsPatchScalarField& patchPhi = bphi[patchI];
                const fvsPatchScalarField& patchPhiU =
                    phiU.boundaryField()[patchI];

                label faceI = patchPhi.patch().start();

                forAll(patchPhi, i)
                {
                    label oldFaceI = faceMap[faceI];

                    if (oldFaceI == -1)
                    {
                        // Inflated/appended
                        patchPhi[i] = patchPhiU[i];
                    }
                    else if (reverseFaceMap[oldFaceI] != faceI)
                    {
                        // face-from-masterface
                        patchPhi[i] = patchPhiU[i];
                    }

                    faceI++;
                }
            }

            // Update master faces
            forAllConstIter(labelHashSet, masterFaces, iter)
            {
                label faceI = iter.key();

                if (isInternalFace(faceI))
                {
                    phi[faceI] = phiU[faceI];
                }
                else
                {
                    label patchI = boundaryMesh().whichPatch(faceI);
                    label i = faceI - boundaryMesh()[patchI].start();

                    const fvsPatchScalarField& patchPhiU =
                        phiU.boundaryField()[patchI];

                    fvsPatchScalarField& patchPhi = bphi[patchI];

                    if ((patchPhi.size() == 0) || (patchPhiU.size() == 0))
                    {
                        WarningInFunction
                            << "Patch " << patchI << ":" << nl
                            << "  patchPhi size: " << patchPhi.size()
                            << "  patchPhiU size: " << patchPhiU.size()
                            << nl;
                    }
                    else
                    {
                        patchPhi[i] = patchPhiU[i];
                    }
                }
            }
        }
    }



//    PackedBoolList newRefineCell(map().cellMap().size(), false);
//    PackedBoolList newUnrefineCell(map().cellMap().size(), false);
//    PackedBoolList newIsCellMRADetailed(map().cellMap().size(), false);
//    forAll(refineCell, i)
//    {
//        newRefineCell[map().reverseCellMap()[i]] = refineCell[i];
//        newUnrefineCell[map().reverseCellMap()[i]] = unrefineCell[i];
//        newIsCellMRADetailed[map().reverseCellMap()[i]] = isCellMRADetailed_[i];
//    }
//    refineCell.transfer(newRefineCell);
//    unrefineCell.transfer(newUnrefineCell);
//    isCellMRADetailed_.transfer(newIsCellMRADetailed);



    // Update numbering of cells/vertices.
    meshCutter_.updateMesh(map);

    // Update numbering of protectedCell_
    if (protectedCell_.size())
    {
        PackedBoolList newProtectedCell(nCells());

        forAll(newProtectedCell, cellI)
        {
            label oldCellI = map().cellMap()[cellI];
            newProtectedCell.set(cellI, protectedCell_.get(oldCellI));
        }
        protectedCell_.transfer(newProtectedCell);
    }

    // Debug: Check refinement levels (across faces only)
    meshCutter_.checkRefinementLevels(-1, labelList(0));

//    gettimeofday(&timeStop, NULL);
//    timersub(&timeStop, &timeStart, &diff);
//    sec = diff.tv_sec;
//    usec = diff.tv_usec;
//    Info << "refTime updateMeshRefinement: " << sec << "." << usec << nl;

    return map;
}


// Combines previously split cells, maps fields and recalculates
// (an approximate) flux
Foam::autoPtr<Foam::mapPolyMesh>
Foam::dynamicRefineMRAFvMesh::unrefine
(
    const labelList& splitPoints,
    const labelList& unrefineTetCellsParent
)
{
    polyTopoChange meshMod(*this);

//    struct timeval timeStart, timeStop, diff;
//    label sec, usec;
//    gettimeofday(&timeStart, NULL);

    // Play refinement commands into mesh changer.
    Map<label> faceToCenterFace =
    meshCutter_.setUnrefinement(splitPoints, unrefineTetCellsParent, meshMod);

//    gettimeofday(&timeStop, NULL);
//    timersub(&timeStop, &timeStart, &diff);
//    sec = diff.tv_sec;
//    usec = diff.tv_usec;
//    Info << "refTime setUnrefinement time:      " << sec << "." << usec << nl;



//    gettimeofday(&timeStart, NULL);

    // Save information on faces that will be combined
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Find the faceMidPoints on cells to be combined.
    // for each face resulting of split of face into four store the
    // midpoint
    Map<label> faceToSplitPoint(3*splitPoints.size());

    {
        forAll(splitPoints, i)
        {
            label pointI = splitPoints[i];

            const labelList& pEdges = pointEdges()[pointI];

            forAll(pEdges, j)
            {
                label otherPointI = edges()[pEdges[j]].otherVertex(pointI);

                const labelList& pFaces = pointFaces()[otherPointI];

                forAll(pFaces, pFaceI)
                {
                    faceToSplitPoint.insert(pFaces[pFaceI], otherPointI);
                }
            }
        }
    }


    // Change mesh and generate map.
    //autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);

    Info<< "Unrefined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << globalData().nTotalCells() << " cells."
        << endl;

    // Update fields
    updateMesh(map);


    // Move mesh
    /*
    pointField newPoints;
    if (map().hasMotionPoints())
    {
        newPoints = map().preMotionPoints();
    }
    else
    {
        newPoints = points();
    }
    movePoints(newPoints);
    */

    // Correct the flux for modified faces.
    {
        const labelList& reversePointMap = map().reversePointMap();
        const labelList& reverseFaceMap = map().reverseFaceMap();

        HashTable<surfaceScalarField*> fluxes
        (
            lookupClass<surfaceScalarField>()
        );
        forAllIter(HashTable<surfaceScalarField*>, fluxes, iter)
        {
            if (!correctFluxes_.found(iter.key()))
            {
                WarningInFunction
                    << "Cannot find surfaceScalarField " << iter.key()
                    << " in user-provided flux mapping table "
                    << correctFluxes_ << endl
                    << "    The flux mapping table is used to recreate the"
                    << " flux on newly created faces." << endl
                    << "    Either add the entry if it is a flux or use ("
                    << iter.key() << " none) to suppress this warning."
                    << endl;
                continue;
            }

            const word& UName = correctFluxes_[iter.key()];

            if (UName == "none")
            {
                continue;
            }

            if (debug)
            {
                Info<< "Mapping flux " << iter.key()
                    << " using interpolated flux " << UName
                    << endl;
            }

            surfaceScalarField& phi = *iter();
            surfaceScalarField::GeometricBoundaryField& bphi =
                phi.boundaryField();

            const surfaceScalarField phiU
            (
                fvc::interpolate
                (
                    lookupObject<volVectorField>(UName)
                )
              & Sf()
            );


            forAllConstIter(Map<label>, faceToSplitPoint, iter)
            {
                label oldFaceI = iter.key();
                label oldPointI = iter();

                if (reversePointMap[oldPointI] < 0)
                {
                    // midpoint was removed. See if face still exists.
                    label faceI = reverseFaceMap[oldFaceI];

                    if (faceI >= 0)
                    {
                        if (isInternalFace(faceI))
                        {
                            phi[faceI] = phiU[faceI];
                        }
                        else
                        {
                            label patchI = boundaryMesh().whichPatch(faceI);
                            label i = faceI - boundaryMesh()[patchI].start();

                            const fvsPatchScalarField& patchPhiU =
                                phiU.boundaryField()[patchI];
                            fvsPatchScalarField& patchPhi = bphi[patchI];
                            patchPhi[i] = patchPhiU[i];
                        }
                    }
                }
            }

            forAllConstIter(Map<label>, faceToCenterFace, iter)
            {
                label oldFaceI = iter.key();
                label oldCenterFaceI = iter();
                label centerFaceI = reverseFaceMap[oldCenterFaceI];

                if (centerFaceI < 0)
                {
                    // centerFace was removed. See if face still exists.
                    label faceI = reverseFaceMap[oldFaceI];

                    if (faceI >= 0)
                    {
                        if (isInternalFace(faceI))
                        {
                            phi[faceI] = phiU[faceI];
                        }
                        else
                        {
                            label patchI = boundaryMesh().whichPatch(faceI);
                            label i = faceI - boundaryMesh()[patchI].start();

                            const fvsPatchScalarField& patchPhiU =
                                phiU.boundaryField()[patchI];
                            fvsPatchScalarField& patchPhi = bphi[patchI];
                            patchPhi[i] = patchPhiU[i];
                        }
                    }
                }
                else
                {
                    if (isInternalFace(centerFaceI))
                    {
                        phi[centerFaceI] = phiU[centerFaceI];
                    }
                    else
                    {
                        label patchI = boundaryMesh().whichPatch(centerFaceI);
                        label i = centerFaceI - boundaryMesh()[patchI].start();

                        const fvsPatchScalarField& patchPhiU =
                            phiU.boundaryField()[patchI];
                        fvsPatchScalarField& patchPhi = bphi[patchI];
                        patchPhi[i] = patchPhiU[i];
                    }
                }
            }
        }
    }


//    PackedBoolList newRefineCell(map().cellMap().size(), false);
//    PackedBoolList newUnrefineCell(map().cellMap().size(), false);
//    PackedBoolList newIsCellMRADetailed(map().cellMap().size(), false);
//    forAll(newRefineCell, i)
//    {
//        newRefineCell[i] = refineCell[map().cellMap()[i]];
//        newUnrefineCell[i] = unrefineCell[map().cellMap()[i]];
//        newIsCellMRADetailed[i] = isCellMRADetailed_[map().cellMap()[i]];
//    }
//    refineCell.transfer(newRefineCell);
//    unrefineCell.transfer(newUnrefineCell);
//    isCellMRADetailed_.transfer(newIsCellMRADetailed);


    // Update numbering of cells/vertices.
    meshCutter_.updateMesh(map);

    // Update numbering of protectedCell_
    if (protectedCell_.size())
    {
        PackedBoolList newProtectedCell(nCells());

        forAll(newProtectedCell, cellI)
        {
            label oldCellI = map().cellMap()[cellI];
            if (oldCellI >= 0)
            {
                newProtectedCell.set(cellI, protectedCell_.get(oldCellI));
            }
        }
        protectedCell_.transfer(newProtectedCell);
    }

    // Debug: Check refinement levels (across faces only)
    meshCutter_.checkRefinementLevels(-1, labelList(0));

//    gettimeofday(&timeStop, NULL);
//    timersub(&timeStop, &timeStart, &diff);
//    sec = diff.tv_sec;
//    usec = diff.tv_usec;
//    Info << "refTime updateMeshUnrefine time:   " << sec << "." << usec << nl;

    return map;
}


// Get max of connected point
Foam::scalarField
Foam::dynamicRefineMRAFvMesh::maxPointField(const scalarField& pFld) const
{
    scalarField vFld(nCells(), -GREAT);

    forAll(pointCells(), pointI)
    {
        const labelList& pCells = pointCells()[pointI];

        forAll(pCells, i)
        {
            vFld[pCells[i]] = max(vFld[pCells[i]], pFld[pointI]);
        }
    }
    return vFld;
}


// Get max of connected cell
Foam::scalarField
Foam::dynamicRefineMRAFvMesh::maxCellField(const volScalarField& vFld) const
{
    scalarField pFld(nPoints(), -GREAT);

    forAll(pointCells(), pointI)
    {
        const labelList& pCells = pointCells()[pointI];

        forAll(pCells, i)
        {
            pFld[pointI] = max(pFld[pointI], vFld[pCells[i]]);
        }
    }
    return pFld;
}


// Simple (non-parallel) interpolation by averaging.
Foam::scalarField
Foam::dynamicRefineMRAFvMesh::cellToPoint(const scalarField& vFld) const
{
    scalarField pFld(nPoints());

    forAll(pointCells(), pointI)
    {
        const labelList& pCells = pointCells()[pointI];

        scalar sum = 0.0;
        forAll(pCells, i)
        {
            sum += vFld[pCells[i]];
        }
        pFld[pointI] = sum/pCells.size();
    }
    return pFld;
}


// Calculate error. Is < 0 or distance to minLevel, maxLevel
Foam::scalarField Foam::dynamicRefineMRAFvMesh::error
(
    const scalarField& fld,
    const scalar minLevel,
    const scalar maxLevel
) const
{
    scalarField c(fld.size(), -1);

    forAll(fld, i)
    {
        scalar err = min(fld[i]-minLevel, maxLevel-fld[i]);

        if (err >= 0)
        {
            c[i] = err;
        }
    }
    return c;
}


void Foam::dynamicRefineMRAFvMesh::selectRefineCandidates
(
    const scalar lowerRefineLevel,
    const scalar upperRefineLevel,
    const scalarField& vFld,
    PackedBoolList& candidateCell
) const
{
    // Get error per cell. Is -1 (not to be refined) to >0 (to be refined,
    // higher more desirable to be refined).
    scalarField cellError
    (
        maxPointField
        (
            error
            (
                cellToPoint(vFld),
                lowerRefineLevel,
                upperRefineLevel
            )
        )
    );

    // Mark cells that are candidates for refinement.
    forAll(cellError, cellI)
    {
        if (cellError[cellI] > 0)
        {
            candidateCell.set(cellI, 1);
        }
    }
}


Foam::labelList Foam::dynamicRefineMRAFvMesh::selectRefineCells
(
    const label maxCells,
    const label maxRefinement,
    const PackedBoolList& candidateCell,
    const bool maxSet
) const
{
    // Every refined cell causes 7 extra cells
    label nTotToRefine = (maxCells - globalData().nTotalCells()) / 7;

    const labelList& cellLevel = meshCutter_.cellLevel();

    // Mark cells that cannot be refined since they would trigger refinement
    // of protected cells (since 2:1 cascade)
    PackedBoolList unrefineableCell;

    // Count current selection
    label nLocalCandidates = count(candidateCell, 1);
    label nCandidates = returnReduce(nLocalCandidates, sumOp<label>());

    // Collect all cells
    DynamicList<label> candidates(nLocalCandidates);

    if (nCandidates < nTotToRefine)
    {
        forAll(candidateCell, cellI)
        {
            if
            (
                cellLevel[cellI] < maxRefinement
             && candidateCell.get(cellI)
             && (
                    unrefineableCell.empty()
                 || !unrefineableCell.get(cellI)
                )
            )
            {
                candidates.append(cellI);
            }
        }
    }
    else
    {
        // Sort by error? For now just truncate.
        for (label level = 0; level < maxRefinement; level++)
        {
            forAll(candidateCell, cellI)
            {
                if
                (
                    cellLevel[cellI] == level
                 && candidateCell.get(cellI)
                 && (
                        unrefineableCell.empty()
                     || !unrefineableCell.get(cellI)
                    )
                )
                {
                    candidates.append(cellI);
                }
            }

            if (returnReduce(candidates.size(), sumOp<label>()) > nTotToRefine)
            {
                break;
            }
        }
    }

    // Guarantee 2:1 refinement after refinement
    labelList consistentSet
    (
        meshCutter_.consistentRefinement
        (
            candidates.shrink(),
            maxSet               // Add to set to guarantee 2:1
        )
    );

    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " cells for refinement out of " << globalData().nTotalCells()
        << "." << endl;

    return consistentSet;
}


void Foam::dynamicRefineMRAFvMesh::extendMarkedCells
(
    PackedBoolList& markedCell
) const
{
    // Mark faces using any marked cell
    boolList markedFace(nFaces(), false);

    forAll(markedCell, cellI)
    {
        if (markedCell.get(cellI))
        {
            const cell& cFaces = cells()[cellI];

            forAll(cFaces, i)
            {
                markedFace[cFaces[i]] = true;
            }
        }
    }

    syncTools::syncFaceList(*this, markedFace, orEqOp<bool>());

    // Update cells using any markedFace
    for (label faceI = 0; faceI < nInternalFaces(); faceI++)
    {
        if (markedFace[faceI])
        {
            markedCell.set(faceOwner()[faceI], 1);
            markedCell.set(faceNeighbour()[faceI], 1);
        }
    }
    for (label faceI = nInternalFaces(); faceI < nFaces(); faceI++)
    {
        if (markedFace[faceI])
        {
            markedCell.set(faceOwner()[faceI], 1);
        }
    }
}


void Foam::dynamicRefineMRAFvMesh::checkEightAnchorPoints
(
    PackedBoolList& protectedCell,
    label& nProtected
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();

    labelList nAnchorPoints(nCells(), 0);

    forAll(pointLevel, pointI)
    {
        const labelList& pCells = pointCells(pointI);

        forAll(pCells, pCellI)
        {
            label cellI = pCells[pCellI];

            if (pointLevel[pointI] <= cellLevel[cellI])
            {
                // Check if cell has already 8 anchor points -> protect cell
                if (nAnchorPoints[cellI] == 8)
                {
                    if (protectedCell.set(cellI, true))
                    {
                        nProtected++;
                    }
                }

                if (!protectedCell[cellI])
                {
                    nAnchorPoints[cellI]++;
                }
            }
        }
    }

    forAll(protectedCell, cellI)
    {
        if
        (
            !protectedCell[cellI]
         && nAnchorPoints[cellI] != 8
         && nAnchorPoints[cellI] != 4
        )
        {
            protectedCell.set(cellI, true);
            nProtected++;
        }
    }
}


Foam::label Foam::dynamicRefineMRAFvMesh::identifyProtectedCells()
{
    protectedCell_.reset();

    protectedCell_.resize(nCells());

    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();

    // Set cells that should not be refined.
    // This is currently any cell which does not have 8 anchor points or
    // uses any face which does not have 4 anchor points.
    // Note: do not use cellPoint addressing

    // Count number of points <= cellLevel
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList nAnchors(nCells(), 0);

    label nProtected = 0;

    forAll(pointCells(), pointI)
    {
        const labelList& pCells = pointCells()[pointI];

        forAll(pCells, i)
        {
            label cellI = pCells[i];

            if (!protectedCell_.get(cellI))
            {
                if (pointLevel[pointI] <= cellLevel[cellI])
                {
                    nAnchors[cellI]++;

                    if (nAnchors[cellI] > 8)
                    {
                        protectedCell_.set(cellI, 1);
                        nProtected++;
                    }
                }
            }
        }
    }


    // Count number of points <= faceLevel
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Bit tricky since proc face might be one more refined than the owner since
    // the coupled one is refined.

    {
        labelList neiLevel(nFaces());

        for (label faceI = 0; faceI < nInternalFaces(); faceI++)
        {
            neiLevel[faceI] = cellLevel[faceNeighbour()[faceI]];
        }
        for (label faceI = nInternalFaces(); faceI < nFaces(); faceI++)
        {
            neiLevel[faceI] = cellLevel[faceOwner()[faceI]];
        }
        syncTools::swapFaceList(*this, neiLevel);


        boolList protectedFace(nFaces(), false);

        forAll(faceOwner(), faceI)
        {
            label faceLevel = max
            (
                cellLevel[faceOwner()[faceI]],
                neiLevel[faceI]
            );

            const face& f = faces()[faceI];

            label nAnchors = 0;

            forAll(f, fp)
            {
                if (pointLevel[f[fp]] <= faceLevel)
                {
                    nAnchors++;

                    if (nAnchors > 4)
                    {
                        protectedFace[faceI] = true;
                        break;
                    }
                }
            }
        }

        syncTools::syncFaceList(*this, protectedFace, orEqOp<bool>());

        for (label faceI = 0; faceI < nInternalFaces(); faceI++)
        {
            if (protectedFace[faceI])
            {
                protectedCell_.set(faceOwner()[faceI], 1);
                nProtected++;
                protectedCell_.set(faceNeighbour()[faceI], 1);
                nProtected++;
            }
        }
        for (label faceI = nInternalFaces(); faceI < nFaces(); faceI++)
        {
            if (protectedFace[faceI])
            {
                protectedCell_.set(faceOwner()[faceI], 1);
                nProtected++;
            }
        }

        // Check cells for 8 corner points
        checkEightAnchorPoints(protectedCell_, nProtected);
    }

    cellSet protectedCells(*this, "protectedCells", nProtected);
    forAll(protectedCell_, cellI)
    {
        if (protectedCell_[cellI])
        {
            protectedCells.insert(cellI);
        }
    }

    Info<< "Detected " << returnReduce(nProtected, sumOp<label>())
        << " cells that are protected from refinement."
        << " Writing these to cellSet "
        << protectedCells.name()
        << "." << endl;

    protectedCells.write();

    return nProtected;
}

//- Check if cell has a parent
bool Foam::dynamicRefineMRAFvMesh::isCellTreatable
(
    label cellI
)
{
    label visibleCell = meshCutter_.history().visibleCells()[cellI];
    label splitCellParent = -2;

    if (visibleCell >= 0)
    {
        splitCellParent = meshCutter_.history().splitCells()[visibleCell].parent_;
        if (splitCellParent < 0)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
    else
    {
        return false;
    }
}


//- load balance the mesh
void Foam::dynamicRefineMRAFvMesh::balanceMesh
(
    const Foam::floatScalar allowableImbalance
)
{
    // source: https://github.com/tgvoskuilen/meshBalancing

    //First determine current level of imbalance - do this for all
    // parallel runs with a changing mesh, even if balancing is disabled
    label nGlobalCells =
        globalData().nTotalCells();

    floatScalar idealNCells =
        floatScalar(nGlobalCells)/floatScalar(Pstream::nProcs());

    floatScalar localImbalance =
        mag(nCells() - idealNCells);

    Foam::reduce(localImbalance, maxOp<floatScalar>());

    floatScalar maxImbalance = localImbalance/idealNCells;

    Info << "Maximum imbalance = " << 100 * maxImbalance << " %" << endl;

    // If imbalanced, construct weighted coarse graph (level 0) with node
    // weights equal to their number of subcells. This partitioning works
    // as long as the number of level 0 cells is several times greater than
    // the number of processors.

    if
    (
        (maxImbalance > allowableImbalance)
    )
    {
        Info << "Re-balancing dynamically refined mesh" << endl;

        const labelIOList& cellLevel = meshCutter_.cellLevel();
        const refinementMRAHistory &history = meshCutter_.history();

        Map<label> coarseIDmap(100);
        labelList uniqueIndex(nCells(), 0);

        label nCoarse = 0;

        const labelList &visibleCells = history.visibleCells();

        forAll(cells(), cellI)
        {
            if
            (
                (visibleCells[cellI] < 0)
             || (
                    visibleCells[cellI] >= 0
                 && history.parentIndex(cellI) < 0
                )
            )
            {
                uniqueIndex[cellI] = cellI;
            }
            else
            {
                uniqueIndex[cellI] = nCells() + topParentID
                (
                    history.parentIndex(cellI)
                );
            }

            if( coarseIDmap.insert(uniqueIndex[cellI], nCoarse) )
            {
                ++nCoarse;
            }
        }

        // Convert to local sequential indexing and calculate coarse
        // points and weights
        labelList localIndex(nCells(),0);
        pointField coarsePoints(nCoarse,vector::zero);
        scalarField coarseWeights(nCoarse,0.0);

        forAll(uniqueIndex, cellI)
        {
            localIndex[cellI] = coarseIDmap[uniqueIndex[cellI]];

            // If 2D refinement (quadtree) is ever implemented, this '3'
            // should be set in general as the number of refinement
            // dimensions.
            label w = (1 << (3*cellLevel[cellI]));

            coarseWeights[localIndex[cellI]] += 1.0;
            coarsePoints[localIndex[cellI]] += C()[cellI]/w;
        }

        //Set up decomposer - a separate dictionary is used here so
        // you can use a simple partitioning for decomposePar and
        // ptscotch for the rebalancing (or any chosen algorithms)
        autoPtr<decompositionMethod> decomposer
        (
            decompositionMethod::New
            (
                IOdictionary
                (
                    IOobject
                    (
                        "balanceParDict",
                        time().system(),
                        *this,
                        IOobject::MUST_READ_IF_MODIFIED,
                        IOobject::NO_WRITE
                    )
                )
            )
        );

        labelList finalDecomp = decomposer().decompose
        (
            *this,
            localIndex,
            coarsePoints,
            coarseWeights
        );

        scalar tolDim = globalMeshData::matchTol_ * bounds().mag();


        fvMeshDistribute distributor(*this, tolDim);

        autoPtr<mapDistributePolyMesh> map =
              distributor.distribute(finalDecomp);

        meshCutter_.distribute(map);

        //Correct values on all cyclic patches

        correctBoundaries<scalar>();
        correctBoundaries<vector>();
        correctBoundaries<sphericalTensor>();
        correctBoundaries<symmTensor>();
        correctBoundaries<tensor>();
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicRefineMRAFvMesh::dynamicRefineMRAFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    meshCutter_(*this),
    dumpLevel_(false),
    nRefinementIterations_(0),
    protectedCell_(nCells(), 0),
    isCellTreated_(nCells(), 0),
    isCellMRADetailed_(nCells(), 0),
    maxCells_(0),
    maxRefinement_(0),
    isInitialized_(false)
{
    // Read static part of dictionary
    readDict();

    label nProtected = identifyProtectedCells();

    if (nProtected > 0)
    {
        // Check if there are any cells that are not hex and tet. If so
        // then decompose them into tets.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        decomposePolyhedra();
    }
}


bool Foam::dynamicRefineMRAFvMesh::conductRefinement
(
    PackedBoolList &refineCell,
    PackedBoolList &unrefineCell
)
{
//    struct timeval timeStart, timeStop, diff;
//    label sec, usec;
//    gettimeofday(&timeStart, NULL);

    labelList cellsToRefine
    (
        selectRefineCells
        (
            maxCells_,
            maxRefinement_,
            refineCell
        )
    );

//    gettimeofday(&timeStop, NULL);
//    timersub(&timeStop, &timeStart, &diff);
//    sec = diff.tv_sec;
//    usec = diff.tv_usec;
//    Info << "refTime consistentRefinement:      " << sec << "." << usec << nl;


    label nCellsToRefine = returnReduce
    (
        cellsToRefine.size(), sumOp<label>()
    );

    if (nCellsToRefine > 0)
    {
        // Refine/update mesh and map fields
        autoPtr<mapPolyMesh> map =
            refine(cellsToRefine);

        const labelList& cellMap = map().cellMap();
        const labelList& reverseCellMap = map().reverseCellMap();

        // Update refineCell. Note that some of the marked ones have
        // not been refined due to constraints.
        {
            PackedBoolList newRefineCell(cellMap.size());

            forAll(cellMap, cellI)
            {
                label oldCellI = cellMap[cellI];

                if (oldCellI < 0)
                {
                    newRefineCell.set(cellI, 1);
                }
                else if (reverseCellMap[oldCellI] != cellI)
                {
                    newRefineCell.set(cellI, 1);
                }
                else
                {
                    newRefineCell.set(cellI, refineCell.get(oldCellI));
                }
            }
            refineCell.transfer(newRefineCell);
        }

        // Update unrefineCell.
        if (unrefineCell.size() > 0)
        {
            const labelList& cellMap = map().cellMap();
            const labelList& reverseCellMap = map().reverseCellMap();

            PackedBoolList newUnrefineCell(cellMap.size());

            forAll(cellMap, cellI)
            {
                label oldCellI = cellMap[cellI];

                if (oldCellI < 0)
                {
                    newUnrefineCell.set(cellI, 0);
                }
                else if (reverseCellMap[oldCellI] != cellI)
                {
                    newUnrefineCell.set(cellI, 0);
                }
                else
                {
                    newUnrefineCell.set(cellI, unrefineCell.get(oldCellI));
                }
            }
            unrefineCell.transfer(newUnrefineCell);
        }

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::dynamicRefineMRAFvMesh::conductUnrefinement
(
    PackedBoolList &unrefineCell
)
{
    labelList unrefineHexSplitPoints;
    labelList unrefineTetCellsParent;
    label nTotalSplitPoints;

    {
        labelList pointsToUnrefine;
        labelList tetCellsToUnrefine;

        {
            PackedBoolList unrefineHexCell(nCells());
            PackedBoolList unrefineTetCell(nCells());

            for(label cellI = 0; cellI < nCells(); cellI++)
            {
                if (unrefineCell.get(cellI))
                {
                    label nAnchors = 0;

                    const labelList &cPoints = cellPoints()[cellI];

                    forAll(cPoints, i)
                    {
                        label pointI = cPoints[i];

                        if
                        (
                            meshCutter_.pointLevel()[pointI]
                                <= meshCutter_.cellLevel()[cellI]
                        )
                        {
                            nAnchors++;
                        }
                    }

                    if (nAnchors == 4)
                    {
                        unrefineTetCell.set(cellI, 1);
                    }
                    else if (nAnchors == 8)
                    {
                        unrefineHexCell.set(cellI, 1);
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "Cell " << cellI << " is neither hex nor tet. "
                            << "It has " << nAnchors << " anchors."
                            << abort(FatalError);
                    }
                }
            }

            // HEX: get candidate split points
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            const labelList splitPoints
                (meshCutter_.getSplitPoints());

            nTotalSplitPoints =
                returnReduce(splitPoints.size(), sumOp<label>());

            const labelListList& pointCells =
                this->pointCells();

            DynamicList<label> splitPointCandidates(splitPoints.size());

            forAll(splitPoints, i)
            {
                label pointI = splitPoints[i];

                // Check that all cells are marked
                const labelList& pCells = pointCells[pointI];

                bool isNotMarked = false;

                forAll(pCells, pCellI)
                {
                    if (!unrefineHexCell.get(pCells[pCellI]))
                    {
                        isNotMarked = true;
                        break;
                    }
                }

                if (!isNotMarked)
                {
                    splitPointCandidates.append(pointI);
                }
            }

            splitPointCandidates.shrink();

            pointsToUnrefine.transfer
            (
                splitPointCandidates
            );



            // TET: get list of tet cells candidateds to unrefine
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            label nUnrefineTetCells = 0;

            forAll(unrefineTetCell, i)
            {
                if(unrefineTetCell.get(i))
                {
                    nUnrefineTetCells++;
                }
            }

            tetCellsToUnrefine.resize(nUnrefineTetCells);
            label k = 0;

            forAll(unrefineTetCell, i)
            {
                if(unrefineTetCell.get(i))
                {
                    tetCellsToUnrefine[k++] = i;
                }
            }
        }

        // We now have a list of candidates:
        // 1. split points for hex cells
        // 2. tet cells
        // feed these into consistency checker

//        struct timeval timeStart, timeStop, diff;
//        label sec, usec;
//        gettimeofday(&timeStart, NULL);

        meshCutter_.consistentUnrefinement
        (
            pointsToUnrefine,
            tetCellsToUnrefine,
            unrefineHexSplitPoints,
            unrefineTetCellsParent,
            false
        );

//        gettimeofday(&timeStop, NULL);
//        timersub(&timeStop, &timeStart, &diff);
//        sec = diff.tv_sec;
//        usec = diff.tv_usec;
//        Info << "refTime consistentUnrefinement:    " << sec << "." << usec << nl;
    }

    label nSplitPoints = returnReduce
    (
        unrefineHexSplitPoints.size(),
        sumOp<label>()
    );

    label nTetCellsToUnrefine = returnReduce
    (
        unrefineTetCellsParent.size(),
        sumOp<label>()
    );

    Info << "Selected " << nSplitPoints << " split points out of a possible "
         << nTotalSplitPoints << ", and " << 8*nTetCellsToUnrefine
         << " tet cells for unrefinement." << endl;

    if (nTetCellsToUnrefine + nSplitPoints > 0)
    {
        // Unrefine/update mesh
        unrefine
        (
            unrefineHexSplitPoints,
            unrefineTetCellsParent
        );

        return true;
    }
    else
    {
        return false;
    }
}


//- Compute the cell average at the requested difference in level from
//  the input cell
double Foam::dynamicRefineMRAFvMesh::computeAverage
(
    label cellI,
    label refLevel
)
{
    // splitCell ID of the cell
    label mySCID
        = meshCutter_.history().visibleCells()[cellI];

    // splitCell ID of the parent
    label targetAncestorSCID
        = meshCutter_.history().splitCells()[mySCID].parent_;

    // if average needed one level higher, then get the grandparent splitCellID
    if (refLevel == 1)
    {
        targetAncestorSCID =
            meshCutter_.history().splitCells()[targetAncestorSCID].parent_;
    }

    // return the average value of the target ancestor
    return meshCutter_.history().splitCells()[targetAncestorSCID].cellAvg_;
}


void Foam::dynamicRefineMRAFvMesh::treatCell
(
    label cellI,
    const volScalarField &vFld,
    PackedBoolList &refineCell,
    PackedBoolList &unrefineCell
)
{
    if (isCellTreatable(cellI))
    {
        // refinement history
        const refinementMRAHistory &history
            = meshCutter_.history();

        // visible cells
        const labelList &visibleCells
            = history.visibleCells();

        // split cells
        const DynamicList<refinementMRAHistory::splitCell8> &splitCells
            = history.splitCells();

        // level of cellI
        const labelList &cellLevel
            = meshCutter_.cellLevel();

        // splitCellID of cellI
        label cellISCID
            = visibleCells[cellI];

        // splitCellID of cellI's parent
        label parentSCID
            = splitCells[cellISCID].parent_;

        // cells that were generated with cellI - note that these may have
        // been refined further
        FixedList<label, 8> cellIGenCells
            = splitCells[parentSCID].addedCellIDsPtr_();

        // cellLevel of cellI
        label levelOwn =
            cellLevel[cellI];

        // this list holds the ID of cellI's brothers. If some of them were
        // refined further, the list at that cell's position holds -1
        labelList brothers(8, -1);
        forAll (brothers, i)
        {
            if (cellLevel[cellIGenCells[i]] == levelOwn)
            {
                brothers[i] = cellIGenCells[i];
            }
        }

        //
        // do an MRA. We assume a CDF11 MRA - i.e., we get the coarse
        // cell value as the average of the 8 bros, and details as the
        // difference for each brother
        //

        // average of all cellI brothers (or its nephews)
        double averageVal = computeAverage(cellI, 0);

        // local threshold for MRA: Roussel et al 2003
        double localThreshold = 0.0;

        if (levelOwn < relativeAccuracyLevel_)
        {
            localThreshold =
                globalAccuracyThreshold_
              * pow
                (
                    0.5,
                    -3*(levelOwn - relativeAccuracyLevel_)
                );
        }
        else
        {
            localThreshold = globalAccuracyThreshold_;
        }

        forAll(brothers, i)
        {
            label bI = brothers[i];

            if (isCellTreated_.get(bI))
            {
                continue;
            }

            if (bI < 0)
            {
                // This means, that one of the brothers was refined further.
                // Thus the unrefinement of cellI and its brothers is not
                // possible

                unrefineCell.unset(cellI);
            }
            else
            {
                // this brother has not been treated yet (see below). We need to
                // check its detail to see if it needs refinement or not.
                label levelBI = cellLevel[bI];

                // detail is computed as the difference between the cell value
                // and the average one level higher.
                // TODO: should we use relative or absolute?
                double detailI =
                    fabs(vFld[bI] - averageVal) / normalizingParameter_;

                if (detailI >= localThreshold)
                {
                    // the detail is above threshold. This means that the cell
                    // definitely cannot be unrefined. Mark the flag as such.
                    // Also specify the the cell is MRADetailed

                    unrefineCell.unset(cellI);

                    unrefineCell.unset(bI);
                    isCellMRADetailed_.set(bI);
                    isCellTreated_.set(bI);

                    if (levelBI == maxRefinement_)
                    {
                        // there is no possibility of further refining this cell
                        // So unmark its refinement.

                        refineCell.unset(bI);
                    } // end refLevel if

                    // also mark all faceNeighbours of same level or less
                    // for supportive refinement. If it was marked for
                    // unrefinement, cancel the unrefinement

                    const labelList &cPoints = cellPoints()[bI];
                    forAll(cPoints, k)
                    {
                        label cP = cPoints[k];
                        const labelList &nbCells = pointCells()[cP];

                        forAll(nbCells, j)
                        {
                            label nbJ = nbCells[j];
                            label levelNbJ = cellLevel[nbJ];

                            if (levelNbJ < levelBI)
                            {
                                unrefineCell.unset(nbJ);
                                refineCell.set(nbJ);
                                isCellTreated_.set(nbJ);
                            }
                            else if (levelNbJ == levelBI)
                            {
                                unrefineCell.unset(nbJ);

                                if (refineCell.get(bI))
                                {
                                    refineCell.set(nbJ);
                                    isCellTreated_.set(nbJ);
                                }
                            }
                            else if (levelNbJ == levelBI+1)
                            {
                                if (refineCell.get(bI))
                                {
                                    unrefineCell.unset(nbJ);
                                }
                            }

                        } // end pointCells
                    } // end cellPoints

                } // end detail if

                else if (detailI >= unrefTriggerFrac_*localThreshold)
                {
                    // check if the brother's detail is not low enough for
                    // unrefinement. If not, mark off unrefinement. Also
                    // specify that the cell is MRADetailed

                    unrefineCell.unset(cellI);

                    refineCell.unset(bI);
                    unrefineCell.unset(bI);
                    isCellMRADetailed_.set(bI);
                    isCellTreated_.set(bI);

                    // also mark all faceNeighbours of less levels
                    // for supportive refinement. If it was marked for
                    // unrefinement, cancel the unrefinement

                    const labelList &cPoints = cellPoints()[bI];
                    forAll(cPoints, k)
                    {
                        label cP = cPoints[k];
                        const labelList &nbCells = pointCells()[cP];

                        forAll(nbCells, j)
                        {
                            label nbJ = nbCells[j];
                            label levelNbJ = cellLevel[nbJ];

                            if (levelNbJ < levelBI)
                            {
                                refineCell.set(nbJ);
                                unrefineCell.unset(nbJ);
                                isCellTreated_.set(nbJ);
                            }
                            else if (levelNbJ == levelBI)
                            {
                                unrefineCell.unset(nbJ);
                            }

                        } // end pointCells
                    } // end cellPoints
                }
                else
                {
                    refineCell.unset(bI);
                }
            } // end brother > -1 if
        } // end brother loop


        // check if some pointNeighbours need MRA-refinement

        const labelList &cPoints = cellPoints()[cellI];
        forAll(cPoints, k)
        {
            label cP = cPoints[k];
            const labelList &nbCells = pointCells()[cP];

            forAll(nbCells, j)
            {
                label nbJ = nbCells[j];
                label levelNbJ = cellLevel[nbJ];

                if (isCellMRADetailed_.get(nbJ))
                {
                    if (refineCell.get(nbJ))
                    {
                        if (levelNbJ >= levelOwn)
                        {
                            refineCell.set(cellI);
                            unrefineCell.unset(cellI);
                            isCellTreated_.set(cellI);
                            // can add a break here
                        }
                        else if (levelNbJ == levelOwn-1)
                        {
                            unrefineCell.unset(cellI);
                            isCellTreated_.set(cellI);
                        }
                    }
                    else
                    {
                        if (levelNbJ > levelOwn)
                        {
                            refineCell.set(cellI);
                            unrefineCell.unset(cellI);
                            isCellTreated_.set(cellI);
                            // can add a break here
                        }
                        else if (levelNbJ == levelOwn)
                        {
                            unrefineCell.unset(cellI);
                            isCellTreated_.set(cellI);
                        }
                    }
                }
            } // end pointCells
        } // end cellPoints


        // ---------------------------------------------------------------------
        // At this point, we have decided if the cell and its brothers need
        // refinement or not. In case the unrefinement of cellI and its brothers
        // is still possible, we need to make sure that they have a grandparent,
        // otherwise MRA will not be possible on them
        // ---------------------------------------------------------------------

        if (unrefineCell.get(cellI))
        {
            // check if it has a grandparent

            label parentOwn = history.parentIndex(cellI);
            label grandparentOwn = splitCells[parentOwn].parent_;

            if (grandparentOwn == -1)
            {
                // don't unrefine if the cell has no grandparent

                unrefineCell.unset(cellI);
                isCellTreated_.set(cellI, 1);
            }

            // If after all these checks the cell is capable of being unrefined
            // then do an MRA one cellLevel above it. If the detail in its
            // compartment is under the appropriate threshold, then mark cellI
            // and its brothers for unrefinement.
            // If not, cancel the unrefinement of same-level cells around them.
            // Also, mark any lower levelled cells around them for refinement.

            if (unrefineCell.get(cellI))
            {
                // compute the average value one level higher, the coarse level
                // threshold, and finally the detail.

                double coarseAverageVal = computeAverage(cellI, 1);

                double coarseThreshold = 0.0;
                if (levelOwn <= relativeAccuracyLevel_)
                {
                    coarseThreshold =
                        globalAccuracyThreshold_
                      * pow
                        (
                            0.5,
                            -3*(levelOwn-1 - relativeAccuracyLevel_)
                        );
                }
                else
                {
                    coarseThreshold = globalAccuracyThreshold_;
                }

                double coarseDetail =
                    fabs(averageVal - coarseAverageVal) / normalizingParameter_;

                forAll(brothers, i)
                {
                    label bI = brothers[i];

                    if (isCellTreated_.get(bI))
                    {
                        continue;
                    }

                    label levelBI = cellLevel[bI];

                    if (levelBI != levelOwn)
                    {
                        FatalErrorInFunction
                            << "At cell " << cellI << ". Brothers okay for "
                            << "unrefinement, but they have different levels. "
                            << "This is not supposed to happen."
                            << abort(FatalError);
                    }

                    if
                    (
                         (coarseDetail < unrefTriggerFrac_ * coarseThreshold
                      || levelBI > maxRefinement_)
                    )
                    {
                    }
                    else
                    {
                        unrefineCell.unset(cellI);

                        isCellMRADetailed_.set(bI);
                        unrefineCell.unset(bI);

                        const labelList &cPoints = cellPoints()[bI];
                        forAll(cPoints, k)
                        {
                            label cP = cPoints[k];
                            const labelList &nbCells = pointCells()[cP];

                            forAll(nbCells, j)
                            {
                                label nbJ = nbCells[j];
                                label levelNbJ = cellLevel[nbJ];

                                if (levelNbJ < levelBI)
                                {
                                    refineCell.set(nbJ);
                                    unrefineCell.unset(nbJ);
                                    isCellTreated_.set(nbJ);
                                }
                                else if (levelNbJ == levelBI)
                                {
                                    unrefineCell.unset(nbJ);
                                }
                            } // end pointCells
                        } // end cellPoints
                    } // end threshold & level if

                    isCellTreated_.set(bI);

                } // end brothers loop
            } // end isMyUnrefStillOkay
        } // end isMyUnrefOkay
        else
        {
            isCellTreated_.set(cellI, 1);
        }
    } // end isCellTreatable if
    else
    {
        refineCell.unset(cellI);
        unrefineCell.unset(cellI);
        isCellTreated_.set(cellI, 1);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicRefineMRAFvMesh::~dynamicRefineMRAFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicRefineMRAFvMesh::update()
{
    // Re-read dictionary. Choosen since usually -small so trivial amount
    // of time compared to actual refinement. Also very useful to be able
    // to modify on-the-fly.
    dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict(typeName + "Coeffs")
    );

    // -------------------------------------------------------------------------

    /*
     * Parameters to be specified by the user
     */

    // -------------------------------------------------------------------------


    // interval at which we should do mesh update for refinement
    label refineInterval =
        readLabel(refineDict.lookup("refineInterval"));

    // maximum number of cells allowed
    maxCells_ =
        readLabel(refineDict.lookup("maxCells"));

    // maximum number of refinement levels
    maxRefinement_ =
        readLabel(refineDict.lookup("maxRefinement"));

    // refinement-driving field
    const word fieldName(refineDict.lookup("field"));

    // fraction of the local threshold under which cells can be unrefined
    unrefTriggerFrac_ =
        readFloatScalar(refineDict.lookup("unrefTriggerFrac"));

    // cell level under and until which the accuracy level is scaled
    // based on Harten's formula
    relativeAccuracyLevel_ =
        readLabel(refineDict.lookup("accuracyLevel"));

    // Accuracy of the solution desired by the user
    globalAccuracyThreshold_ =
        readFloatScalar(refineDict.lookup("fverror"));


    bool hasChanged = false;
    bool globalHasChanged = false;
    bool hasChangedInitRef = false;
    bool hasChangedRefinement = false;
    bool hasChangedUnrefinement = false;

    if (refineInterval == 0)
    {
        topoChanging(hasChanged);

        return false;
    }
    else if (refineInterval < 0)
    {
        FatalErrorInFunction
            << "Illegal refineInterval " << refineInterval << nl
            << "The refineInterval setting in the dynamicMeshDict should"
            << " be >= 1." << nl
            << exit(FatalError);
    }

    // -------------------------------------------------------------------------

    /*
     * This is where we initialize the mesh.
     * This means refining all cells once.
     */

    // -------------------------------------------------------------------------

    if (!isInitialized_)
    {
        // Check if the mesh is initial-refined. If not, then refine each
        // cell so that it has a parent - this is necessary for the MRA.

        isInitialized_ =
            readBool(refineDict.lookup("isMeshInitialRefined"));

        if (!isInitialized_)
        {
            Info << "Ensuring each refinable cell has a father." << endl;
            PackedBoolList refineCell;

            for(label cellI = 0; cellI < nCells(); cellI++)
            {
                if
                (
                    meshCutter_.cellLevel()[cellI] >= 0
                 && meshCutter_.cellLevel()[cellI] < maxRefinement_
                )
                {
                    refineCell.set(cellI, 1);
                }
            }

            // conduct initial refinement

            {
                PackedBoolList dummy;
                hasChangedInitRef = conductRefinement(refineCell, dummy);
            }

            topoChanging(hasChangedInitRef);

            if (hasChangedInitRef)
            {
                // Reset moving flag (if any). If not using inflation we'll
                // not move, if are using inflation any follow on movePoints
                // will set it
                moving(false);
            }

            // call to initialize mesh for phase interface in multi-phase sims
            bool initSetFields = refineDict.lookupOrDefault<bool>("initSetFields", false);
            if (initSetFields)
            {
                initializeSetFields();
            }

            isInitialized_ = true;
        }
    }


    // Note: cannot refine at time 0 since no V0 present since mesh not
    //       moved yet.

    if (time().timeIndex() > 0 && time().timeIndex() % refineInterval == 0)
    {
        // reset the isCellTreated and isCellMRADetailed arrays

        isCellTreated_.reset();
        isCellTreated_.resize(nCells());

        isCellMRADetailed_.reset();
        isCellMRADetailed_.resize(nCells());

        if (maxCells_ <= 0)
        {
            FatalErrorInFunction
                << "Illegal maximum number of cells " << maxCells_ << nl
                << "The maxCells setting in the dynamicMeshDict should"
                << " be > 0." << nl
                << exit(FatalError);
        }

        if (maxRefinement_ <= 0)
        {
            FatalErrorInFunction
                << "Illegal maximum refinement level " << maxRefinement_ << nl
                << "The maxCells setting in the dynamicMeshDict should"
                << " be > 0." << nl
                << exit(FatalError);
        }

        // read the field based on which refinement/unrefinement is to be
        // done
        const volScalarField& vFld = lookupObject<volScalarField>(fieldName);

        // read the normalizing parameter as the maximum value of the field to
        // be refined
        bool normalize = refineDict.lookupOrDefault<bool>("normalize", false);

        if (normalize)
        {
            normalizingParameter_ = max(vFld).value();
        }
        else
        {
            normalizingParameter_ = 1.0;
        }

        // Cells marked for refinement or otherwise protected from unrefinement.
        PackedBoolList refineCell(nCells(), true);

        // Cells marked for unrefinement
        PackedBoolList unrefineCell(nCells(), true);

        forAll(refineCell, i)
        {
            refineCell.set(i);
            unrefineCell.set(i);

            if (meshCutter_.cellLevel()[i] == maxRefinement_)
            {
                refineCell.unset(i);
            }
        }

        // ---------------------------------------------------------------------

        /*
         * This is the actual procedure for adaptive mesh refinement. We first
         * select cells to be refined based on Multiresolution Analysis. Then,
         * we feed the refinement information into the hexTetRef8 engine, transform
         * the mesh, and correct for fluxes.
         */

        // ---------------------------------------------------------------------

//        struct timeval timeStart, timeStop, diff;
//        label sec, usec;
//        gettimeofday(&timeStart, NULL);

        // fill the history for the new field values
        meshCutter_.fillHistory(vFld);

        // treat all cells
        forAll(cells(), cellI)
        {
            if (!isCellTreated_[cellI])
            {
                treatCell(cellI, vFld, refineCell, unrefineCell);
            }
        }

        forAll(cells(), cellI)
        {
            if (!isCellTreated_[cellI])
            {
                FatalErrorInFunction
                     << "cell " << cellI << " not treated." << nl
                     << exit(FatalError);
            }
            if (refineCell.get(cellI) && unrefineCell.get(cellI))
            {
                FatalErrorInFunction
                     << "cell " << cellI << " needs both ref and unref." << nl
                     << exit(FatalError);
            }
        }

//        gettimeofday(&timeStop, NULL);
//        timersub(&timeStop, &timeStart, &diff);
//        sec = diff.tv_sec;
//        usec = diff.tv_usec;
//        Info << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << nl;
//        Info << "refTime MRA:                       " << sec << "." << usec << nl;

        // Do the refinement and unrefinement
        hasChangedRefinement = conductRefinement
        (
            refineCell,
            unrefineCell
        );

        hasChangedUnrefinement = conductUnrefinement
        (
            unrefineCell
        );

        if ((nRefinementIterations_ % 10) == 0)
        {
            // Compact refinement history occassionally (how often?).
            // Unrefinement causes holes in the refinementMRAHistory.
            const_cast<refinementMRAHistory&>(meshCutter().history()).compact();
        }
        nRefinementIterations_++;

        hasChanged = hasChangedRefinement
                  || hasChangedUnrefinement
                  || hasChangedInitRef;

        globalHasChanged = returnReduce(hasChanged, orOp<bool>());


        // Check if mesh needs load balancing.

        bool enableBalancing = readBool(refineDict.lookup("enableBalancing"));

        if
        (
            Pstream::parRun()       // if parallel run
         && globalHasChanged        // if the mesh has changed
         && enableBalancing
        )
        {
            const floatScalar allowableImbalance =
                readFloatScalar(refineDict.lookup("allowableImbalance"));

            balanceMesh(allowableImbalance);
        }
    }

    topoChanging(globalHasChanged);

    if (globalHasChanged)
    {
        // Reset moving flag (if any). If not using inflation we'll not move,
        // if are using inflation any follow on movePoints will set it.
        moving(false);
    }

    return globalHasChanged;
}


void Foam::dynamicRefineMRAFvMesh::initializeSetFields()
{
    Info << "Initializing phase.." << nl;

    // read the setFields dictionary
    dictionary setFieldsDict
    (
        IOdictionary
        (
            IOobject
            (
                "setFieldsDict",
                time().system(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        )
    );

    PtrList<entry> regions(setFieldsDict.lookup("regions"));

    bool globalHasChanged = true;

    while (globalHasChanged)
    {
        bool changing = false;

        List<volumeType> overallPointStatus(nPoints(), volumeType::OUTSIDE);

        forAll(regions, regionI)
        {
            const entry& region = regions[regionI];

            autoPtr<topoSetSource> source =
                topoSetSource::New(region.keyword(), *this, region.dict());

            List<volumeType> localPointStatus(nPoints(), volumeType::OUTSIDE);

            IOobject dummy
            (
                "abc",                         // dummy name
                time().constant(),             // directory
                "triSurface",                  // instance
                time()                         // registry
            );

            if (dynamic_cast<boxToCell *>(&source()))
            {
                boxToCell boxSource(*this, region.dict());
                treeBoundBoxList bbs = boxSource.getBoundingBoxList();

                forAll(bbs, i)
                {
                    treeBoundBox bb = bbs[i];
                    searchableBox sBox
                    (
                        dummy,
                        bb
                    );

                    List<volumeType> superLocalPointStatus(nPoints(), volumeType::OUTSIDE);
                    sBox.getVolumeType(points(), superLocalPointStatus);

                    forAll(superLocalPointStatus, i)
                    {
                        if (superLocalPointStatus[i] == volumeType::INSIDE)
                        {
                            localPointStatus[i] = volumeType::INSIDE;
                        }
                    }
                }
            }
            else if (dynamic_cast<sphereToCell *>(&source()))
            {
                sphereToCell sphereSource(*this, region.dict());
                point centre = sphereSource.getCentre();
                scalar radius = sphereSource.getRadius();

                searchableSphere sSphere
                (
                    dummy,
                    centre,
                    radius
                );

                sSphere.getVolumeType(points(), localPointStatus);
            }
            else if (dynamic_cast<cylinderToCell *>(&source()))
            {
                cylinderToCell cylinderSource(*this, region.dict());
                point p1 = cylinderSource.getPoint1();
                point p2 = cylinderSource.getPoint2();
                scalar radius = cylinderSource.getRadius();

                searchableCylinder sCylinder
                (
                    dummy,
                    p1,
                    p2,
                    radius
                );

                sCylinder.getVolumeType(points(), localPointStatus);
            }
            else
            {
                WarningInFunction
                    << "This topoSetSource is not yet supported."
                    << nl;
            }

            forAll(localPointStatus, i)
            {
                if (localPointStatus[i] == volumeType::INSIDE)
                {
                    overallPointStatus[i] = volumeType::INSIDE;
                }
            }
        }

        // Cells marked for refinement or otherwise protected from unrefinement.
        PackedBoolList refineCell(nCells(), false);

        // Cells marked for unrefinement - dummy
        PackedBoolList unrefineCell(0);

        forAll(cells(), i)
        {
            labelList cPoints = cellPoints()[i];

            forAll(cPoints, j)
            {
                label pJ = cPoints[j];
                label pNext = cPoints[cPoints.fcIndex(j)];

                if (
                    (overallPointStatus[pJ] != overallPointStatus[pNext])
                 && (meshCutter_.cellLevel()[i] < maxRefinement_)
                )
                {
                    refineCell.set(i);
                    changing = true;
                    break;
                }
            }
        }

        conductRefinement(refineCell, unrefineCell);
        globalHasChanged = returnReduce(changing, orOp<bool>());

        if
        (
            Pstream::parRun()       // if parallel run
         && globalHasChanged        // if the mesh has changed
        )
        {
//            balanceMesh(0.5);
        }
    }

    // initialize
    setFieldsDynamic(setFieldsDict, *this);
}


bool Foam::dynamicRefineMRAFvMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    // Force refinement data to go to the current time directory.
    const_cast<hexTetRef8&>(meshCutter_).setInstance(time().timeName());

    bool writeOk =
    (
        dynamicFvMesh::writeObjects(fmt, ver, cmp)
     && meshCutter_.write()
    );

    if (dumpLevel_)
    {
        volScalarField scalarCellLevel
        (
            IOobject
            (
                "cellLevel",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            *this,
            dimensionedScalar("level", dimless, 0)
        );

        const labelList& cellLevel = meshCutter_.cellLevel();

        forAll(cellLevel, cellI)
        {
            scalarCellLevel[cellI] = cellLevel[cellI];
        }

        writeOk = writeOk && scalarCellLevel.write();
    }

    if (dumpID_)
    {
        volScalarField scalarCellID
        (
            IOobject
            (
                "cellID",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            *this,
            dimensionedScalar("level", dimless, 0)
        );

        forAll(cells(), cellI)
        {
            scalarCellID[cellI] = cellI;
        }

        writeOk = writeOk && scalarCellID.write();
    }

    return writeOk;
}

// ************************************************************************* //
