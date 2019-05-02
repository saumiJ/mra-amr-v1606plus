/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "hexTetRef8.H"

#include "polyMesh.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "polyAddFace.H"
#include "polyAddPoint.H"
#include "polyAddCell.H"
#include "polyModifyFace.H"
#include "syncTools.H"
#include "faceSet.H"
#include "cellSet.H"
#include "pointSet.H"
#include "OFstream.H"
#include "Time.H"
#include "FaceCellWave.H"
#include "mapDistributePolyMesh.H"
#include "refinementData.H"
#include "refinementDistanceData.H"
#include "degenerateMatcher.H"
#include "face.H"
#include "SortableList.H"
#include "plane.H"
#include "vectorTools.H"

#include <set>
#include <algorithm>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hexTetRef8, 0);

    //- Reduction class. If x and y are not equal assign value.
    template<label value>
    class ifEqEqOp
    {
        public:
        void operator()(label& x, const label y) const
        {
            x = (x == y) ? x : value;
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::hexTetRef8::reorder
(
    const labelList& map,
    const label len,
    const label null,
    labelList& elems
)
{
    labelList newElems(len, null);

    forAll(elems, i)
    {
        label newI = map[i];

        if (newI >= len)
        {
            FatalErrorInFunction << abort(FatalError);
        }

        if (newI >= 0)
        {
            newElems[newI] = elems[i];
        }
    }

    elems.transfer(newElems);
}


void Foam::hexTetRef8::getFaceInfo
(
    const label faceI,
    label& patchID,
    label& zoneID,
    label& zoneFlip
) const
{
    patchID = -1;

    if (!mesh_.isInternalFace(faceI))
    {
        patchID = mesh_.boundaryMesh().whichPatch(faceI);
    }

    zoneID = mesh_.faceZones().whichZone(faceI);

    zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh_.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
    }
}

//void Foam::hexTetRef8::getNewFaceInfo
//(
//    polyTopoChange& meshMod,
//    const label faceI,
//    label& patchID,
//    label& zoneID,
//    label& zoneFlip
//) const
//{
//    patchID = -2;
//    zoneID  = -2;
//
//    label originalFace = meshMod.faceMap()[faceI];
//
//    patchID = mesh_.boundaryMesh().whichPatch(originalFace);
//
//    zoneID = mesh_.faceZones().whichZone(originalFace);
//
//    zoneFlip = false;
//
//    if (zoneID >= 0)
//    {
//        const faceZone& fZone = mesh_.faceZones()[zoneID];
//
//        zoneFlip = fZone.flipMap()[fZone.whichFace(originalFace)];
//    }
//
//    if (patchID == -2 || zoneID == -2)
//    {
//        FatalErrorInFunction
//            << "Unable to locate patchID/zoneID"
//            << abort(FatalError);
//    }
//}


// Adds a face on top of existing faceI.
Foam::label Foam::hexTetRef8::addFace
(
    polyTopoChange& meshMod,
    const label faceI,
    const face& newFace,
    const label own,
    const label nei
) const
{
    label patchID, zoneID, zoneFlip;

    getFaceInfo(faceI, patchID, zoneID, zoneFlip);

    label newFaceI = -1;

    if ((nei == -1) || (own < nei))
    {
        // Ordering ok.
        newFaceI = meshMod.setAction
        (
            polyAddFace
            (
                newFace,                    // face
                own,                        // owner
                nei,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                faceI,                      // master face for addition
                false,                      // flux flip
                patchID,                    // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face zone flip
            )
        );
    }
    else
    {
        // Reverse owner/neighbour
        newFaceI = meshMod.setAction
        (
            polyAddFace
            (
                newFace.reverseFace(),      // face
                nei,                        // owner
                own,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                faceI,                      // master face for addition
                false,                      // flux flip
                patchID,                    // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face zone flip
            )
        );
    }
    return newFaceI;
}


// Adds an internal face from an edge. Assumes orientation correct.
// Problem is that the face is between four new vertices. So what do we provide
// as master? The only existing mesh item we have is the edge we have split.
// Have to be careful in only using it if it has internal faces since otherwise
// polyMeshMorph will complain (because it cannot generate a sensible mapping
// for the face)
Foam::label Foam::hexTetRef8::addInternalFace
(
    polyTopoChange& meshMod,
    const label meshFaceI,
    const label meshPointI,
    const face& newFace,
    const label own,
    const label nei
) const
{
    if (mesh_.isInternalFace(meshFaceI))
    {
        return meshMod.setAction
        (
            polyAddFace
            (
                newFace,                    // face
                own,                        // owner
                nei,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                meshFaceI,                  // master face for addition
                false,                      // flux flip
                -1,                         // patch for face
                -1,                         // zone for face
                false                       // face zone flip
            )
        );
    }
    else
    {
        // Two choices:
        // - append (i.e. create out of nothing - will not be mapped)
        //   problem: field does not get mapped.
        // - inflate from point.
        //   problem: does interpolative mapping which constructs full
        //   volPointInterpolation!

        // For now create out of nothing

        return meshMod.setAction
        (
            polyAddFace
            (
                newFace,                    // face
                own,                        // owner
                nei,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                -1,                         // master face for addition
                false,                      // flux flip
                -1,                         // patch for face
                -1,                         // zone for face
                false                       // face zone flip
            )
        );


        ////- Inflate-from-point:
        //// Check if point has any internal faces we can use.
        //label masterPointI = -1;
        //
        //const labelList& pFaces = mesh_.pointFaces()[meshPointI];
        //
        //forAll(pFaces, i)
        //{
        //    if (mesh_.isInternalFace(pFaces[i]))
        //    {
        //        // meshPoint uses internal faces so ok to inflate from it
        //        masterPointI = meshPointI;
        //
        //        break;
        //    }
        //}
        //
        //return meshMod.setAction
        //(
        //    polyAddFace
        //    (
        //        newFace,                    // face
        //        own,                        // owner
        //        nei,                        // neighbour
        //        masterPointI,               // master point
        //        -1,                         // master edge
        //        -1,                         // master face for addition
        //        false,                      // flux flip
        //        -1,                         // patch for face
        //        -1,                         // zone for face
        //        false                       // face zone flip
        //    )
        //);
    }
}


// Modifies existing faceI for either new owner/neighbour or new face points.
void Foam::hexTetRef8::modFace
(
    polyTopoChange& meshMod,
    const label faceI,
    const face& newFace,
    const label own,
    const label nei
) const
{
    label patchID, zoneID, zoneFlip;

    getFaceInfo(faceI, patchID, zoneID, zoneFlip);

    if
    (
        (own != mesh_.faceOwner()[faceI])
     || (
            mesh_.isInternalFace(faceI)
         && (nei != mesh_.faceNeighbour()[faceI])
        )
     || (newFace != mesh_.faces()[faceI])
    )
    {
        if ((nei == -1) || (own < nei))
        {
            meshMod.setAction
            (
                polyModifyFace
                (
                    newFace,            // modified face
                    faceI,              // label of face being modified
                    own,                // owner
                    nei,                // neighbour
                    false,              // face flip
                    patchID,            // patch for face
                    false,              // remove from zone
                    zoneID,             // zone for face
                    zoneFlip            // face flip in zone
                )
            );
        }
        else
        {
            meshMod.setAction
            (
                polyModifyFace
                (
                    newFace.reverseFace(),  // modified face
                    faceI,                  // label of face being modified
                    nei,                    // owner
                    own,                    // neighbour
                    false,                  // face flip
                    patchID,                // patch for face
                    false,                  // remove from zone
                    zoneID,                 // zone for face
                    zoneFlip                // face flip in zone
                )
            );
        }
    }
}

// Modifies new faceI for either new owner/neighbour or new face points.
void Foam::hexTetRef8::modNewFace
(
    polyTopoChange& meshMod,
    const label faceI,
    const label parentFace,
    const face& newFace,
    const label own,
    const label nei
) const
{
    label patchID, zoneID, zoneFlip;

//    getNewFaceInfo(meshMod, faceI, patchID, zoneID, zoneFlip);
    getFaceInfo(parentFace, patchID, zoneID, zoneFlip);

    {
        if ((nei == -1) || (own < nei))
        {
            meshMod.setAction
            (
                polyModifyFace
                (
                    newFace,            // modified face
                    faceI,              // label of face being modified
                    own,                // owner
                    nei,                // neighbour
                    false,              // face flip
                    patchID,            // patch for face
                    false,              // remove from zone
                    zoneID,             // zone for face
                    zoneFlip            // face flip in zone
                )
            );
        }
        else
        {
            meshMod.setAction
            (
                polyModifyFace
                (
                    newFace.reverseFace(),  // modified face
                    faceI,                  // label of face being modified
                    nei,                    // owner
                    own,                    // neighbour
                    false,                  // face flip
                    patchID,                // patch for face
                    false,                  // remove from zone
                    zoneID,                 // zone for face
                    zoneFlip                // face flip in zone
                )
            );
        }
    }
}

// Bit complex way to determine the unrefined edge length.
Foam::scalar Foam::hexTetRef8::getLevel0EdgeLength() const
{
    if (cellLevel_.size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size()
            << endl
            << "This might be because of a restart with inconsistent cellLevel."
            << abort(FatalError);
    }

    // Determine minimum edge length per refinement level
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const scalar GREAT2 = sqr(GREAT);

    label nLevels = gMax(cellLevel_)+1;

    scalarField typEdgeLenSqr(nLevels, GREAT2);


    // 1. Look only at edges surrounded by cellLevel cells only.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // Per edge the cellLevel of connected cells. -1 if not set,
        // labelMax if different levels, otherwise levels of connected cells.
        labelList edgeLevel(mesh_.nEdges(), -1);

        forAll(cellLevel_, cellI)
        {
            const label cLevel = cellLevel_[cellI];

            const labelList& cEdges = mesh_.cellEdges(cellI);

            forAll(cEdges, i)
            {
                label edgeI = cEdges[i];

                if (edgeLevel[edgeI] == -1)
                {
                    edgeLevel[edgeI] = cLevel;
                }
                else if (edgeLevel[edgeI] == labelMax)
                {
                    // Already marked as on different cellLevels
                }
                else if (edgeLevel[edgeI] != cLevel)
                {
                    edgeLevel[edgeI] = labelMax;
                }
            }
        }

        // Make sure that edges with different levels on different processors
        // are also marked. Do the same test (edgeLevel != cLevel) on coupled
        // edges.
        syncTools::syncEdgeList
        (
            mesh_,
            edgeLevel,
            ifEqEqOp<labelMax>(),
            labelMin
        );

        // Now use the edgeLevel with a valid value to determine the
        // length per level.
        forAll(edgeLevel, edgeI)
        {
            const label eLevel = edgeLevel[edgeI];

            if (eLevel >= 0 && eLevel < labelMax)
            {
                const edge& e = mesh_.edges()[edgeI];

                scalar edgeLenSqr = magSqr(e.vec(mesh_.points()));

                typEdgeLenSqr[eLevel] = min(typEdgeLenSqr[eLevel], edgeLenSqr);
            }
        }
    }

    // Get the minimum per level over all processors. Note minimum so if
    // cells are not cubic we use the smallest edge side.
    Pstream::listCombineGather(typEdgeLenSqr, minEqOp<scalar>());
    Pstream::listCombineScatter(typEdgeLenSqr);

    if (debug)
    {
        Pout<< "hexTetRef8::getLevel0EdgeLength() :"
            << " After phase1: Edgelengths (squared) per refinementlevel:"
            << typEdgeLenSqr << endl;
    }


    // 2. For any levels where we haven't determined a valid length yet
    //    use any surrounding cell level. Here we use the max so we don't
    //    pick up levels between celllevel and higher celllevel (will have
    //    edges sized according to highest celllevel)
    //    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    scalarField maxEdgeLenSqr(nLevels, -GREAT2);

    forAll(cellLevel_, cellI)
    {
        const label cLevel = cellLevel_[cellI];

        const labelList& cEdges = mesh_.cellEdges(cellI);

        forAll(cEdges, i)
        {
            const edge& e = mesh_.edges()[cEdges[i]];

            scalar edgeLenSqr = magSqr(e.vec(mesh_.points()));

            maxEdgeLenSqr[cLevel] = max(maxEdgeLenSqr[cLevel], edgeLenSqr);
        }
    }

    Pstream::listCombineGather(maxEdgeLenSqr, maxEqOp<scalar>());
    Pstream::listCombineScatter(maxEdgeLenSqr);

    if (debug)
    {
        Pout<< "hexTetRef8::getLevel0EdgeLength() :"
            << " Crappy Edgelengths (squared) per refinementlevel:"
            << maxEdgeLenSqr << endl;
    }


    // 3. Combine the two sets of lengths
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(typEdgeLenSqr, levelI)
    {
        if (typEdgeLenSqr[levelI] == GREAT2 && maxEdgeLenSqr[levelI] >= 0)
        {
            typEdgeLenSqr[levelI] = maxEdgeLenSqr[levelI];
        }
    }

    if (debug)
    {
        Pout<< "hexTetRef8::getLevel0EdgeLength() :"
            << " Final Edgelengths (squared) per refinementlevel:"
            << typEdgeLenSqr << endl;
    }

    // Find lowest level present
    scalar level0Size = -1;

    forAll(typEdgeLenSqr, levelI)
    {
        scalar lenSqr = typEdgeLenSqr[levelI];

        if (lenSqr < GREAT2)
        {
            level0Size = Foam::sqrt(lenSqr)*(1<<levelI);

            if (debug)
            {
                Pout<< "hexTetRef8::getLevel0EdgeLength() :"
                    << " For level:" << levelI
                    << " have edgeLen:" << Foam::sqrt(lenSqr)
                    << " with equivalent level0 len:" << level0Size
                    << endl;
            }
            break;
        }
    }

    if (level0Size == -1)
    {
        FatalErrorInFunction
            << "Problem : typEdgeLenSqr:" << typEdgeLenSqr << abort(FatalError);
    }

    return level0Size;
}


// Check whether pointI is an anchor on cellI.
// If it is not check whether any other point on the face is an anchor cell.
Foam::label Foam::hexTetRef8::getAnchorCell
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const label cellI,
    const label faceI,
    const label pointI
) const
{
    if (cellAnchorPoints[cellI].size())
    {
        label index = findIndex(cellAnchorPoints[cellI], pointI);

        if (index != -1)
        {
            return cellAddedCells[cellI][index];
        }


        // pointI is not an anchor cell.
        // Maybe we are already a refined face so check all the face
        // vertices.
        const face& f = mesh_.faces()[faceI];

        forAll(f, fp)
        {
            label index = findIndex(cellAnchorPoints[cellI], f[fp]);

            if (index != -1)
            {
                return cellAddedCells[cellI][index];
            }
        }

        // Problem.
        dumpCell(cellI);
        Perr<< "cell:" << cellI << " anchorPoints:" << cellAnchorPoints[cellI]
            << endl;

        FatalErrorInFunction
            << "Could not find point " << pointI
            << " in the anchorPoints for cell " << cellI << endl
            << "Does your original mesh obey the 2:1 constraint and"
            << " did you use consistentRefinement to make your cells to refine"
            << " obey this constraint as well?"
            << abort(FatalError);

        return -1;
    }
    else
    {
        return cellI;
    }
}


// Get new owner and neighbour
void Foam::hexTetRef8::getFaceNeighbours
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const label faceI,
    const label pointI,

    label& own,
    label& nei
) const
{
    // Is owner split?
    own = getAnchorCell
    (
        cellAnchorPoints,
        cellAddedCells,
        mesh_.faceOwner()[faceI],
        faceI,
        pointI
    );

    if (mesh_.isInternalFace(faceI))
    {
        nei = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            mesh_.faceNeighbour()[faceI],
            faceI,
            pointI
        );
    }
    else
    {
        nei = -1;
    }
}


// Get point with the lowest pointLevel
Foam::label Foam::hexTetRef8::findMinLevel(const labelList& f) const
{
    label minLevel = labelMax;
    label minFp = -1;

    forAll(f, fp)
    {
        label level = pointLevel_[f[fp]];

        if (level < minLevel)
        {
            minLevel = level;
            minFp = fp;
        }
    }

    return minFp;
}


// Get point with the highest pointLevel
Foam::label Foam::hexTetRef8::findMaxLevel(const labelList& f) const
{
    label maxLevel = labelMin;
    label maxFp = -1;

    forAll(f, fp)
    {
        label level = pointLevel_[f[fp]];

        if (level > maxLevel)
        {
            maxLevel = level;
            maxFp = fp;
        }
    }

    return maxFp;
}


Foam::label Foam::hexTetRef8::countAnchors
(
    const labelList& f,
    const label anchorLevel
) const
{
    label nAnchors = 0;

    forAll(f, fp)
    {
        if (pointLevel_[f[fp]] <= anchorLevel)
        {
            nAnchors++;
        }
    }
    return nAnchors;
}


Foam::label Foam::hexTetRef8::countAnchors
(
    const labelList& f,
    const label anchorLevel,
    const DynamicList<label> &newPointLevel
) const
{
    label nAnchors = 0;

    forAll(f, fp)
    {
        if (newPointLevel[f[fp]] <= anchorLevel)
        {
            nAnchors++;
        }
    }
    return nAnchors;
}


void Foam::hexTetRef8::dumpCell(const label cellI) const
{
    OFstream str(mesh_.time().path()/"cell_" + Foam::name(cellI) + ".obj");
    Pout<< "hexTetRef8 : Dumping cell as obj to " << str.name() << endl;

    const cell& cFaces = mesh_.cells()[cellI];

    Map<label> pointToObjVert;
    label objVertI = 0;

    forAll(cFaces, i)
    {
        const face& f = mesh_.faces()[cFaces[i]];

        forAll(f, fp)
        {
            if (pointToObjVert.insert(f[fp], objVertI))
            {
                meshTools::writeOBJ(str, mesh_.points()[f[fp]]);
                objVertI++;
            }
        }
    }

    forAll(cFaces, i)
    {
        const face& f = mesh_.faces()[cFaces[i]];

        forAll(f, fp)
        {
            label pointI = f[fp];
            label nexPointI = f[f.fcIndex(fp)];

            str << "l " << pointToObjVert[pointI]+1
                << ' ' << pointToObjVert[nexPointI]+1 << nl;
        }
    }
}


// Find point with certain pointLevel. Skip any higher levels.
Foam::label Foam::hexTetRef8::findLevel
(
    const label faceI,
    const face& f,
    const label startFp,
    const bool searchForward,
    const label wantedLevel
) const
{
    label fp = startFp;

    forAll(f, i)
    {
        label pointI = f[fp];

        if (pointLevel_[pointI] < wantedLevel)
        {
            dumpCell(mesh_.faceOwner()[faceI]);
            if (mesh_.isInternalFace(faceI))
            {
                dumpCell(mesh_.faceNeighbour()[faceI]);
            }

            FatalErrorInFunction
                << "face:" << f
                << " level:" << UIndirectList<label>(pointLevel_, f)()
                << " startFp:" << startFp
                << " wantedLevel:" << wantedLevel
                << abort(FatalError);
        }
        else if (pointLevel_[pointI] == wantedLevel)
        {
            return fp;
        }

        if (searchForward)
        {
            fp = f.fcIndex(fp);
        }
        else
        {
            fp = f.rcIndex(fp);
        }
    }

    dumpCell(mesh_.faceOwner()[faceI]);
    if (mesh_.isInternalFace(faceI))
    {
        dumpCell(mesh_.faceNeighbour()[faceI]);
    }

    FatalErrorInFunction
        << "face:" << f
        << " level:" << UIndirectList<label>(pointLevel_, f)()
        << " startFp:" << startFp
        << " wantedLevel:" << wantedLevel
        << abort(FatalError);

    return -1;
}


//- Fill cellAvg of all history cells
void Foam::hexTetRef8::fillHistory
(
    const volScalarField& vFld
)
{
    history_.fillCellAvgs(vFld);
}


//- For the triangulated cell surface, find the faceID
//  of the face that shares edge e with face ownFace
//  and return the cellID of that face's cell
Foam::label Foam::hexTetRef8::edgeFaceNeighbour
(
    Foam::edge e,
    Foam::label ownFace,
    const Foam::DynamicList<Foam::DynamicList<Foam::label> > &triListList,
    const Foam::DynamicList<Foam::DynamicList<Foam::label> > &newCellListList,
    const Foam::polyTopoChange &meshMod
)
{
    forAll(triListList, i)
    {
        const DynamicList<label> &triListI = triListList[i];

        forAll(triListI, j)
        {
            label triJID = triListI[j];

            if (triJID != ownFace)
            {
                face triJ = meshMod.faces()[triJID];

                edgeList triJEdges = triJ.edges();

                forAll(triJEdges, k)
                {
                    edge eK = triJEdges[k];

                    if (eK.compare(eK, e) != 0)
                    {
                        return newCellListList[i][j];
                    }
                }
            }
        }
    }

    return -1;
}


void Foam::hexTetRef8::triangulateCellFaces
(
    label                            parentCell,
    PackedBoolList                   &isFaceSplit,
    DynamicList<DynamicList<label> > &triListListLocal,
    DynamicList<DynamicList<label> > &triListListGlobal,
    labelList                        &indexToSplitFaces,
    DynamicList<DynamicList<label> > &newCellListListLocal,
    DynamicList<label>               &newCellLevel,
    polyTopoChange                   &meshMod
)
{
    // get cell faces
    const cell &cellIFaces = mesh_.cells()[parentCell];

    // for all faces of the cell
    forAll(cellIFaces, j)
    {
        label parentFaceLabel = cellIFaces[j];

        // lists for new sub-triangles and sub-cells
        DynamicList<label> triList;
        DynamicList<label> newCellList;

        face  parentFace    = mesh_.faces()[parentFaceLabel];
        label ownParentFace = mesh_.faceOwner()[parentFaceLabel];
        label neiParentFace = -1;
        if (mesh_.isInternalFace(parentFaceLabel))
        {
            neiParentFace   = mesh_.faceNeighbour()[parentFaceLabel];
        }

        if (!isFaceSplit[parentFaceLabel])
        {
            // if face is not split, then compute a triangulation of the face
            // using the functionality provided in the face class

            label nSubFaces = parentFace.nTriangles(mesh_.points());
            faceList subFaces(nSubFaces);
            label subFacesI = 0;    // not used

            parentFace.triangles(mesh_.points(), subFacesI, subFaces);

            if (subFacesI != nSubFaces)
            {
                FatalErrorInFunction
                    << "Number of sub-triangles " << subFacesI
                    << " does not match estimate " << nSubFaces
                    << abort(FatalError);
            }

            // for each subFace of faceJ, we create a new triangular face on top
            // of it, except for the last triangular face over faceJ. This face
            // is recycled as the value of the old faceJ.
            // for every faceJ, its subFaces get a new cell, except the last
            // faceJ. For this last faceJ, all its subFaces get the recycled old
            // cell.

            forAll(subFaces, k)
            {
                face &subFace      = subFaces[k];  // kth sub-face
                label faceLabel    = -1;           // subFacek label
                label newCellLabel = -1;           // subFacek tet label
                label ownNew       = -1;           // subFacek owner
                label neiNew       = -2;           // subFacek neighbour

                if (k < nSubFaces - 1)
                {
                    newCellLabel =
                        meshMod.setAction
                        (
                            polyAddCell
                            (
                                -1,           // master point
                                -1,           // master edge
                                -1,           // master face
                                parentCell,   // master cell
                                mesh_.cellZones().whichZone(parentCell)
                            )
                        );

                    newCellLevel(newCellLabel) = cellLevel_[parentCell];

                    if (ownParentFace == parentCell)
                    {
                        // if parent cell was the owner of parent face

                        if (neiParentFace == -1)
                        {
                            neiNew = neiParentFace;
                            ownNew = newCellLabel;
                        }
                        else
                        {
                            neiNew = newCellLabel;
                            ownNew = neiParentFace;
                            subFace = subFace.reverseFace();
                        }
                    }
                    else
                    {
                        neiNew = newCellLabel;
                        ownNew = ownParentFace;
                    }

                    faceLabel = addFace
                    (
                        meshMod,
                        parentFaceLabel,
                        subFace,
                        ownNew,
                        neiNew
                    );
                }
                else if (k == nSubFaces - 1)
                {
                    // for the last sub-face of parent face

                    if (j == cellIFaces.size() - 1)
                    {
                        // if parent face was the last face of parent cell,
                        // recycle the parent cell as a tet.
                        newCellLabel = parentCell;
                        ownNew = ownParentFace;
                        neiNew = neiParentFace;
                    }
                    else
                    {
                        newCellLabel =
                            meshMod.setAction
                            (
                                polyAddCell
                                (
                                    -1,           // master point
                                    -1,           // master edge
                                    -1,           // master face
                                    parentCell,   // master cell
                                    mesh_.cellZones().whichZone(parentCell)
                                )
                            );

                        newCellLevel(newCellLabel) = cellLevel_[parentCell];

                        if (ownParentFace == parentCell)
                        {
                            if (neiParentFace == -1)
                            {
                                neiNew = neiParentFace;
                                ownNew = newCellLabel;
                            }
                            else
                            {
                                neiNew = newCellLabel;
                                ownNew = neiParentFace;
                                subFace = subFace.reverseFace();
                            }
                        }
                        else
                        {
                            neiNew = newCellLabel;
                            ownNew = ownParentFace;
                        }
                    }

                    faceLabel = parentFaceLabel;

                    modFace
                    (
                        meshMod,
                        parentFaceLabel,
                        subFace,
                        ownNew,
                        neiNew
                    );
                } // end which subFace if

                triList.append(faceLabel);          // add sub-face to list
                newCellList.append(newCellLabel);   // add sub-cell to list
            }

            // store index to sub-faces of parent face
            indexToSplitFaces[parentFaceLabel] = triListListGlobal.size();

            // store sub-face list
            triListListGlobal.append(triList);

            // set parent face as split
            isFaceSplit.set(parentFaceLabel, 1);
        }
        else    // face already split
        {
            // get labels of sub-faces for parent face
            triList = triListListGlobal [ indexToSplitFaces[parentFaceLabel] ];

            forAll(triList, k)
            {
                // label of kth sub-face
                label subFaceLabel = triList[k];

                // kth sub-face
                face subFace       = meshMod.faces()[subFaceLabel];

                // subFace previous owner
                label ownOld       = meshMod.faceOwner()[subFaceLabel];

                // subFace previous neighbour
                label neiOld       = meshMod.faceNeighbour()[subFaceLabel];

                label newCellLabel = -1;    // subFace tet label
                label ownNew = -1;          // subFace updated owner
                label neiNew = -1;          // subFace updated neighbour

                if
                (
                    (k == triList.size() - 1)
                 && (j == cellIFaces.size() - 1)
                )
                {
                    // if parent face is last face of parent cell
                    // AND subFace is last subFace of parent face
                    // recycle parent cell

                    newCellLabel = parentCell;
                    ownNew = ownOld;
                    neiNew = neiOld;
                }
                else
                {
                    newCellLabel =
                        meshMod.setAction
                        (
                            polyAddCell
                            (
                                -1,           // master point
                                -1,           // master edge
                                -1,           // master face
                                parentCell,   // master cell
                                mesh_.cellZones().whichZone(parentCell)
                            )
                        );  // add new tet

                    newCellLevel(newCellLabel) = cellLevel_[parentCell];

                    if (ownOld == parentCell)
                    {
                        if (neiOld != -1)
                        {
                            ownNew = neiOld;
                            neiNew = newCellLabel;
                            subFace = subFace.reverseFace();
                        }
                        else
                        {
                            ownNew = newCellLabel;
                            neiNew = neiOld;
                        }
                    }
                    else
                    {
                        ownNew = ownOld;
                        neiNew = newCellLabel;
                    }
                }

                modNewFace
                (
                    meshMod,
                    subFaceLabel,
                    parentFaceLabel,
                    subFace,
                    ownNew,
                    neiNew
                );

                newCellList.append(newCellLabel);
            }
        }

        triListListLocal.append(triList);
        newCellListListLocal.append(newCellList);
    }

    if (triListListLocal.size() != cellIFaces.size())
    {
         FatalErrorInFunction
            << "Faces the underwent decomposition not equal to original faces."
            << abort(FatalError);
    }
}


void Foam::hexTetRef8::addInternalFacesForTets
(
    label                            cellI,
    label                            cellICenter,
    DynamicList<DynamicList<label> > &triListListLocal,
    DynamicList<DynamicList<label> > &newCellListListLocal,
    DynamicList<DynamicList<label> > &triListListGlobal,
    labelList                        &indexToSplitFaces,
    polyTopoChange                   &meshMod
)
{
    // edges from points on faces to cellCenters
    DynamicList<edge> edgesConnectedToCenter;

    const cell &cellIFaces = mesh_.cells()[cellI];

    forAll(cellIFaces, j)
    {
        label subFacesIndex = indexToSplitFaces[cellIFaces[j]];
        DynamicList<label> subFaces = triListListGlobal[subFacesIndex];

        // for all subfaces (triangles) of faceJ
        forAll(subFaces, k)
        {
            face subFaceK = meshMod.faces()[subFaces[k]];
            edgeList kEdges = subFaceK.edges();

            // for all edges in subFaceK
            forAll(kEdges, l)
            {
                edge e = kEdges[l];

                bool isFacedWithCenter = false;

                // check if a face has already been created for this face
                forAll(edgesConnectedToCenter, m)
                {
                    if (e.compare(e, edgesConnectedToCenter[m]) != 0)
                    {
                        isFacedWithCenter = true;
                        break;
                    }
                }

                // if a face has not been created for this edge
                if (!isFacedWithCenter)
                {
                    edgesConnectedToCenter.append(e);

                    label nei = edgeFaceNeighbour
                        (
                            e,
                            subFaces[k],
                            triListListLocal,
                            newCellListListLocal,
                            meshMod
                        );  // get neighbour cell for new face

                    if (nei == -1)
                    {
                        FatalErrorInFunction
                            << "cell: " << cellI
                            << " triListList: " << triListListLocal
                            << " newCellListList: " << nl
                            << "No edge-sharing sub-cell found. "
                            << "Check if all cells have been added correctly."
                            << abort(FatalError);
                    }

                    label own = newCellListListLocal[j][k];
                    if (own > nei)
                    {
                        label temp = own;
                        own = nei;
                        nei = temp;
                    }

                    labelList newInternalFaceVertices(3);
                    if
                    (
                        meshMod.faceOwner()[subFaces[k]] == own
                     || (
                            newCellListListLocal[j][k] != own
                         && newCellListListLocal[j][k] !=
                                meshMod.faceOwner()[subFaces[k]]
                        )
                    )
                    {
                        // this ordering ensures correct orientation
                        newInternalFaceVertices[0] = e[0];
                        newInternalFaceVertices[1] = cellICenter;
                        newInternalFaceVertices[2] = e[1];
                    }
                    else
                    {
                        // this ordering ensures correct orientation
                        newInternalFaceVertices[0] = e[1];
                        newInternalFaceVertices[1] = cellICenter;
                        newInternalFaceVertices[2] = e[0];
                    }

                    face newInternalFace(newInternalFaceVertices);
                    label newInternalFaceLabel = meshMod.setAction
                    (
                        polyAddFace
                        (
                            newInternalFace,  // face
                            own,              // owner
                            nei,              // neighbour
                            e[0],             // master point
                            -1,               // master edge
                            subFaces[k],      // master face for addition
                            false,            // flux flip
                            -1,               // patch for face
                            -1,               // zone for face
                            false             // face zone flip
                        )
                    );

                    if (newInternalFaceLabel < 0)
                    {
                        FatalErrorInFunction
                            << "Failed to add face."
                            << abort(FatalError);
                    }
                }
            }
        }
    }
}

//- Break up polyhedra into tetrahedra
void Foam::hexTetRef8::decomposePolyhedra
(
    polyTopoChange &meshMod,
    PackedBoolList &protectedCells
)
{
    DynamicList<label> newCellLevel(cellLevel_.size());
    forAll(cellLevel_, cellI)
    {
        newCellLevel.append(cellLevel_[cellI]);
    }
    DynamicList<label> newPointLevel(pointLevel_.size());
    forAll(pointLevel_, pointI)
    {
        newPointLevel.append(pointLevel_[pointI]);
    }

    const cellList &cells = mesh_.cells();
    PackedBoolList isFaceSplit(mesh_.nFaces());
    labelList indexToSplitFaces(mesh_.nFaces(), -1);
    DynamicList<DynamicList<label> > triListListGlobal;

    forAll(cells, cellI)
    {
        if (protectedCells[cellI])    // if it is a non-hex cell
        {
            /*
             * Add cell mid point
             */

            label cellMidPoint = meshMod.setAction
            (
                polyAddPoint
                (
                    mesh_.cellCentres()[cellI],     // point
                    -1,                             // master point
                    -1,                             // zone for point
                    true                            // supports a cell
                )
            );
            newPointLevel(cellMidPoint) = cellLevel_[cellI];

            /*
             * Split cell faces into tetrahedra
             */

            DynamicList<DynamicList<label> > triListListLocal;
            DynamicList<DynamicList<label> > newCellListListLocal;

            triangulateCellFaces
            (
                cellI,
                isFaceSplit,
                triListListLocal,
                triListListGlobal,
                indexToSplitFaces,
                newCellListListLocal,
                newCellLevel,
                meshMod
            );

            /*
             * Add new internal faces
             */

            addInternalFacesForTets
            (
                cellI,
                cellMidPoint,
                triListListLocal,
                newCellListListLocal,
                triListListGlobal,
                indexToSplitFaces,
                meshMod
            );

        }
    }

    pointLevel_.transfer(newPointLevel);
    cellLevel_.transfer(newCellLevel);

    // Mark files as changed
    setInstance(mesh_.facesInstance());


    // Update the live split cells tree.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // New unrefinement structure
    if (history_.active())
    {
        if (debug)
        {
            Pout<< "hexTetRef8::setRefinement :"
                << " Updating refinement history to " << cellLevel_.size()
                << " cells" << endl;
        }

        // Extend refinement history for new cells
        history_.resize(cellLevel_.size());
    }
}


// Gets cell level such that the face has four points <= level.
Foam::label Foam::hexTetRef8::faceLevel
(
    const label faceI
) const
{
    const face& f = mesh_.faces()[faceI];

    label ownLevel = cellLevel_[mesh_.faceOwner()[faceI]];

    if (countAnchors(f, ownLevel) >= 3)
    {
        return ownLevel;
    }
    else if (countAnchors(f, ownLevel+1) >= 3)
    {
        return ownLevel+1;
    }
    else
    {
        return -1;
    }
}


void Foam::hexTetRef8::checkInternalOrientation
(
    polyTopoChange& meshMod,
    const label cellI,
    const label faceI,
    const point& ownPt,
    const point& neiPt,
    const face& newFace
)
{
    face compactFace(identity(newFace.size()));
    pointField compactPoints(meshMod.points(), newFace);

    vector n(compactFace.normal(compactPoints));

    vector dir(neiPt - ownPt);

    if ((dir & n) < 0)
    {
        FatalErrorInFunction
            << "cell:" << cellI << " old face:" << faceI
            << " newFace:" << newFace << endl
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " neiPt:" << neiPt
            << abort(FatalError);
    }

    vector fcToOwn(compactFace.centre(compactPoints) - ownPt);

    scalar s = (fcToOwn&n) / (dir&n);

    if (s < 0.1 || s > 0.9)
    {
        FatalErrorInFunction
            << "cell:" << cellI << " old face:" << faceI
            << " newFace:" << newFace << endl
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " neiPt:" << neiPt
            << " s:" << s
            << abort(FatalError);
    }
}


void Foam::hexTetRef8::checkBoundaryOrientation
(
    polyTopoChange& meshMod,
    const label cellI,
    const label faceI,
    const point& ownPt,
    const point& boundaryPt,
    const face& newFace
)
{
    face compactFace(identity(newFace.size()));
    pointField compactPoints(meshMod.points(), newFace);

    vector n(compactFace.normal(compactPoints));

    vector dir(boundaryPt - ownPt);

    if ((dir & n) < 0)
    {
        FatalErrorInFunction
            << "cell:" << cellI << " old face:" << faceI
            << " newFace:" << newFace
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " boundaryPt:" << boundaryPt
            << abort(FatalError);
    }

    vector fcToOwn(compactFace.centre(compactPoints) - ownPt);

    scalar s = (fcToOwn&dir) / magSqr(dir);

    if (s < 0.7 || s > 1.3)
    {
        WarningInFunction
            << "cell:" << cellI << " old face:" << faceI
            << " newFace:" << newFace
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " boundaryPt:" << boundaryPt
            << " s:" << s
            << endl;
            //<< abort(FatalError);
    }
}


// If p0 and p1 are existing vertices check if edge is split and insert
// splitPoint.
void Foam::hexTetRef8::insertEdgeSplit
(
    const labelList& edgeMidPoint,
    const label p0,
    const label p1,
    DynamicList<label>& verts
) const
{
    if (p0 < mesh_.nPoints() && p1 < mesh_.nPoints())
    {
        label edgeI = meshTools::findEdge(mesh_, p0, p1);

        if (edgeI != -1 && edgeMidPoint[edgeI] != -1)
        {
            verts.append(edgeMidPoint[edgeI]);
        }
    }
}


// Internal faces are one per edge between anchor points. So one per midPoint
// between the anchor points. Here we store the information on the midPoint
// and if we have enough information:
// - two anchors
// - two face mid points
// we add the face. Note that this routine can get called anywhere from
// two times (two unrefined faces) to four times (two refined faces) so
// the first call that adds the information creates the face.
Foam::label Foam::hexTetRef8::storeMidPointInfo
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& cellMidPoint,
    const labelList& edgeMidPoint,
    const label cellI,
    const label faceI,
    const bool faceOrder,
    const label edgeMidPointI,
    const label anchorPointI,
    const label faceMidPointI,

    Map<edge>& midPointToAnchors,
    Map<edge>& midPointToFaceMids,
    polyTopoChange& meshMod
) const
{
    // See if need to store anchors.

    bool changed = false;
    bool haveTwoAnchors = false;

    Map<edge>::iterator edgeMidFnd = midPointToAnchors.find(edgeMidPointI);

    if (edgeMidFnd == midPointToAnchors.end())
    {
        midPointToAnchors.insert(edgeMidPointI, edge(anchorPointI, -1));
    }
    else
    {
        edge& e = edgeMidFnd();

        if (anchorPointI != e[0])
        {
            if (e[1] == -1)
            {
                e[1] = anchorPointI;
                changed = true;
            }
        }

        if (e[0] != -1 && e[1] != -1)
        {
            haveTwoAnchors = true;
        }
    }

    bool haveTwoFaceMids = false;

    Map<edge>::iterator faceMidFnd = midPointToFaceMids.find(edgeMidPointI);

    if (faceMidFnd == midPointToFaceMids.end())
    {
        midPointToFaceMids.insert(edgeMidPointI, edge(faceMidPointI, -1));
    }
    else
    {
        edge& e = faceMidFnd();

        if (faceMidPointI != e[0])
        {
            if (e[1] == -1)
            {
                e[1] = faceMidPointI;
                changed = true;
            }
        }

        if (e[0] != -1 && e[1] != -1)
        {
            haveTwoFaceMids = true;
        }
    }

    // Check if this call of storeMidPointInfo is the one that completed all
    // the necessary information.

    if (changed && haveTwoAnchors && haveTwoFaceMids)
    {
        const edge& anchors = midPointToAnchors[edgeMidPointI];
        const edge& faceMids = midPointToFaceMids[edgeMidPointI];

        label otherFaceMidPointI = faceMids.otherVertex(faceMidPointI);

        // Create face consistent with anchorI being the owner.
        // Note that the edges between the edge mid point and the face mids
        // might be marked for splitting. Note that these edge splits cannot
        // be between cellMid and face mids.

        DynamicList<label> newFaceVerts(4);
        if (faceOrder == (mesh_.faceOwner()[faceI] == cellI))
        {
            newFaceVerts.append(faceMidPointI);

            // Check & insert edge split if any
            insertEdgeSplit
            (
                edgeMidPoint,
                faceMidPointI,  // edge between faceMid and
                edgeMidPointI,  // edgeMid
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointI);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointI,
                otherFaceMidPointI,
                newFaceVerts
            );

            newFaceVerts.append(otherFaceMidPointI);
            newFaceVerts.append(cellMidPoint[cellI]);
        }
        else
        {
            newFaceVerts.append(otherFaceMidPointI);

            insertEdgeSplit
            (
                edgeMidPoint,
                otherFaceMidPointI,
                edgeMidPointI,
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointI);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointI,
                faceMidPointI,
                newFaceVerts
            );

            newFaceVerts.append(faceMidPointI);
            newFaceVerts.append(cellMidPoint[cellI]);
        }

        face newFace;
        newFace.transfer(newFaceVerts);

        label anchorCell0 = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            cellI,
            faceI,
            anchorPointI
        );
        label anchorCell1 = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            cellI,
            faceI,
            anchors.otherVertex(anchorPointI)
        );


        label own, nei;
        point ownPt, neiPt;

        if (anchorCell0 < anchorCell1)
        {
            own = anchorCell0;
            nei = anchorCell1;

            ownPt = mesh_.points()[anchorPointI];
            neiPt = mesh_.points()[anchors.otherVertex(anchorPointI)];

        }
        else
        {
            own = anchorCell1;
            nei = anchorCell0;
            newFace.flip();

            ownPt = mesh_.points()[anchors.otherVertex(anchorPointI)];
            neiPt = mesh_.points()[anchorPointI];
        }

        if (debug)
        {
            point ownPt, neiPt;

            if (anchorCell0 < anchorCell1)
            {
                ownPt = mesh_.points()[anchorPointI];
                neiPt = mesh_.points()[anchors.otherVertex(anchorPointI)];
            }
            else
            {
                ownPt = mesh_.points()[anchors.otherVertex(anchorPointI)];
                neiPt = mesh_.points()[anchorPointI];
            }

            checkInternalOrientation
            (
                meshMod,
                cellI,
                faceI,
                ownPt,
                neiPt,
                newFace
            );
        }

        return addInternalFace
        (
            meshMod,
            faceI,
            anchorPointI,
            newFace,
            own,
            nei
        );
    }
    else
    {
        return -1;
    }
}


//- Get the relevant anchor for a tet-face on a hex-cell
void Foam::hexTetRef8::getRelevantAnchorAndMid
(
    label &relevantAnchor,
    label &relevantMid,
    const label &faceI,
    const cell &faceList,
    const label &referenceLevel,
    const DynamicList<label> &newPointLevel,
    const labelList &edgeMidPoint,
    const polyTopoChange &meshMod,
    const floatScalar &maxAngleBetweenFaces
) const
{
    // get face anchors
    labelList fAnchors(3, -1);
    {
        labelList fMids(3, -1);
        getFaceAnchorsAndMids
        (
            faceI,
            fAnchors,
            fMids,
            referenceLevel,
            newPointLevel,
            meshMod,
            false
        );
    }

    const pointField &points = mesh_.points();

    // plane of the face
    plane p1
    (
        points[fAnchors[0]],
        points[fAnchors[1]],
        points[fAnchors[2]]
    );

    label pairFace = -1;
    labelList commonAnchors(2, -1);

    // loop over all cell faces
    forAll(faceList, j)
    {
        const face &cF = mesh_.faces()[faceList[j]];

        // if the face is a tri-face
        if (faceList[j] != faceI && countAnchors(cF, referenceLevel) == 3)
        {
            label nCommonAnchors = 0;
            labelList commonAnchorsPotential(2, -1);

            forAll(cF, k)
            {
                if (whichAnchor(cF[k], fAnchors) != -1)
                {
                    commonAnchorsPotential[nCommonAnchors] = cF[k];
                    nCommonAnchors++;
                }
                if (nCommonAnchors == 2)
                {
                    break;
                }
            }

            // .. and if it has 2 anchors in common with
            // the face
            if (nCommonAnchors == 2)
            {
                // get face anchors
                labelList fCandidateAnchors(3, -1);
                {
                    labelList fMids(3, -1);
                    getFaceAnchorsAndMids
                    (
                        faceList[j],
                        fCandidateAnchors,
                        fMids,
                        referenceLevel,
                        newPointLevel,
                        meshMod,
                        false
                    );
                }

                plane p2
                (
                    points[fCandidateAnchors[0]],
                    points[fCandidateAnchors[1]],
                    points[fCandidateAnchors[2]]
                );

                // compute cos of angle between planes
                floatScalar cosOfAngle =
                    fabs(vectorTools::cosPhi(p1.normal(), p2.normal()));

                // if the angle b/w the faces is under a certain
                if (cosOfAngle >= cos(maxAngleBetweenFaces))
                {
                    // we have found the pairFace
                    pairFace = faceList[j];
                    commonAnchors = commonAnchorsPotential;
                }
            }
        }

        if (pairFace >= 0)
        {
            break;
        }
    }

    forAll(fAnchors, i)
    {
        if
        (
            fAnchors[i] != commonAnchors[0] &&
            fAnchors[i] != commonAnchors[1]
        )
        {
            relevantAnchor = fAnchors[i];
            break;
        }
    }

    const labelList &fEdges = mesh_.faceEdges()[faceI];
    forAll(fEdges, i)
    {
        const edge &e = mesh_.edges()[fEdges[i]];
        if
        (
            (e[0] == commonAnchors[0] && e[1] == commonAnchors[1]) ||
            (e[0] == commonAnchors[1] && e[1] == commonAnchors[0])
        )
        {
            relevantMid = edgeMidPoint[fEdges[i]];
            break;
        }
    }

    if (relevantAnchor == -1 || relevantMid == -1)
    {
        FatalErrorInFunction
            << "Relevant anchor / mid not found!"
            << abort(FatalError);
    }
}


// Creates all the 12 internal faces for cellI.
void Foam::hexTetRef8::createInternalFaces
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& cellMidPoint,
    const labelList& faceMidPoint,
    const labelList& faceAnchorLevel,
    const labelList& edgeMidPoint,
    const DynamicList<label> &newPointLevel,
    const label cellI,

    polyTopoChange& meshMod
) const
{
    // Find in every face the cellLevel+1 points (from edge subdivision)
    // and the anchor points.

    const cell& cFaces = mesh_.cells()[cellI];
    const label cLevel = cellLevel_[cellI];

    // From edge mid to anchor points
    Map<edge> midPointToAnchors(24);
    // From edge mid to face mids
    Map<edge> midPointToFaceMids(24);

    // Storage for on-the-fly addressing
    DynamicList<label> storage;


    // Running count of number of internal faces added so far.
    label nFacesAdded = 0;

    forAll(cFaces, i)
    {
        label faceI = cFaces[i];

        const face& f = mesh_.faces()[faceI];
        const labelList& fEdges = mesh_.faceEdges(faceI, storage);

        // We are on the cellI side of face f. The face will have 1, 3 or 4
        // cLevel points and lots of higher numbered ones.

        label faceMidPointI = -1;

        label nAnchors = countAnchors(f, cLevel);

        if (nAnchors == 1)
        {
            // Only one anchor point. So the other side of the face has already
            // been split using cLevel+1 and cLevel+2 points.

            // Find the one anchor.
            label anchorFp = -1;

            forAll(f, fp)
            {
                if (pointLevel_[f[fp]] <= cLevel)
                {
                    anchorFp = fp;
                    break;
                }
            }

            // Now the face mid point is the second cLevel+1 point
            label edgeMid = findLevel
            (
                faceI,
                f,
                f.fcIndex(anchorFp),
                true,
                cLevel+1
            );
            label faceMid = findLevel
            (
                faceI,
                f,
                f.fcIndex(edgeMid),
                true,
                cLevel+1
            );

            faceMidPointI = f[faceMid];
        }
        else if (nAnchors == 3)
        {
            // This is a triangular face. Find the corresp. other triangular face
            // this is done by looping over all cellIFaces -> if the face has 3
            // anchors then see the angle between this face and the other face.
            // If the angle is under a certain threshold, then this is the face.

            label relevantAnchor = -1; // not used
            getRelevantAnchorAndMid
            (
                relevantAnchor,
                faceMidPointI,
                faceI,
                mesh_.cells()[cellI],
                cellLevel_[cellI],
                newPointLevel,
                edgeMidPoint,
                meshMod,
                0.174533
            );
        }
        else if (nAnchors == 4)
        {
            // There is no face middle yet but the face will be marked for
            // splitting.

            faceMidPointI = faceMidPoint[faceI];
        }
        else
        {
            dumpCell(mesh_.faceOwner()[faceI]);
            if (mesh_.isInternalFace(faceI))
            {
                dumpCell(mesh_.faceNeighbour()[faceI]);
            }

            FatalErrorInFunction
                << "nAnchors:" << nAnchors
                << " faceI:" << faceI
                << abort(FatalError);
        }



        // Now loop over all the relevant anchors (might be just one) and
        // store the edge mids connected to it. storeMidPointInfo will collect
        // all the info and combine it all.

        forAll(f, fp0)
        {
            label point0 = f[fp0];

            if (pointLevel_[point0] <= cLevel)
            {
                // anchor

                // Walk forward
                // ~~~~~~~~~~~~
                // to cLevel+1 or edgeMidPoint of this level.


                label edgeMidPointI = -1;

                label fp1 = f.fcIndex(fp0);

                if
                (
                    (pointLevel_[f[fp1]] <= cLevel) &&
                    (f[fp1] != faceMidPointI)
                )
                {
                    // Anchor. Edge will be split.
                    label edgeI = fEdges[fp0];

                    edgeMidPointI = edgeMidPoint[edgeI];

                    if (edgeMidPointI == -1)
                    {
                        dumpCell(cellI);

                        const labelList& cPoints = mesh_.cellPoints(cellI);

                        FatalErrorInFunction
                            << "cell:" << cellI << " cLevel:" << cLevel
                            << " cell points:" << cPoints
                            << " pointLevel:"
                            << UIndirectList<label>(pointLevel_, cPoints)()
                            << " face:" << faceI
                            << " f:" << f
                            << " pointLevel:"
                            << UIndirectList<label>(pointLevel_, f)()
                            << " faceAnchorLevel:" << faceAnchorLevel[faceI]
                            << " faceMidPoint:" << faceMidPoint[faceI]
                            << " faceMidPointI:" << faceMidPointI
                            << " fp:" << fp0
                            << abort(FatalError);
                    }
                }
                else
                {
                    // Search forward in face to clevel+1
                    label edgeMid = findLevel(faceI, f, fp1, true, cLevel+1);

                    edgeMidPointI = f[edgeMid];
                }

                label newFaceI = storeMidPointInfo
                (
                    cellAnchorPoints,
                    cellAddedCells,
                    cellMidPoint,
                    edgeMidPoint,

                    cellI,
                    faceI,
                    true,                   // mid point after anchor
                    edgeMidPointI,          // edgemid
                    point0,                 // anchor
                    faceMidPointI,

                    midPointToAnchors,
                    midPointToFaceMids,
                    meshMod
                );

                if (newFaceI != -1)
                {
                    nFacesAdded++;

                    if (nFacesAdded == 12)
                    {
                        break;
                    }
                }



                // Walk backward
                // ~~~~~~~~~~~~~

                label fpMin1 = f.rcIndex(fp0);

                if
                (
                    (pointLevel_[f[fpMin1]] <= cLevel) &&
                    (f[fpMin1] != faceMidPointI)
                )
                {
                    // Anchor. Edge will be split.
                    label edgeI = fEdges[fpMin1];

                    edgeMidPointI = edgeMidPoint[edgeI];

                    if (edgeMidPointI == -1)
                    {
                        dumpCell(cellI);

                        const labelList& cPoints = mesh_.cellPoints(cellI);

                        FatalErrorInFunction
                            << "cell:" << cellI << " cLevel:" << cLevel
                            << " cell points:" << cPoints
                            << " pointLevel:"
                            << UIndirectList<label>(pointLevel_, cPoints)()
                            << " face:" << faceI
                            << " f:" << f
                            << " pointLevel:"
                            << UIndirectList<label>(pointLevel_, f)()
                            << " faceAnchorLevel:" << faceAnchorLevel[faceI]
                            << " faceMidPoint:" << faceMidPoint[faceI]
                            << " faceMidPointI:" << faceMidPointI
                            << " fp:" << fp0
                            << abort(FatalError);
                    }
                }
                else
                {
                    // Search back to clevel+1
                    label edgeMid = findLevel
                    (
                        faceI,
                        f,
                        fpMin1,
                        false,
                        cLevel+1
                    );

                    edgeMidPointI = f[edgeMid];
                }

                newFaceI = storeMidPointInfo
                (
                    cellAnchorPoints,
                    cellAddedCells,
                    cellMidPoint,
                    edgeMidPoint,

                    cellI,
                    faceI,
                    false,                  // mid point before anchor
                    edgeMidPointI,          // edgemid
                    point0,                 // anchor
                    faceMidPointI,

                    midPointToAnchors,
                    midPointToFaceMids,
                    meshMod
                );

                if (newFaceI != -1)
                {
                    nFacesAdded++;

                    if (nFacesAdded == 12)
                    {
                        break;
                    }
                }
            }   // done anchor
        }   // done face

        if (nFacesAdded == 12)
        {
            break;
        }
    }
}


void Foam::hexTetRef8::walkFaceToMid
(
    const labelList& edgeMidPoint,
    const label cLevel,
    const label faceI,
    const label startFp,
    DynamicList<label>& faceVerts
) const
{
    const face& f = mesh_.faces()[faceI];
    const labelList& fEdges = mesh_.faceEdges(faceI);

    label fp = startFp;

    // Starting from fp store all (1 or 2) vertices until where the face
    // gets split

    while (true)
    {
        if (edgeMidPoint[fEdges[fp]] >= 0)
        {
            faceVerts.append(edgeMidPoint[fEdges[fp]]);
        }

        fp = f.fcIndex(fp);

        if (pointLevel_[f[fp]] <= cLevel)
        {
            // Next anchor. Have already append split point on edge in code
            // above.
            return;
        }
        else if (pointLevel_[f[fp]] == cLevel+1)
        {
            // Mid level
            faceVerts.append(f[fp]);

            return;
        }
        else if (pointLevel_[f[fp]] == cLevel+2)
        {
            // Store and continue to cLevel+1.
            faceVerts.append(f[fp]);
        }
    }
}


// Same as walkFaceToMid but now walk back.
void Foam::hexTetRef8::walkFaceFromMid
(
    const labelList& edgeMidPoint,
    const label cLevel,
    const label faceI,
    const label startFp,
    DynamicList<label>& faceVerts
) const
{
    const face& f = mesh_.faces()[faceI];
    const labelList& fEdges = mesh_.faceEdges(faceI);

    label fp = f.rcIndex(startFp);

    while (true)
    {
        if (pointLevel_[f[fp]] <= cLevel)
        {
            // anchor.
            break;
        }
        else if (pointLevel_[f[fp]] == cLevel+1)
        {
            // Mid level
            faceVerts.append(f[fp]);
            break;
        }
        else if (pointLevel_[f[fp]] == cLevel+2)
        {
            // Continue to cLevel+1.
        }
        fp = f.rcIndex(fp);
    }

    // Store
    while (true)
    {
        if (edgeMidPoint[fEdges[fp]] >= 0)
        {
            faceVerts.append(edgeMidPoint[fEdges[fp]]);
        }

        fp = f.fcIndex(fp);

        if (fp == startFp)
        {
            break;
        }
        faceVerts.append(f[fp]);
    }
}


// Updates refineCell (cells marked for refinement) so across all faces
// there will be 2:1 consistency after refinement.
Foam::label Foam::hexTetRef8::faceConsistentRefinement
(
    const bool maxSet,
    PackedBoolList& refineCell
) const
{
    label nChanged = 0;

    // Internal faces.
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label own = mesh_.faceOwner()[faceI];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        label nei = mesh_.faceNeighbour()[faceI];
        label neiLevel = cellLevel_[nei] + refineCell.get(nei);

        label nOwnAnchors
            = countAnchors(mesh_.cellPoints()[own], cellLevel_[own]);

        label nNeiAnchors
            = countAnchors(mesh_.cellPoints()[nei], cellLevel_[nei]);

        bool isHexTetInterface = !(nNeiAnchors == nOwnAnchors);

        if
        (
            (ownLevel > (neiLevel+1))
         || ((ownLevel > neiLevel) && isHexTetInterface)
        )
        {
            if (maxSet)
            {
                refineCell.set(nei);
            }
            else
            {
                refineCell.unset(own);
            }

            nChanged++;
        }
        else if
        (
            (neiLevel > (ownLevel+1))
         || ((neiLevel > ownLevel) && isHexTetInterface)
        )
        {
            if (maxSet)
            {
                refineCell.set(own);
            }
            else
            {
                refineCell.unset(nei);
            }
            nChanged++;
        }
    }


    // Coupled faces. Swap owner level to get neighbouring cell level.
    // (only boundary faces of neiLevel used)
    // TODO: above fix here!
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    labelList nNeiAnchors(mesh_.nFaces()-mesh_.nInternalFaces());

    forAll(neiLevel, i)
    {
        label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

        neiLevel[i] = cellLevel_[own] + refineCell.get(own);

        nNeiAnchors[i]
            = countAnchors(mesh_.cellPoints()[own], cellLevel_[own]);
    }

    // Swap to neighbour
    syncTools::swapBoundaryFaceList(mesh_, neiLevel);
    syncTools::swapBoundaryFaceList(mesh_, nNeiAnchors);

    // Now we have neighbour value see which cells need refinement
    forAll(neiLevel, i)
    {
        label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        label nOwnAnchors
            = countAnchors(mesh_.cellPoints()[own], cellLevel_[own]);

        bool isHexTetInterface = !(nNeiAnchors[i] == nOwnAnchors);

        if ((ownLevel > (neiLevel[i]+1)) ||
            ((ownLevel > neiLevel[i]) && isHexTetInterface))
        {
            if (!maxSet)
            {
                refineCell.unset(own);
                nChanged++;
            }
        }
        else if ((neiLevel[i] > (ownLevel+1)) ||
                ((neiLevel[i] > ownLevel) && isHexTetInterface))
        {
            if (maxSet)
            {
                refineCell.set(own);
                nChanged++;
            }
        }
    }

    return nChanged;
}


// Debug: check if wanted refinement is compatible with 2:1
void Foam::hexTetRef8::checkWantedRefinementLevels
(
    const labelList& cellsToRefine
) const
{
    PackedBoolList refineCell(mesh_.nCells());
    forAll(cellsToRefine, i)
    {
        refineCell.set(cellsToRefine[i]);
    }

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label own = mesh_.faceOwner()[faceI];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        label nei = mesh_.faceNeighbour()[faceI];
        label neiLevel = cellLevel_[nei] + refineCell.get(nei);

        if (mag(ownLevel-neiLevel) > 1)
        {
            dumpCell(own);
            dumpCell(nei);
            FatalErrorInFunction
                << "cell:" << own
                << " current level:" << cellLevel_[own]
                << " level after refinement:" << ownLevel
                << nl
                << "neighbour cell:" << nei
                << " current level:" << cellLevel_[nei]
                << " level after refinement:" << neiLevel
                << nl
                << "which does not satisfy 2:1 constraints anymore."
                << abort(FatalError);
        }
    }

    // Coupled faces. Swap owner level to get neighbouring cell level.
    // (only boundary faces of neiLevel used)
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

    forAll(neiLevel, i)
    {
        label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

        neiLevel[i] = cellLevel_[own] + refineCell.get(own);
    }

    // Swap to neighbour
    syncTools::swapBoundaryFaceList(mesh_, neiLevel);

    // Now we have neighbour value see which cells need refinement
    forAll(neiLevel, i)
    {
        label faceI = i + mesh_.nInternalFaces();

        label own = mesh_.faceOwner()[faceI];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        if (mag(ownLevel - neiLevel[i]) > 1)
        {
            label patchI = mesh_.boundaryMesh().whichPatch(faceI);

            dumpCell(own);
            FatalErrorInFunction
                << "Celllevel does not satisfy 2:1 constraint."
                << " On coupled face "
                << faceI
                << " on patch " << patchI << " "
                << mesh_.boundaryMesh()[patchI].name()
                << " owner cell " << own
                << " current level:" << cellLevel_[own]
                << " level after refinement:" << ownLevel
                << nl
                << " (coupled) neighbour cell will get refinement "
                << neiLevel[i]
                << abort(FatalError);
        }
    }
}


// Set instance for mesh files
void Foam::hexTetRef8::setInstance(const fileName& inst)
{
    if (debug)
    {
        Pout<< "hexTetRef8::setInstance(const fileName& inst) : "
            << "Resetting file instance to " << inst << endl;
    }

    cellLevel_.instance() = inst;
    pointLevel_.instance() = inst;
    level0Edge_.instance() = inst;
    history_.instance() = inst;
}


void Foam::hexTetRef8::collectLevelPoints
(
    const labelList& f,
    const label level,
    DynamicList<label>& points
) const
{
    forAll(f, fp)
    {
        if (pointLevel_[f[fp]] <= level)
        {
            points.append(f[fp]);
        }
    }
}


void Foam::hexTetRef8::collectLevelPoints
(
    const labelList& meshPoints,
    const labelList& f,
    const label level,
    DynamicList<label>& points
) const
{
    forAll(f, fp)
    {
        label pointI = meshPoints[f[fp]];
        if (pointLevel_[pointI] <= level)
        {
            points.append(pointI);
        }
    }
}


// Return true if we've found 6 quads. faces guaranteed to be outwards pointing.
bool Foam::hexTetRef8::matchHexShape
(
    const label cellI,
    const label cellLevel,
    DynamicList<face>& quads
) const
{
    const cell& cFaces = mesh_.cells()[cellI];

    // Work arrays
    DynamicList<label> verts(4);
    quads.clear();


    // 1. pick up any faces with four cellLevel points

    forAll(cFaces, i)
    {
        label faceI = cFaces[i];
        const face& f = mesh_.faces()[faceI];

        verts.clear();
        collectLevelPoints(f, cellLevel, verts);
        if (verts.size() == 4)
        {
            if (mesh_.faceOwner()[faceI] != cellI)
            {
                reverse(verts);
            }
            quads.append(face(0));
            labelList& quadVerts = quads.last();
            quadVerts.transfer(verts);
        }
    }


    if (quads.size() < 6)
    {
        Map<labelList> pointFaces(2*cFaces.size());

        forAll(cFaces, i)
        {
            label faceI = cFaces[i];
            const face& f = mesh_.faces()[faceI];

            // Pick up any faces with only one level point.
            // See if there are four of these where the commont point
            // is a level+1 point. This common point is then the mid of
            // a split face.

            verts.clear();
            collectLevelPoints(f, cellLevel, verts);
            if (verts.size() == 1)
            {
                // Add to pointFaces for any level+1 point (this might be
                // a midpoint of a split face)
                forAll(f, fp)
                {
                    label pointI = f[fp];
                    if (pointLevel_[pointI] == cellLevel+1)
                    {
                        Map<labelList>::iterator iter =
                            pointFaces.find(pointI);
                        if (iter != pointFaces.end())
                        {
                            labelList& pFaces = iter();
                            if (findIndex(pFaces, faceI) == -1)
                            {
                                pFaces.append(faceI);
                            }
                        }
                        else
                        {
                            pointFaces.insert
                            (
                                pointI,
                                labelList(1, faceI)
                            );
                        }
                    }
                }
            }
        }

        // 2. Check if we've collected any midPoints.
        forAllConstIter(Map<labelList>, pointFaces, iter)
        {
            const labelList& pFaces = iter();

            if (pFaces.size() == 4)
            {
                // Collect and orient.
                faceList fourFaces(pFaces.size());
                forAll(pFaces, pFaceI)
                {
                    label faceI = pFaces[pFaceI];
                    const face& f = mesh_.faces()[faceI];
                    if (mesh_.faceOwner()[faceI] == cellI)
                    {
                        fourFaces[pFaceI] = f;
                    }
                    else
                    {
                        fourFaces[pFaceI] = f.reverseFace();
                    }
                }

                primitivePatch bigFace
                (
                    SubList<face>(fourFaces, fourFaces.size()),
                    mesh_.points()
                );
                const labelListList& edgeLoops = bigFace.edgeLoops();

                if (edgeLoops.size() == 1)
                {
                    // Collect the 4 cellLevel points
                    verts.clear();
                    collectLevelPoints
                    (
                        bigFace.meshPoints(),
                        bigFace.edgeLoops()[0],
                        cellLevel,
                        verts
                    );

                    if (verts.size() == 4)
                    {
                        quads.append(face(0));
                        labelList& quadVerts = quads.last();
                        quadVerts.transfer(verts);
                    }
                }
            }
        }
    }

    return (quads.size() == 6);
}

// ========================
// TET Refinement Functions
// ========================


//- Get anchors for given cell. The anchors are stored
//  such that the first three are in an outward pointing
//  order
Foam::labelList Foam::hexTetRef8::getAnchors
(
    label cell
)
{
    DynamicList<label> anchorPoints;
    label cellLevel  = cellLevel_[cell];

    labelList points = mesh_.cellPoints()[cell];

    forAll(points, i)
    {
        label pI = points[i];

        if (pointLevel_[pI] <= cellLevel)
        {
            anchorPoints.append(pI);
        }
    }

    if (anchorPoints.size() != 4)
    {
        FatalErrorInFunction
            << "Tet cells must have 4 anchor points. "
            << anchorPoints.size() << " found."
            << "cell: " << cell
            << "points: " << points
            << abort(FatalError);
    }

    {
        point p0 = mesh_.points()[anchorPoints[0]];
        point p1 = mesh_.points()[anchorPoints[1]];
        point p2 = mesh_.points()[anchorPoints[2]];
        point p3 = mesh_.points()[anchorPoints[3]];

        plane p(p0, p1, p2);

        if (p.sideOfPlane(p3) == 0)
        {
            label temp = anchorPoints[0];
            anchorPoints[0] = anchorPoints[1];
            anchorPoints[1] = temp;
        }
    }

    return anchorPoints.shrink();
}


//- Find the mid point between given anchors
Foam::label Foam::hexTetRef8::pointBetweenAnchors
(
    Pair<label> anchorPair,
    const labelList &anchors,
    label cell,
    const labelList &edgeMidPoint
)
{
    label midPoint = -1;

    // CASE 1 : do AS and AD form an undivided edge
    //          in the original mesh?

    const labelList &cellEdges = mesh_.cellEdges()[cell];
    forAll(cellEdges, i)
    {
        label cellEdgeI = cellEdges[i];
        edge e          = mesh_.edges()[cellEdgeI];

        if
        (
            (e[0] == anchorPair.first() && e[1] == anchorPair.second())
         || (e[1] == anchorPair.first() && e[0] == anchorPair.second())
        )
        {
            // edge found
            midPoint = edgeMidPoint[cellEdgeI];
            return midPoint;
        }
    }

    // CASE 2 : Do they form a part of a common face?

    const labelList &cellFaces = mesh_.cells()[cell];
    forAll(cellFaces, i)
    {
        label faceI = cellFaces[i];
        face f      = mesh_.faces()[faceI];

        forAll(f, j)
        {
            label pJ = f[j];

            if (anchorPair[0] == pJ || anchorPair[1] == pJ)
            {
                label next1Index = f.fcIndex(j);
                label next2Index = f.fcIndex(next1Index);
                label next3Index = f.fcIndex(next2Index);
                label next4Index = f.fcIndex(next3Index);

                if (debug)
                {
                    if (f[next1Index] == anchorPair.other(pJ))
                    {
                        FatalErrorInFunction
                        << "This case should have been covered in Case II"
                        << abort(FatalError);
                    }
                }

                if (f[next2Index] == anchorPair.other(pJ))
                {
                    // next1index can be the mid point
                    // check for level
                    if (pointLevel_[f[next1Index]] == cellLevel_[cell]+1)
                    {
                        return f[next1Index];
                    }
                }
                else if (f[next3Index] == anchorPair.other(pJ))
                {
                    if
                    (
                        (whichAnchor(f[next1Index], anchors) == -1) &&
                        (whichAnchor(f[next2Index], anchors) == -1)
                    )
                    {
                        // Two possibilities:
                        // 1. next1 is mid
                        // 2. next2 is mid

                        if (pointLevel_[f[next1Index]] == cellLevel_[cell]+1)
                        {
                            return f[next1Index];
                        }
                        else if
                        (
                            pointLevel_[f[next2Index]] == cellLevel_[cell]+1
                        )
                        {
                            return f[next2Index];
                        }
                    }
                }
                else if (f[next4Index] == anchorPair.other(pJ))
                {
                    if
                    (
                        (whichAnchor(f[next1Index], anchors) == -1) &&
                        (whichAnchor(f[next2Index], anchors) == -1) &&
                        (whichAnchor(f[next3Index], anchors) == -1)
                    )
                    // check if next2 is mid
                    if (pointLevel_[f[next2Index]] == cellLevel_[cell]+1)
                    {
                        return f[next2Index];
                    }
                }
            }
        }
    }

    // CASE 3 : Do the faces they are in have a point in common?

    {
        // gather facePoints for AnchorA and AnchorB
        FixedList<std::set<label>, 2> cellPointsConnectedToAnchor;

        forAll(cellFaces, i)
        {
            label faceI = cellFaces[i];
            face f      = mesh_.faces()[faceI];
            label nAnchors = countAnchors(f, cellLevel_[cell]);

            if (nAnchors < 3)
            {
                forAll(f, j)
                {
                    label pJ = f[j];

                    if (pJ == anchorPair[0])
                    {
                        forAll(f, k)
                        {
                            // avoid other anchors
                            if (whichAnchor(f[k], anchors) == -1)
                            {
                                cellPointsConnectedToAnchor[0].insert(f[k]);
                            }
                        }
                        break;
                    }
                    else if (pJ == anchorPair[1])
                    {
                        forAll(f, k)
                        {
                            // avoid other anchors
                            if (whichAnchor(f[k], anchors) == -1)
                            {
                                cellPointsConnectedToAnchor[1].insert(f[k]);
                            }
                        }
                        break;
                    }
                }
            }
        }

        // find common point between them
        std::set<label> commonPoint;

        std::set_intersection
        (
            cellPointsConnectedToAnchor[0].begin(),
            cellPointsConnectedToAnchor[0].end(),
            cellPointsConnectedToAnchor[1].begin(),
            cellPointsConnectedToAnchor[1].end(),
            std::inserter(commonPoint, commonPoint.begin())
        );

        if (commonPoint.size() != 1)
        {
            FatalErrorInFunction
            << "Point sets for anchorA and anchorB have "
            << "more than one / no point in common"
            << abort(FatalError);
        }
        else
        {
            return *commonPoint.begin();
        }
    }

    FatalErrorInFunction
    << "Failed to find mid-point between anchors."
    << abort(FatalError);

    return -1;
}


//- Find anchor toes: these are those points (new or old) that,
//  for each anchor, are the midPoints of the edges that connect
//  that anchor to other anchors.
Foam::labelListList Foam::hexTetRef8::getAnchorToes
(
    label cell,
    const labelList& anchors,
    const labelList& edgeMidPoint
)
{

    labelListList toes(4, labelList(4, -1));

    forAll(anchors, i)
    {
        label AS = anchors[i];  // source anchor
        label ad = anchors.fcIndex(i);  // destination anchor index

        for (int j = i+1; j < 4; j++)
        {
            label AD = anchors[ad];

            Pair<label> anchorPair(AS, AD);
            label mid = pointBetweenAnchors
                        (
                            anchorPair,
                            anchors,
                            cell,
                            edgeMidPoint
                        );

            if (whichAnchor(mid, anchors) != -1)
            {
                FatalErrorInFunction
                    << "Mid point is an anchor. This must not happen."
                    << abort(FatalError);
            }

            toes[i][j] = mid;
            toes[j][i] = mid;

            ad = anchors.fcIndex(ad);
        }
    }

    return toes;
}


// Gets cell level such that the face has three points <= level.
Foam::label Foam::hexTetRef8::faceLevelTet
(
    const label faceI
) const
{
    const face& f = mesh_.faces()[faceI];

    if (f.size() < 3)
    {
        FatalErrorInFunction
        << "Face " << faceI << " has less than 3 points."
        << abort(FatalError);

        return -1;
    }
    else if (f.size() == 3)
    {
        return pointLevel_[f[findMaxLevel(f)]];
    }
    else
    {
        label ownLevel = cellLevel_[mesh_.faceOwner()[faceI]];

        if (countAnchors(f, ownLevel) == 3)
        {
            return ownLevel;
        }
        else if (countAnchors(f, ownLevel+1) == 3)
        {
            return ownLevel+1;
        }
        else
        {
            return -1;
        }
    }
}


//- Find if point is anchor. If yes, get anchorIndex.
//  If not, return -1
Foam::label Foam::hexTetRef8::whichAnchor
(
    label pointID,
    const labelList &anchors
) const
{
    forAll(anchors, i)
    {
        if (anchors[i] == pointID)
        {
            return i;
        }
    }

    return -1;
}


//- Get the centerFace that face f represents.
//  CellAddedCells contains 8 sub-cells for each cell:
//      -> the first four are cells adjacent to the
//         anchors (in order)
//      -> the next four are the cells adjacent to the
//         center faces of the cellFaces
//              -- b/w anchors 0-1-2
//              -- b/w anchors 0-1-3
//              -- b/w anchors 0-2-3
//              -- b/w anchors 1-2-3
//      -> the last centerFace to be modded recycles
//         the original cell
Foam::label Foam::hexTetRef8::whichCenterFace
(
    const face &f,
    const labelListList &anchorToes,
    label referenceLevel,
    const DynamicList<label> &newPointLevel
)
{
    labelHashSet faceAnchors;
    label nAnchorsFound = 0;
    forAll(f, i)
    {
        if (newPointLevel[f[i]] <= referenceLevel)
        {
            faceAnchors.insert(f[i]);
            nAnchorsFound++;
        }
    }

    if (nAnchorsFound != 3)
    {
        FatalErrorInFunction
            << "Face must have 3 toes, but "
            << "it has only " << nAnchorsFound
            << abort(FatalError);
    }

    if
    (
        faceAnchors[anchorToes[0][1]] &&
        faceAnchors[anchorToes[1][2]] &&
        faceAnchors[anchorToes[2][0]]
    )
    {
        return 0;
    }

    else if
    (
        faceAnchors[anchorToes[0][1]] &&
        faceAnchors[anchorToes[1][3]] &&
        faceAnchors[anchorToes[3][0]]
    )
    {
        return 1;
    }

    else if
    (
        faceAnchors[anchorToes[0][2]] &&
        faceAnchors[anchorToes[2][3]] &&
        faceAnchors[anchorToes[3][0]]
    )
    {
        return 2;
    }

    else if
    (
        faceAnchors[anchorToes[1][2]] &&
        faceAnchors[anchorToes[2][3]] &&
        faceAnchors[anchorToes[3][1]]
    )
    {
        return 3;
    }

    else
    {
        return -1;
    }
}


//- For given faceID and cellID, find the subCell that
//  that face lies next to.
//  CellAddedCells contains 8 sub-cells for each cell:
//      -> the first four are cells adjacent to the
//         anchors (in order)
//      -> the next four are the cells adjacent to the
//         center faces of the cellFaces
//              -- b/w anchors 0-1-2
//              -- b/w anchors 0-1-3
//              -- b/w anchors 0-2-3
//              -- b/w anchors 1-2-3
//      -> the last centerFace to be modded recycles
//         the original cell
Foam::label Foam::hexTetRef8::getSubCellID
(
    const face &f,
    label parentFaceID,
    label referenceLevel,
    const DynamicList<label> &newPointLevel,
    const labelListList &anchorToes,
    label cellID,
    const labelList &anchors,
    const labelListList &cellAddedCells,
    polyTopoChange &meshMod
)
{
    // get number of sub-anchors in face
    label nSubAnchorsInF = countAnchors(f, referenceLevel, newPointLevel);

    face fInFocus;

    if (nSubAnchorsInF <= 1)
    {
        // it is a gen-II face. check if parent face
        // has any cellAnchors
        fInFocus = meshMod.faces()[parentFaceID];
    }
    else if (nSubAnchorsInF >= 3)
    {
        // it is a gen-I face. check if it has any cellAnchors
        fInFocus = f;
    }
    else
    {
        FatalErrorInFunction
            << "Face has " << nSubAnchorsInF << " anchors of level "
            << referenceLevel << " , which must be 0, 1 or 3."
            << abort(FatalError);
    }

    label nAnchorsInFocusFace
        = countAnchors(fInFocus, referenceLevel-1, newPointLevel);

    if (nAnchorsInFocusFace == 1)
    {
        // it is an anchorFace. find anchor index and
        // return the anchor cell corresp. to it
        forAll(fInFocus, i)
        {
            label anchorIndex = whichAnchor(fInFocus[i], anchors);
            if (anchorIndex != -1)
            {
                return cellAddedCells[cellID][anchorIndex];
            }
        }
    }
    else if (nAnchorsInFocusFace == 0)
    {
        // return the centerCell corresponding to the
        // parent face
        label centerFaceIndex = whichCenterFace
        (
            fInFocus,
            anchorToes,
            cellLevel_[cellID]+1,
            newPointLevel
        );

        if (centerFaceIndex != -1)
        {
            return cellAddedCells[cellID][4+centerFaceIndex];
        }
        else
        {
            FatalErrorInFunction
                << "Center face not matched for face "
                << fInFocus << nl
                << "Face Anchors: " << anchorToes
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Face has " << nAnchorsInFocusFace << " anchors of level "
            << referenceLevel-1 << " , which must be 1 or 0."
            << abort(FatalError);
    }

    return -1;
}


//- for given face, find its anchors and the mid points
//  between them
void Foam::hexTetRef8::getFaceAnchorsAndMids
(
    label faceID,
    labelList &fAnchors,
    labelList &fMids,
    label fLevel,
    const DynamicList<label> &newPointLevel,
    const polyTopoChange &meshMod,
    const bool midsWanted
) const
{
    face f;
    if (midsWanted)
    {
        f = meshMod.faces()[faceID];
    }
    else
    {
        f = mesh_.faces()[faceID];
    }

    // find face anchors
    label firstAnchorIndex = -1;
    label nAnchorsFound = 0;

    forAll(f, j)
    {
        if (newPointLevel[f[j]] <= fLevel)
        {
            fAnchors[nAnchorsFound++] = f[j];
            if (firstAnchorIndex == -1)
            {
                firstAnchorIndex = j;
            }
        }
    }

    if (firstAnchorIndex == -1)
    {
        FatalErrorInFunction
            << "No anchors with level " << fLevel
            << " found for face " << faceID << " " << f
            << abort(FatalError);
    }

    if (midsWanted)
    {
        // find face mids
        label nMidsFound = 0;
        label k = firstAnchorIndex;

        forAll(f, j)
        {
            if (newPointLevel[f[k]] == fLevel+1)
            {
                fMids[nMidsFound++] = f[k];
            }
            k = f.fcIndex(k);
        }

        if (nAnchorsFound != 3 || nMidsFound != 3)
        {
            FatalErrorInFunction
            << "Number of anchors/mids for face " << faceID << " is not 3."
            << abort(FatalError);
        }
    }
    else
    {
        if (nAnchorsFound != 3)
        {
            FatalErrorInFunction
            << "Number of anchors for face " << faceID << " is not 3."
            << abort(FatalError);
        }
    }
}


//- get the right own/nei
void Foam::hexTetRef8::getNeighbours
(
    const face &f,
    label parentFaceID,
    label cellID,
    const labelList &anchors,
    const labelListList &anchorToes,
    const labelListList &cellAddedCells,
    label &ownNew,
    label &neiNew,
    const labelList &edgeMidPoint,
    const DynamicList<label> &newPointLevel,
    polyTopoChange &meshMod
)
{
    label ownOld = meshMod.faceOwner()[parentFaceID];
    label neiOld = meshMod.faceNeighbour()[parentFaceID];

    label subCell = getSubCellID
    (
        f,
        parentFaceID,
        cellLevel_[cellID]+1,
        newPointLevel,
        anchorToes,
        cellID,
        anchors,
        cellAddedCells,
        meshMod
    );

    if (neiOld == -1)
    {
        ownNew = subCell;
        neiNew = neiOld;
        return;
    }
    else if (ownOld != cellID && neiOld != cellID)
    {
        ownNew = ownOld;
        neiNew = neiOld;
        return;
    }

    // get adjacent cell

    label adjacentCell = -1;

    if (ownOld == cellID)
    {
        adjacentCell = neiOld;
    }
    else
    {
        adjacentCell = ownOld;
    }

    // find the ID of the correct adjacent subCell

    label adjacentSubCellID = -1;

    if (adjacentCell >= mesh_.nCells())
    {
        adjacentSubCellID = adjacentCell;   // already looking at a subcell
    }
    else
    {
        const labelList &adjacentCellAddedCells
            = cellAddedCells[adjacentCell];

        if (adjacentCellAddedCells.size() != 8) // adjacent cell not split
        {
            adjacentSubCellID = adjacentCell;
        }
        else
        {
            const labelList &cPoints = mesh_.cellPoints()[adjacentCell];
            label nCAnchors = countAnchors(cPoints, cellLevel_[adjacentCell]);

            if (nCAnchors == 4)
            {
                // tet cell

                labelList adjacentCellAnchors
                    = getAnchors(adjacentCell);

                labelListList adjacentCellAnchorToes
                    = getAnchorToes
                      (
                          adjacentCell,
                          adjacentCellAnchors,
                          edgeMidPoint
                      );

                adjacentSubCellID = getSubCellID
                (
                    f,
                    parentFaceID,
                    cellLevel_[adjacentCell]+1,
                    newPointLevel,
                    adjacentCellAnchorToes,
                    adjacentCell,
                    adjacentCellAnchors,
                    cellAddedCells,
                    meshMod
                );
            }
            else if (nCAnchors == 8)
            {
                // find if face has an anchor from the neighbouring cell
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                labelList cAnchors(8, -1);
                label nCAnchorsFound = 0;

                forAll(cPoints, i)
                {
                    if (pointLevel_[cPoints[i]] <= cellLevel_[adjacentCell])
                    {
                        cAnchors[nCAnchorsFound++] = cPoints[i];
                    }
                }

                forAll(f, i)
                {
                    forAll(cAnchors, j)
                    {
                        if (cAnchors[j] == f[i])
                        {
                            adjacentSubCellID = cellAddedCells[adjacentCell][j];
                        }
                    }
                }

                // if not, find the relevant anchor point in the parent face
                // and return the cell corresponding to it
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                if (adjacentSubCellID == -1)
                {
                    label relevantAnchor = -1;
                    label relevantMid = -1;

                    getRelevantAnchorAndMid
                    (
                        relevantAnchor,
                        relevantMid,
                        parentFaceID,
                        mesh_.cells()[adjacentCell],
                        cellLevel_[cellID],
                        newPointLevel,
                        edgeMidPoint,
                        meshMod,
                        0.174533
                    );

                    forAll(cAnchors, j)
                    {
                        if (cAnchors[j] == relevantAnchor)
                        {
                            adjacentSubCellID = cellAddedCells[adjacentCell][j];
                        }
                    }
                }
            }
        }
    }

    if (ownOld == cellID)
    {
        ownNew = subCell;
        neiNew = adjacentSubCellID;
    }
    else
    {
        ownNew = adjacentSubCellID;
        neiNew = subCell;
    }
}


//- Compute the cartesian distance between
//  two points
Foam::doubleScalar Foam::hexTetRef8::cartDistance
(
    point p1,
    point p2
)
{
    doubleScalar dist =
        sqrt
        (
            (p2.x() - p1.x()) * (p2.x() - p1.x()) +
            (p2.y() - p1.y()) * (p2.y() - p1.y()) +
            (p2.z() - p1.z()) * (p2.z() - p1.z())
        );

    return dist;
}


//- Choose which two opposite mids to connect
//  such that tets are of best quality. Ref
//  Bey 1995
Foam::Pair<Foam::label> Foam::hexTetRef8::chooseDiagonal
(
    const labelListList &anchorToes,
    polyTopoChange &meshMod
)
{
    // Code should be compacted, but okay for now

    doubleScalar distMin = doubleScalarGREAT;
    Pair<label> diag(-1, -1);
    doubleScalar dist;
    label S, F;

    // Case 1: [A1, A2] <-> [A3, A4]
    {
        S = anchorToes[0][1];
        F = anchorToes[2][3];

        point pS = meshMod.points()[S];
        point pF = meshMod.points()[F];

        dist = cartDistance(pS, pF);
    }

    if (dist < distMin)
    {
        distMin = dist;
        diag.first() = F;
        diag.second() = S;
    }

    // Case 2: [A1, A3] <-> [A2, A4]
    {
        S = anchorToes[0][2];
        F = anchorToes[1][3];

        point pS = meshMod.points()[S];
        point pF = meshMod.points()[F];

        dist = cartDistance(pS, pF);
    }

    if (dist < distMin)
    {
        distMin = dist;
        diag.first() = F;
        diag.second() = S;
    }

    // Case 3: [A1, A4] <-> [A2, A3]
    {
        S = anchorToes[0][3];
        F = anchorToes[1][2];

        point pS = meshMod.points()[S];
        point pF = meshMod.points()[F];

        dist = cartDistance(pS, pF);
    }

    if (dist < distMin)
    {
        distMin = dist;
        diag.first() = F;
        diag.second() = S;
    }

    return diag;
}


//- From the list of centerFaces, return the indices
//  of those that contain mid. Also give the corresp.
//  centerCells and face points
void Foam::hexTetRef8::getCenterFacesWithPoint
(
    const faceList &centerFaces,
    label mid,
    const labelList &cellIAddedCells,
    faceList &mCF,
    labelList &mCFIndices,
    labelList &mCFCells
)
{
    label nCenterFacesWithPointFound = 0;

    forAll(centerFaces, i)
    {
        face f = centerFaces[i];

        if (f.which(mid) != -1)
        {
            mCFIndices[nCenterFacesWithPointFound] = i;
            mCF[nCenterFacesWithPointFound] = f;
            nCenterFacesWithPointFound++;
        }
    }

    if (nCenterFacesWithPointFound != 2)
    {
        FatalErrorInFunction
            << "Anchor mid point " << mid << " is in "
            << nCenterFacesWithPointFound << " faces."
            << " Should be in 2 faces."
            << abort(FatalError);
    }

    forAll(mCFIndices, i)
    {
        mCFCells[i] = cellIAddedCells[4 + mCFIndices[i]];
    }
}


//- Given two points, find if there is a point between
//  them using the list of faces provided.
//  Note: referenceLevel is the level of the toes
Foam::label Foam::hexTetRef8::walkBetweenPoints
(
    label p1,
    label p2,
    const faceList &searchFaceList,
    const DynamicList<label> &newPointLevel,
    label referenceLevel,
    DynamicList<label> &newFace,
    bool isOrderImportantForSearch,
    bool isOrderImportantForAddition
)
{
    DynamicList<label> points;
    bool pointsFound = false;

    forAll(searchFaceList, i)
    {
        face f = searchFaceList[i];

        forAll (f, j)
        {
            if
            (
                isOrderImportantForSearch
                ? (f[j] == p1)
                : (f[j] == p1 || f[j] == p2)
            )
            {
                label other = (f[j] == p1) ? p2 : p1;
                label next1Index = f.fcIndex(j);
                label next2Index = f.fcIndex(next1Index);
                label next3Index = f.fcIndex(next2Index);
                label next4Index = f.fcIndex(next3Index);

                // check if next/previous point is p1/p2
                if (f[next1Index] == other)
                {
                    points.append(f[j]);
                    points.append(f[next1Index]);
                    pointsFound = true;
                    break;
                }
                else if (f[next2Index] == other)
                {
                    if (newPointLevel[f[next1Index]] > referenceLevel)
                    {
                        points.append(f[j]);
                        points.append(f[next1Index]);
                        points.append(f[next2Index]);
                        pointsFound = true;
                        break;
                    }
                }
                else if (f[next3Index] == other)
                {
                    if
                    (
                        newPointLevel[f[next1Index]] > referenceLevel &&
                        newPointLevel[f[next2Index]] > referenceLevel
                    )
                    {
                        points.append(f[j]);
                        points.append(f[next1Index]);
                        points.append(f[next2Index]);
                        points.append(f[next3Index]);
                        pointsFound = true;
                        break;
                    }
                }
                else if (f[next4Index] == other)
                {
                    if
                    (
                        newPointLevel[f[next1Index]] > referenceLevel &&
                        newPointLevel[f[next2Index]] > referenceLevel &&
                        newPointLevel[f[next3Index]] > referenceLevel
                    )
                    {
                        points.append(f[j]);
                        points.append(f[next1Index]);
                        points.append(f[next2Index]);
                        points.append(f[next3Index]);
                        points.append(f[next4Index]);
                        pointsFound = true;
                        break;
                    }
                }
            }
            if (pointsFound)
            {
                break;
            }
        }
        if (pointsFound)
        {
            break;
        }
    }

    if (!pointsFound)
    {
        // this means that the points are not yet connected.
        points.append(p1);
        points.append(p2);
    }

    if (newFace.size() == 0)
    {
        // this is the first time we are adding points to the face
        // thus append all points (including the first and last)
        if (isOrderImportantForAddition)
        {
            // it is imperative that we store the points in the
            // same order as requested in the input
            if (points[0] == p1)
            {
                // order is correct
                for (label i = 0; i < points.size(); i++)
                {
                    newFace.append(points[i]);
                }
            }
            else
            {
                // order is reversed
                for (label i = points.size()-1; i >= 0 ; i--)
                {
                    newFace.append(points[i]);
                }
            }
        }
        else
        {
            // store the points in the order in which
            // they were found
            newFace.append(points);
        }
    }
    else
    {
        // never add the first point, since it is already added
        if (points.size() > 2)
        {
            // there are points in between. add middle points
            if (isOrderImportantForAddition)
            {
                if (points[0] == p1)
                {
                    // order is correct
                    for (label i = 1; i < points.size()-1; i++)
                    {
                        newFace.append(points[i]);
                    }
                }
                else
                {
                    // order is reversed
                    for (label i = points.size()-2; i > 0 ; i--)
                    {
                        newFace.append(points[i]);
                    }
                }
            }
            else
            {
                for (label i = 1; i < points.size()-1; i++)
                {
                    newFace.append(points[i]);
                }
            }
        }

        // if the last point in points is not the same as
        // the first point in the face, add it too.
        if (isOrderImportantForAddition)
        {
            if (newFace[0] != p2)
            {
                newFace.append(p2);
            }
        }
        else
        {
            if (newFace[0] != points.last())
            {
                newFace.append(points.last());
            }
        }
    }

    // return the last point of the face
    return newFace.last();
}


Foam::face Foam::hexTetRef8::getNewInternalFace
(
    const faceList &searchFaces,
    const face &mCF,
    label mCFCell,
    label adjacentCell,
    label mid,
    Pair<label> diag,
    bool faceOutOfMCFCell,
    label &ownNew,
    label &neiNew,
    label option,
    const DynamicList<label> &newPointLevel,
    label level,
    polyTopoChange &meshMod
)
{
    label diagPointInMCF    = mCF.which(diag[0]) == -1 ? diag[1] : diag[0];
    label diagPointNotInMCF = mCF.which(diag[0]) == -1 ? diag[0] : diag[1];

    DynamicList<label> newFaceVerts;

    if (option == 0)
    {
        // add points from mid to diagPointInMFC (right order)
        label lastPointAdded =
        walkBetweenPoints
        (
            mid,
            diagPointInMCF,
            searchFaces,
            newPointLevel,
            level,
            newFaceVerts,
            false,
            false
        );

        // add points from lastPointAdded to diagPointNotInMCF
        walkBetweenPoints
        (
            lastPointAdded,
            diagPointNotInMCF,
            searchFaces,
            newPointLevel,
            level,
            newFaceVerts,
            false,
            true
        );

        // add points from diagPointNotInMCF to otherPoint
        walkBetweenPoints
        (
            diagPointNotInMCF,
            lastPointAdded == mid ? diagPointInMCF : mid,
            searchFaces,
            newPointLevel,
            level,
            newFaceVerts,
            false,
            true
        );
    }
    else if (option == 1)
    {
        // find the other point in mCF that has the same level as mid
        label nonDiagPoint = -1;
        forAll(mCF, i)
        {
            if ((mCF[i] != mid) && (mCF[i] != diagPointInMCF))
            {
                if (newPointLevel[mCF[i]] <= level)
                {
                    nonDiagPoint = mCF[i];
                }
            }
        }
        if (nonDiagPoint == -1)
        {
            FatalErrorInFunction
                << "Other face anchor not found"
                << abort(FatalError);
        }

        // add points from mid to nonDiagPoint (right order)
        label lastPointAdded =
        walkBetweenPoints
        (
            mid,
            nonDiagPoint,
            searchFaces,
            newPointLevel,
            level,
            newFaceVerts,
            false,
            false
        );

        // add points from lastPointAdded to diagPointNotInMFC
        walkBetweenPoints
        (
            lastPointAdded,
            diagPointNotInMCF,
            searchFaces,
            newPointLevel,
            level,
            newFaceVerts,
            false,
            true
        );

        // add points from diagPointNotInMFC to otherPoint
        walkBetweenPoints
        (
            diagPointNotInMCF,
            lastPointAdded == mid ? nonDiagPoint : mid,
            searchFaces,
            newPointLevel,
            level,
            newFaceVerts,
            false,
            true
        );
    }
    else
    {
        FatalErrorInFunction
            << "Option not recognized"
            << abort(FatalError);
    }

    face newFace;
    newFace.transfer(newFaceVerts);

    // point order points in outward direction
    if (faceOutOfMCFCell)   // need to take reverse
    {
        newFace = newFace.reverseFace();

        ownNew = mCFCell;
        neiNew = adjacentCell;
    }
    else    // keep order
    {
        ownNew = adjacentCell;
        neiNew = mCFCell;
    }

    return newFace;
}


//- Return diagonal point that is NOT
//  present in face
Foam::label Foam::hexTetRef8::findOtherDiagPoint
(
    Pair<label> diag,
    face f
)
{
    forAll(f, i)
    {
        if (f[i] == diag.first())
        {
            return diag.second();
        }
    }

    return diag.first();
}


//- Find the right midCenterFace for adding internal
//  anchor face for given anchor
Foam::label Foam::hexTetRef8::findAnchorMCF
(
    const faceList& mCF,
    const labelList& anchorAToes,
    Pair<label> diag
)
{
    // find which diagPoint lies in anchorToes

    label diagPointInAnchorToes = -1;

    forAll(anchorAToes, i)
    {
        if(anchorAToes[i] == diag.first() ||
            anchorAToes[i] == diag.second())
        {
            diagPointInAnchorToes = anchorAToes[i];
            break;
        }
    }

    // check which midCenterFace contains the diagPoint
    // and return the other face

    forAll(mCF, i)
    {
        const face &mCFI = mCF[i];

        forAll(mCFI, j)
        {
            if (mCFI[j] == diagPointInAnchorToes)
            {
                return mCF.fcIndex(i);
            }
        }
    }

    return -1;
}


//- Find the cell adjacent to the given face of the
//  given cell in the mesh. Returns -1 if boundary cell
Foam::label Foam::hexTetRef8::getAdjacentCell
(
    label face,
    label cell
)
{
    label adjacentCell = -1;
    if (mesh_.faceOwner()[face] != cell)
    {
        adjacentCell = mesh_.faceOwner()[face];
    }
    else
    {
        if (mesh_.isInternalFace(face))
        {
            adjacentCell = mesh_.faceNeighbour()[face];
        }
    }

    return adjacentCell;
}


//- For a given cell, get its centerFaces in the right
//  order.
//      -> b/w anchors 0-1-2
//      -> b/w anchors 0-1-3
//      -> b/w anchors 0-2-3
//      -> b/w anchors 1-2-3
Foam::faceList Foam::hexTetRef8::getCenterFaces
(
    const labelListList &anchorToes,
    const faceList &cellFaces,
    const DynamicList<label> &newPointLevel,
    const label referenceLevel
)
{
    faceList centerFaces(4);

    labelListList pOrd(4, labelList(3, -1));
    // for centerFace0: 0,1   1,2   2,0
    pOrd[0][0] = anchorToes[0][1];
    pOrd[0][1] = anchorToes[1][2];
    pOrd[0][2] = anchorToes[2][0];

    // for centerFace1: 0,3   3,1   1,0
    pOrd[1][0] = anchorToes[0][3];
    pOrd[1][1] = anchorToes[3][1];
    pOrd[1][2] = anchorToes[1][0];

    // for centerFace2: 0,2   2,3   3,0
    pOrd[2][0] = anchorToes[0][2];
    pOrd[2][1] = anchorToes[2][3];
    pOrd[2][2] = anchorToes[3][0];

    // for centerFace3: 1,3   3,2   2,1
    pOrd[3][0] = anchorToes[1][3];
    pOrd[3][1] = anchorToes[3][2];
    pOrd[3][2] = anchorToes[2][1];

    forAll(centerFaces, i)
    {
        DynamicList<label> faceVerts;
        forAll(pOrd[i], j)
        {
            walkBetweenPoints
            (
                pOrd[i][j],
                pOrd[i][pOrd[i].fcIndex(j)],
                cellFaces,
                newPointLevel,
                referenceLevel,
                faceVerts,
                false,
                true
            );
        }
        centerFaces[i].transfer(faceVerts);
    }

    return centerFaces;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh, read refinement data
Foam::hexTetRef8::hexTetRef8(const polyMesh& mesh, const bool readHistory)
:
    mesh_(mesh),
    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        labelList(mesh_.nCells(), 0)
    ),
    pointLevel_
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        labelList(mesh_.nPoints(), 0)
    ),
    level0Edge_
    (
        IOobject
        (
            "level0Edge",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar("level0Edge", dimLength, getLevel0EdgeLength())
    ),
    history_
    (
        IOobject
        (
            "refinementMRAHistory",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (readHistory ? mesh_.nCells() : 0)  // All cells visible if not be read
    ),
    faceRemover_(mesh_, GREAT),     // merge boundary faces wherever possible
    savedPointLevel_(0),
    savedCellLevel_(0),
    tetTesting_(false)
{
    if (readHistory)
    {
        // Make sure we don't use the master-only reading. Bit of a hack for
        // now.
        regIOobject::fileCheckTypes oldType =
            regIOobject::fileModificationChecking;
        regIOobject::fileModificationChecking = regIOobject::timeStamp;
        history_.readOpt() = IOobject::READ_IF_PRESENT;
        if (history_.headerOk())
        {
            history_.read();
        }
        regIOobject::fileModificationChecking = oldType;
    }

    if (history_.active() && history_.visibleCells().size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "History enabled but number of visible cells "
            << history_.visibleCells().size() << " in "
            << history_.objectPath()
            << " is not equal to the number of cells in the mesh "
            << mesh_.nCells()
            << abort(FatalError);
    }

    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorInFunction
            << "Restarted from inconsistent cellLevel or pointLevel files."
            << endl
            << "cellLevel file " << cellLevel_.objectPath() << endl
            << "pointLevel file " << pointLevel_.objectPath() << endl
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size() << endl
            << "Number of points in mesh:" << mesh_.nPoints()
            << " does not equal size of pointLevel:" << pointLevel_.size()
            << abort(FatalError);
    }


    // Check refinement levels for consistency
    checkRefinementLevels(-1, labelList(0));


    // Check initial mesh for consistency

    //if (debug)
    {
        checkMesh();
    }
}


// Construct from components
Foam::hexTetRef8::hexTetRef8
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const refinementMRAHistory& history,
    const scalar level0Edge
)
:
    mesh_(mesh),
    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        cellLevel
    ),
    pointLevel_
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointLevel
    ),
    level0Edge_
    (
        IOobject
        (
            "level0Edge",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar
        (
            "level0Edge",
            dimLength,
            (level0Edge >= 0 ? level0Edge : getLevel0EdgeLength())
        )
    ),
    history_
    (
        IOobject
        (
            "refinementMRAHistory",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        history
    ),
    faceRemover_(mesh_, GREAT),     // merge boundary faces wherever possible
    savedPointLevel_(0),
    savedCellLevel_(0),
    tetTesting_(false)
{
    if (history_.active() && history_.visibleCells().size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "History enabled but number of visible cells in it "
            << history_.visibleCells().size()
            << " is not equal to the number of cells in the mesh "
            << mesh_.nCells() << abort(FatalError);
    }

    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorInFunction
            << "Incorrect cellLevel or pointLevel size." << endl
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size() << endl
            << "Number of points in mesh:" << mesh_.nPoints()
            << " does not equal size of pointLevel:" << pointLevel_.size()
            << abort(FatalError);
    }

    // Check refinement levels for consistency
    checkRefinementLevels(-1, labelList(0));


    // Check initial mesh for consistency

    //if (debug)
    {
        checkMesh();
    }
}


// Construct from components
Foam::hexTetRef8::hexTetRef8
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const scalar level0Edge
)
:
    mesh_(mesh),
    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        cellLevel
    ),
    pointLevel_
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointLevel
    ),
    level0Edge_
    (
        IOobject
        (
            "level0Edge",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar
        (
            "level0Edge",
            dimLength,
            (level0Edge >= 0 ? level0Edge : getLevel0EdgeLength())
        )
    ),
    history_
    (
        IOobject
        (
            "refinementMRAHistory",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        List<refinementMRAHistory::splitCell8>(0),
        labelList(0),
        false
    ),
    faceRemover_(mesh_, GREAT),     // merge boundary faces wherever possible
    savedPointLevel_(0),
    savedCellLevel_(0),
    tetTesting_(false)
{
    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorInFunction
            << "Incorrect cellLevel or pointLevel size." << endl
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size() << endl
            << "Number of points in mesh:" << mesh_.nPoints()
            << " does not equal size of pointLevel:" << pointLevel_.size()
            << abort(FatalError);
    }

    // Check refinement levels for consistency
    checkRefinementLevels(-1, labelList(0));

    // Check initial mesh for consistency

    //if (debug)
    {
        checkMesh();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- creates new cells and calls the setHexRefinement and setTetRefinement
//  functions to set the mesh connectivity changes
//  All selected cells will be split into 8.
//  Returns per element in cells the 8 cells they were split into.
//  Guarantees that the 0th element is the original cell label.
Foam::labelListList Foam::hexTetRef8::setRefinement
(
    const labelList& cellsToRefine,
    polyTopoChange& meshMod
)
{
    if (debug)
    {
        Pout<< "hexTetRef8::setRefinement :"
            << " Checking initial mesh just to make sure" << endl;

        checkMesh();
        // Cannot call checkRefinementlevels since hanging points might
        // get triggered by the mesher after subsetting.
        //checkRefinementLevels(-1, labelList(0));
    }

    // Clear any saved point/cell data.
    savedPointLevel_.clear();
    savedCellLevel_.clear();

    // create new cell level and point level arrays
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // New point/cell level. Copy of pointLevel for existing points.
    DynamicList<label> newCellLevel(cellLevel_.size());
    forAll(cellLevel_, cellI)
    {
        newCellLevel.append(cellLevel_[cellI]);
    }
    DynamicList<label> newPointLevel(pointLevel_.size());
    forAll(pointLevel_, pointI)
    {
        newPointLevel.append(pointLevel_[pointI]);
    }

    // true  : cell split
    // false : cells not split
    boolList cellNeedsSplit(mesh_.nCells(), false);

    forAll(cellsToRefine, i)
    {
        label cellI = cellsToRefine[i];
        cellNeedsSplit[cellI] = true;
    }

    // Split edges
    // ~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexTetRef8::setRefinement :"
            << " Allocating edge midpoints."
            << endl;
    }

    // Unrefined edges are ones between cellLevel or lower points.
    // If any cell using this edge gets split then the edge needs to be split.

    // -1  : no need to split edge
    // >=0 : label of introduced mid point
    labelList newEdgeMids(mesh_.nEdges(), -1);
    // true : edge is split
    // false : edge is not split
    boolList edgeNeedsSplit(mesh_.nEdges(), false);

    // Note: Loop over cells to be refined or edges?
    forAll(cellNeedsSplit, cellI)
    {
        if (cellNeedsSplit[cellI])
        {
            const labelList& cEdges = mesh_.cellEdges(cellI);

            forAll(cEdges, i)
            {
                label edgeI = cEdges[i];

                const edge& e = mesh_.edges()[edgeI];

                if
                (
                    pointLevel_[e[0]] <= cellLevel_[cellI]
                 && pointLevel_[e[1]] <= cellLevel_[cellI]
                )
                {
                    edgeNeedsSplit[edgeI] = true;    // mark need for splitting
                }
            }
        }
    }

    // Synchronize edgeNeedsSplit across coupled patches. Take OR so that
    // any split takes precedence.
    syncTools::syncEdgeList
    (
        mesh_,
        edgeNeedsSplit,
        orEqOp<bool>(),
        true
    );


    // Introduce edge points
    // ~~~~~~~~~~~~~~~~~~~~~

    {
        // Phase 1: calculate midpoints and sync.
        // This needs doing for if people do not write binary and we slowly
        // get differences.

        pointField edgeMids(mesh_.nEdges(), point(-GREAT, -GREAT, -GREAT));

        forAll(edgeNeedsSplit, edgeI)
        {
            if (edgeNeedsSplit[edgeI])
            {
                // Edge marked to be split.
                edgeMids[edgeI] = mesh_.edges()[edgeI].centre(mesh_.points());
            }
        }
        syncTools::syncEdgePositions
        (
            mesh_,
            edgeMids,
            maxEqOp<vector>(),
            point(-GREAT, -GREAT, -GREAT)
        );


        // Phase 2: introduce points at the synced locations.
        forAll(edgeNeedsSplit, edgeI)
        {
            if (edgeNeedsSplit[edgeI])
            {
                // Edge marked to be split. Replace newEdgeMids with actual
                // point label.

                const edge& e = mesh_.edges()[edgeI];

                newEdgeMids[edgeI] = meshMod.setAction
                (
                    polyAddPoint
                    (
                        edgeMids[edgeI],            // point
                        e[0],                       // master point
                        -1,                         // zone for point
                        true                        // supports a cell
                    )
                );

                newPointLevel(newEdgeMids[edgeI]) =
                    max
                    (
                        pointLevel_[e[0]],
                        pointLevel_[e[1]]
                    )
                  + 1;
            }
        }
    }

    if (debug)
    {
        OFstream str(mesh_.time().path()/"newEdgeMids.obj");

        forAll(newEdgeMids, edgeI)
        {
            if (newEdgeMids[edgeI] >= 0)
            {
                const edge& e = mesh_.edges()[edgeI];

                meshTools::writeOBJ(str, e.centre(mesh_.points()));
            }
        }

        Pout<< "hexTetRef8::setRefinement :"
            << " Dumping edge centres to split to file " << str.name() << endl;
    }


    // Add the cells
    // ~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexTetRef8::setRefinement :"
            << " Adding cells."
            << endl;
    }

    // Per cell the 7 added cells (+ original cell)
    labelListList cellAddedCells(mesh_.nCells());

    forAll(cellsToRefine, i)
    {
        label cellI = cellsToRefine[i];

        labelList& cAdded = cellAddedCells[cellI];
        cAdded.setSize(8);

        // Original cell at 0
        cAdded[0] = cellI;

        for (label i = 1; i < 8; i++)
        {
            cAdded[i] = meshMod.setAction
            (
                polyAddCell
                (
                    -1,                                 // master point
                    -1,                                 // master edge
                    -1,                                 // master face
                    cellI,                              // master cell
                    mesh_.cellZones().whichZone(cellI)  // zone for cell
                )
            );

            newCellLevel(cAdded[i]) = cellLevel_[cellI]+1;
        }

        // Update cell level
        newCellLevel[cellI] = cellLevel_[cellI]+1;
    }


    // Find which faces need split
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexTetRef8::setRefinement :"
            << " Determining which faces to split."
            << endl;
    }

    // Face anchor level. There are guaranteed 3 OR 4 points with level
    // <= anchorLevel. These are the corner points.
    labelList faceAnchorLevel(mesh_.nFaces());

    for (label faceI = 0; faceI < mesh_.nFaces(); faceI++)
    {
        faceAnchorLevel[faceI] = faceLevel(faceI);
        if (faceAnchorLevel[faceI] == -1)
        {
            faceAnchorLevel[faceI] = faceLevelTet(faceI);
        }
    }

    // true : split face
    // false: don't split face
    boolList faceNeedsSplit(mesh_.nFaces(), false);

    // Internal faces: look at cells on both sides. Uniquely determined since
    // face itself guaranteed to be same level as most refined neighbour.
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        if (faceAnchorLevel[faceI] >= 0)
        {
            label own = mesh_.faceOwner()[faceI];
            label ownLevel = cellLevel_[own];
            label newOwnLevel = ownLevel + (cellNeedsSplit[own] ? 1 : 0);

            label nei = mesh_.faceNeighbour()[faceI];
            label neiLevel = cellLevel_[nei];
            label newNeiLevel = neiLevel + (cellNeedsSplit[nei] ? 1 : 0);

            if
            (
                newOwnLevel > faceAnchorLevel[faceI]
             || newNeiLevel > faceAnchorLevel[faceI]
            )
            {
                faceNeedsSplit[faceI] = true;  // mark to split
            }
        }
    }

    // Coupled patches handled like internal faces except now all information
    // from neighbour comes from across processor.
    // Boundary faces are more complicated since the boundary face can
    // be more refined than its owner (or neighbour for coupled patches)
    // (does not happen if refining/unrefining only, but does e.g. when
    //  refinining and subsetting)

    {
        labelList newNeiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(newNeiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
            label ownLevel = cellLevel_[own];
            label newOwnLevel = ownLevel + (cellNeedsSplit[own] ? 1 : 0);

            newNeiLevel[i] = newOwnLevel;
        }

        // Swap.
        syncTools::swapBoundaryFaceList(mesh_, newNeiLevel);

        // So now we have information on the neighbour.

        forAll(newNeiLevel, i)
        {
            label faceI = i+mesh_.nInternalFaces();

            if (faceAnchorLevel[faceI] >= 0)
            {
                label own = mesh_.faceOwner()[faceI];
                label ownLevel = cellLevel_[own];
                label newOwnLevel = ownLevel + (cellNeedsSplit[own] ? 1 : 0);

                if
                (
                    newOwnLevel > faceAnchorLevel[faceI]
                 || newNeiLevel[i] > faceAnchorLevel[faceI]
                )
                {

                    faceNeedsSplit[faceI] = true;  // mark to split
                }
            }
        }
    }


    // Synchronize faceNeedsSplit across coupled patches. (logical or)
    syncTools::syncFaceList
    (
        mesh_,
        faceNeedsSplit,
        orEqOp<bool>()
    );

    // Add edge splits to faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexTetRef8::setRefinement :"
            << " Adding edge splits to faces"
            << endl;
    }

    DynamicList<label> eFacesStorage;
    DynamicList<label> fEdgesStorage;

    forAll(newEdgeMids, edgeI)
    {
        if (newEdgeMids[edgeI] >= 0)
        {
            // Split edge. Check that face not already handled above.

            const labelList& eFaces = mesh_.edgeFaces(edgeI, eFacesStorage);

            forAll(eFaces, i)
            {
                label faceI = eFaces[i];

                const face& f = mesh_.faces()[faceI];
                const labelList& fEdges = mesh_.faceEdges
                (
                    faceI,
                    fEdgesStorage
                );

                DynamicList<label> newFaceVerts(f.size());

                forAll(f, fp)
                {
                    newFaceVerts.append(f[fp]);

                    label edgeJ = fEdges[fp];

                    if (newEdgeMids[edgeJ] >= 0)
                    {
                        newFaceVerts.append(newEdgeMids[edgeJ]);
                    }
                }

                face newFace;
                newFace.transfer(newFaceVerts);

                label own = mesh_.faceOwner()[faceI];
                label nei = -1;
                if (mesh_.isInternalFace(faceI))
                {
                    nei = mesh_.faceNeighbour()[faceI];
                }

                if (debug)
                {
                    if (mesh_.isInternalFace(faceI))
                    {
                        label oldOwn = mesh_.faceOwner()[faceI];
                        label oldNei = mesh_.faceNeighbour()[faceI];

                        checkInternalOrientation
                        (
                            meshMod,
                            oldOwn,
                            faceI,
                            mesh_.cellCentres()[oldOwn],
                            mesh_.cellCentres()[oldNei],
                            newFace
                        );
                    }
                    else
                    {
                        label oldOwn = mesh_.faceOwner()[faceI];

                        checkBoundaryOrientation
                        (
                            meshMod,
                            oldOwn,
                            faceI,
                            mesh_.cellCentres()[oldOwn],
                            mesh_.faceCentres()[faceI],
                            newFace
                        );
                    }
                }

                modFace(meshMod, faceI, newFace, own, nei);
            }
        }
    }

    // Separate the cell list into hex and tet cells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList hexCells;
    labelList tetCells;

    {
        DynamicList<label> hexCellsTemp;
        DynamicList<label> tetCellsTemp;

        forAll(cellsToRefine, i)
        {
            label cellI = cellsToRefine[i];

            label nAnchors = 0;

            const labelList &cPoints = mesh_.cellPoints()[cellI];

            forAll(cPoints, i)
            {
                label pointI = cPoints[i];

                if (pointLevel_[pointI] <= cellLevel_[cellI])
                {
                    nAnchors++;
                }
            }

            if (nAnchors == 4)
            {
                tetCellsTemp.append(cellI);
            }
            else if (nAnchors == 8)
            {
                hexCellsTemp.append(cellI);
            }
            else
            {
                FatalErrorInFunction
                    << "Cell " << cellI << " is neither hex nor tet. "
                    << "It has " << nAnchors << " anchors."
                    << abort(FatalError);
            }
        }

        hexCells.transfer(hexCellsTemp);
        tetCells.transfer(tetCellsTemp);
    }

    // Play refinement into the meshMod
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexTetRef8::setRefinement :"
            << " Feeding hex-refinements into the engine"
            << endl;
    }

    setHexRefinement
    (
        hexCells,
        faceNeedsSplit,
        newEdgeMids,
        cellAddedCells,
        newPointLevel,
        meshMod
    );

    if (debug)
    {
        Pout<< "hexTetRef8::setRefinement :"
            << " Feeding tet-refinements into the engine"
            << endl;
    }

    setTetRefinement
    (
        tetCells,
        cellNeedsSplit,
        faceNeedsSplit,
        newEdgeMids,
        cellAddedCells,
        newPointLevel,
        meshMod
    );

    pointLevel_.transfer(newPointLevel);
    cellLevel_.transfer(newCellLevel);

    // Mark files as changed
    setInstance(mesh_.facesInstance());


    // Update the live split cells tree.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // New unrefinement structure
    if (history_.active())
    {
        if (debug)
        {
            Pout<< "hexTetRef8::setRefinement :"
                << " Updating refinement history to " << cellLevel_.size()
                << " cells" << endl;
        }

        // Extend refinement history for new cells
        history_.resize(cellLevel_.size());

        forAll(cellAddedCells, cellI)
        {
            const labelList& addedCells = cellAddedCells[cellI];

            if (addedCells.size())
            {
                // Cell was split.
                history_.storeSplit(cellI, addedCells);
            }
        }
    }

    // Compact cellAddedCells.

    labelListList refinedCells(hexCells.size()+tetCells.size());

    forAll(hexCells, i)
    {
        label cellI = hexCells[i];

        refinedCells[i].transfer(cellAddedCells[cellI]);
    }
    forAll(tetCells, i)
    {
        label cellI = tetCells[i];

        refinedCells[i].transfer(cellAddedCells[cellI]);
    }

    return refinedCells;
}


Foam::labelList Foam::hexTetRef8::consistentRefinement
(
    const labelList& cellsToRefine,
    const bool maxSet
) const
{
    // Loop, modifying cellsToRefine, until no more changes to due to 2:1
    // conflicts.
    // maxSet = false : unselect cells to refine
    // maxSet = true  : select cells to refine

    // Go to straight boolList.
    PackedBoolList refineCell(mesh_.nCells());
    forAll(cellsToRefine, i)
    {
        refineCell.set(cellsToRefine[i]);
    }

    while (true)
    {
        label nChanged = faceConsistentRefinement(maxSet, refineCell);

        reduce(nChanged, sumOp<label>());

        if (debug)
        {
            Pout<< "hexTetRef8::consistentRefinement : Changed " << nChanged
                << " refinement levels due to 2:1 conflicts."
                << endl;
        }

        if (nChanged == 0)
        {
            break;
        }
    }


    // Convert back to labelList.
    label nRefined = 0;

    forAll(refineCell, cellI)
    {
        if (refineCell.get(cellI))
        {
            nRefined++;
        }
    }

    labelList newCellsToRefine(nRefined);
    nRefined = 0;

    forAll(refineCell, cellI)
    {
        if (refineCell.get(cellI))
        {
            newCellsToRefine[nRefined++] = cellI;
        }
    }

    if (debug)
    {
        checkWantedRefinementLevels(newCellsToRefine);
    }

    return newCellsToRefine;
}


// Given a list of cells to refine determine additional cells to refine
// such that the overall refinement:
// - satisfies maxFaceDiff (e.g. 2:1) across neighbouring faces
// - satisfies maxPointDiff (e.g. 4:1) across selected point connected
//   cells. This is used to ensure that e.g. cells on the surface are not
//   point connected to cells which are 8 times smaller.
Foam::labelList Foam::hexTetRef8::consistentSlowRefinement
(
    const label maxFaceDiff,
    const labelList& cellsToRefine,
    const labelList& facesToCheck,
    const label maxPointDiff,
    const labelList& pointsToCheck
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();


    if (maxFaceDiff <= 0)
    {
        FatalErrorInFunction
            << "Illegal maxFaceDiff " << maxFaceDiff << nl
            << "Value should be >= 1" << exit(FatalError);
    }


    // Bit tricky. Say we want a distance of three cells between two
    // consecutive refinement levels. This is done by using FaceCellWave to
    // transport out the new refinement level. It gets decremented by one
    // every cell it crosses so if we initialize it to maxFaceDiff
    // we will get a field everywhere that tells us whether an unselected cell
    // needs refining as well.


    // Initial information about (distance to) cellLevel on all cells
    List<refinementData> allCellInfo(mesh_.nCells());

    // Initial information about (distance to) cellLevel on all faces
    List<refinementData> allFaceInfo(mesh_.nFaces());

    forAll(allCellInfo, cellI)
    {
        // maxFaceDiff since refinementData counts both
        // faces and cells.
        allCellInfo[cellI] = refinementData
        (
            maxFaceDiff*(cellLevel_[cellI]+1),// when cell is to be refined
            maxFaceDiff*cellLevel_[cellI]     // current level
        );
    }

    // Cells to be refined will have cellLevel+1
    forAll(cellsToRefine, i)
    {
        label cellI = cellsToRefine[i];

        allCellInfo[cellI].count() = allCellInfo[cellI].refinementCount();
    }


    // Labels of seed faces
    DynamicList<label> seedFaces(mesh_.nFaces()/100);
    // refinementLevel data on seed faces
    DynamicList<refinementData> seedFacesInfo(mesh_.nFaces()/100);

    // Dummy additional info for FaceCellWave
    int dummyTrackData = 0;


    // Additional buffer layer thickness by changing initial count. Usually
    // this happens on boundary faces. Bit tricky. Use allFaceInfo to mark
    // off thus marked faces so they're skipped in the next loop.
    forAll(facesToCheck, i)
    {
        label faceI = facesToCheck[i];

        if (allFaceInfo[faceI].valid(dummyTrackData))
        {
            // Can only occur if face has already gone through loop below.
            FatalErrorInFunction
                << "Argument facesToCheck seems to have duplicate entries!"
                << endl
                << "face:" << faceI << " occurs at positions "
                << findIndices(facesToCheck, faceI)
                << abort(FatalError);
        }


        const refinementData& ownData = allCellInfo[faceOwner[faceI]];

        if (mesh_.isInternalFace(faceI))
        {
            // Seed face if neighbouring cell (after possible refinement)
            // will be refined one more than the current owner or neighbour.

            const refinementData& neiData = allCellInfo[faceNeighbour[faceI]];

            label faceCount;
            label faceRefineCount;
            if (neiData.count() > ownData.count())
            {
                faceCount = neiData.count() + maxFaceDiff;
                faceRefineCount = faceCount + maxFaceDiff;
            }
            else
            {
                faceCount = ownData.count() + maxFaceDiff;
                faceRefineCount = faceCount + maxFaceDiff;
            }

            seedFaces.append(faceI);
            seedFacesInfo.append
            (
                refinementData
                (
                    faceRefineCount,
                    faceCount
                )
            );
            allFaceInfo[faceI] = seedFacesInfo.last();
        }
        else
        {
            label faceCount = ownData.count() + maxFaceDiff;
            label faceRefineCount = faceCount + maxFaceDiff;

            seedFaces.append(faceI);
            seedFacesInfo.append
            (
                refinementData
                (
                    faceRefineCount,
                    faceCount
                )
            );
            allFaceInfo[faceI] = seedFacesInfo.last();
        }
    }


    // Just seed with all faces inbetween different refinement levels for now
    // (alternatively only seed faces on cellsToRefine but that gives problems
    //  if no cells to refine)
    forAll(faceNeighbour, faceI)
    {
        // Check if face already handled in loop above
        if (!allFaceInfo[faceI].valid(dummyTrackData))
        {
            label own = faceOwner[faceI];
            label nei = faceNeighbour[faceI];

            // Seed face with transported data from highest cell.

            if (allCellInfo[own].count() > allCellInfo[nei].count())
            {
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    own,
                    allCellInfo[own],
                    FaceCellWave<refinementData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFaces.append(faceI);
                seedFacesInfo.append(allFaceInfo[faceI]);
            }
            else if (allCellInfo[own].count() < allCellInfo[nei].count())
            {
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    nei,
                    allCellInfo[nei],
                    FaceCellWave<refinementData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFaces.append(faceI);
                seedFacesInfo.append(allFaceInfo[faceI]);
            }
        }
    }

    // Seed all boundary faces with owner value. This is to make sure that
    // they are visited (probably only important for coupled faces since
    // these need to be visited from both sides)
    for (label faceI = mesh_.nInternalFaces(); faceI < mesh_.nFaces(); faceI++)
    {
        // Check if face already handled in loop above
        if (!allFaceInfo[faceI].valid(dummyTrackData))
        {
            label own = faceOwner[faceI];

            // Seed face with transported data from owner.
            refinementData faceData;
            faceData.updateFace
            (
                mesh_,
                faceI,
                own,
                allCellInfo[own],
                FaceCellWave<refinementData, int>::propagationTol(),
                dummyTrackData
            );
            seedFaces.append(faceI);
            seedFacesInfo.append(faceData);
        }
    }


    // face-cell-face transport engine
    FaceCellWave<refinementData, int> levelCalc
    (
        mesh_,
        allFaceInfo,
        allCellInfo,
        dummyTrackData
    );

    while (true)
    {
        if (debug)
        {
            Pout<< "hexTetRef8::consistentSlowRefinement : Seeded "
                << seedFaces.size() << " faces between cells with different"
                << " refinement level." << endl;
        }

        // Set seed faces
        levelCalc.setFaceInfo(seedFaces.shrink(), seedFacesInfo.shrink());
        seedFaces.clear();
        seedFacesInfo.clear();

        // Iterate until no change. Now 2:1 face difference should be satisfied
        levelCalc.iterate(mesh_.globalData().nTotalFaces()+1);


        // Now check point-connected cells (face-connected cells already ok):
        // - get per point max of connected cells
        // - sync across coupled points
        // - check cells against above point max

        if (maxPointDiff == -1)
        {
            // No need to do any point checking.
            break;
        }

        // Determine per point the max cell level. (done as count, not
        // as cell level purely for ease)
        labelList maxPointCount(mesh_.nPoints(), 0);

        forAll(maxPointCount, pointI)
        {
            label& pLevel = maxPointCount[pointI];

            const labelList& pCells = mesh_.pointCells(pointI);

            forAll(pCells, i)
            {
                pLevel = max(pLevel, allCellInfo[pCells[i]].count());
            }
        }

        // Sync maxPointCount to neighbour
        syncTools::syncPointList
        (
            mesh_,
            maxPointCount,
            maxEqOp<label>(),
            labelMin            // null value
        );

        // Update allFaceInfo from maxPointCount for all points to check
        // (usually on boundary faces)

        // Per face the new refinement data
        Map<refinementData> changedFacesInfo(pointsToCheck.size());

        forAll(pointsToCheck, i)
        {
            label pointI = pointsToCheck[i];

            // Loop over all cells using the point and check whether their
            // refinement level is much less than the maximum.

            const labelList& pCells = mesh_.pointCells(pointI);

            forAll(pCells, pCellI)
            {
                label cellI = pCells[pCellI];

                refinementData& cellInfo = allCellInfo[cellI];

                if
                (
                   !cellInfo.isRefined()
                 && (
                        maxPointCount[pointI]
                      > cellInfo.count() + maxFaceDiff*maxPointDiff
                    )
                )
                {
                    // Mark cell for refinement
                    cellInfo.count() = cellInfo.refinementCount();

                    // Insert faces of cell as seed faces.
                    const cell& cFaces = mesh_.cells()[cellI];

                    forAll(cFaces, cFaceI)
                    {
                        label faceI = cFaces[cFaceI];

                        refinementData faceData;
                        faceData.updateFace
                        (
                            mesh_,
                            faceI,
                            cellI,
                            cellInfo,
                            FaceCellWave<refinementData, int>::propagationTol(),
                            dummyTrackData
                        );

                        if (faceData.count() > allFaceInfo[faceI].count())
                        {
                            changedFacesInfo.insert(faceI, faceData);
                        }
                    }
                }
            }
        }

        label nChanged = changedFacesInfo.size();
        reduce(nChanged, sumOp<label>());

        if (nChanged == 0)
        {
            break;
        }


        // Transfer into seedFaces, seedFacesInfo
        seedFaces.setCapacity(changedFacesInfo.size());
        seedFacesInfo.setCapacity(changedFacesInfo.size());

        forAllConstIter(Map<refinementData>, changedFacesInfo, iter)
        {
            seedFaces.append(iter.key());
            seedFacesInfo.append(iter());
        }
    }


    if (debug)
    {
        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            label own = mesh_.faceOwner()[faceI];
            label ownLevel =
                cellLevel_[own]
              + (allCellInfo[own].isRefined() ? 1 : 0);

            label nei = mesh_.faceNeighbour()[faceI];
            label neiLevel =
                cellLevel_[nei]
              + (allCellInfo[nei].isRefined() ? 1 : 0);

            if (mag(ownLevel-neiLevel) > 1)
            {
                dumpCell(own);
                dumpCell(nei);
                FatalErrorInFunction
                    << "cell:" << own
                    << " current level:" << cellLevel_[own]
                    << " current refData:" << allCellInfo[own]
                    << " level after refinement:" << ownLevel
                    << nl
                    << "neighbour cell:" << nei
                    << " current level:" << cellLevel_[nei]
                    << " current refData:" << allCellInfo[nei]
                    << " level after refinement:" << neiLevel
                    << nl
                    << "which does not satisfy 2:1 constraints anymore." << nl
                    << "face:" << faceI << " faceRefData:" << allFaceInfo[faceI]
                    << abort(FatalError);
            }
        }


        // Coupled faces. Swap owner level to get neighbouring cell level.
        // (only boundary faces of neiLevel used)

        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
        labelList neiCount(mesh_.nFaces()-mesh_.nInternalFaces());
        labelList neiRefCount(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(neiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
            neiLevel[i] = cellLevel_[own];
            neiCount[i] = allCellInfo[own].count();
            neiRefCount[i] = allCellInfo[own].refinementCount();
        }

        // Swap to neighbour
        syncTools::swapBoundaryFaceList(mesh_, neiLevel);
        syncTools::swapBoundaryFaceList(mesh_, neiCount);
        syncTools::swapBoundaryFaceList(mesh_, neiRefCount);

        // Now we have neighbour value see which cells need refinement
        forAll(neiLevel, i)
        {
            label faceI = i+mesh_.nInternalFaces();

            label own = mesh_.faceOwner()[faceI];
            label ownLevel =
                cellLevel_[own]
              + (allCellInfo[own].isRefined() ? 1 : 0);

            label nbrLevel =
                neiLevel[i]
              + ((neiCount[i] >= neiRefCount[i]) ? 1 : 0);

            if (mag(ownLevel - nbrLevel) > 1)
            {
                dumpCell(own);
                label patchI = mesh_.boundaryMesh().whichPatch(faceI);

                FatalErrorInFunction
                    << "Celllevel does not satisfy 2:1 constraint."
                    << " On coupled face "
                    << faceI
                    << " refData:" << allFaceInfo[faceI]
                    << " on patch " << patchI << " "
                    << mesh_.boundaryMesh()[patchI].name() << nl
                    << "owner cell " << own
                    << " current level:" << cellLevel_[own]
                    << " current count:" << allCellInfo[own].count()
                    << " current refCount:"
                    << allCellInfo[own].refinementCount()
                    << " level after refinement:" << ownLevel
                    << nl
                    << "(coupled) neighbour cell"
                    << " has current level:" << neiLevel[i]
                    << " current count:" << neiCount[i]
                    << " current refCount:" << neiRefCount[i]
                    << " level after refinement:" << nbrLevel
                    << abort(FatalError);
            }
        }
    }

    // Convert back to labelList of cells to refine.

    label nRefined = 0;

    forAll(allCellInfo, cellI)
    {
        if (allCellInfo[cellI].isRefined())
        {
            nRefined++;
        }
    }

    // Updated list of cells to refine
    labelList newCellsToRefine(nRefined);
    nRefined = 0;

    forAll(allCellInfo, cellI)
    {
        if (allCellInfo[cellI].isRefined())
        {
            newCellsToRefine[nRefined++] = cellI;
        }
    }

    if (debug)
    {
        Pout<< "hexTetRef8::consistentSlowRefinement : From "
            << cellsToRefine.size() << " to " << newCellsToRefine.size()
            << " cells to refine." << endl;
    }

    return newCellsToRefine;
}


Foam::labelList Foam::hexTetRef8::consistentSlowRefinement2
(
    const label maxFaceDiff,
    const labelList& cellsToRefine,
    const labelList& facesToCheck
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    if (maxFaceDiff <= 0)
    {
        FatalErrorInFunction
            << "Illegal maxFaceDiff " << maxFaceDiff << nl
            << "Value should be >= 1" << exit(FatalError);
    }

    const scalar level0Size = 2*maxFaceDiff*level0EdgeLength();


    // Bit tricky. Say we want a distance of three cells between two
    // consecutive refinement levels. This is done by using FaceCellWave to
    // transport out the 'refinement shell'. Anything inside the refinement
    // shell (given by a distance) gets marked for refinement.

    // Initial information about (distance to) cellLevel on all cells
    List<refinementDistanceData> allCellInfo(mesh_.nCells());

    // Initial information about (distance to) cellLevel on all faces
    List<refinementDistanceData> allFaceInfo(mesh_.nFaces());

    // Dummy additional info for FaceCellWave
    int dummyTrackData = 0;


    // Mark cells with wanted refinement level
    forAll(cellsToRefine, i)
    {
        label cellI = cellsToRefine[i];

        allCellInfo[cellI] = refinementDistanceData
        (
            level0Size,
            mesh_.cellCentres()[cellI],
            cellLevel_[cellI]+1             // wanted refinement
        );
    }
    // Mark all others with existing refinement level
    forAll(allCellInfo, cellI)
    {
        if (!allCellInfo[cellI].valid(dummyTrackData))
        {
            allCellInfo[cellI] = refinementDistanceData
            (
                level0Size,
                mesh_.cellCentres()[cellI],
                cellLevel_[cellI]           // wanted refinement
            );
        }
    }


    // Labels of seed faces
    DynamicList<label> seedFaces(mesh_.nFaces()/100);
    // refinementLevel data on seed faces
    DynamicList<refinementDistanceData> seedFacesInfo(mesh_.nFaces()/100);

    const pointField& cc = mesh_.cellCentres();

    forAll(facesToCheck, i)
    {
        label faceI = facesToCheck[i];

        if (allFaceInfo[faceI].valid(dummyTrackData))
        {
            // Can only occur if face has already gone through loop below.
            FatalErrorInFunction
                << "Argument facesToCheck seems to have duplicate entries!"
                << endl
                << "face:" << faceI << " occurs at positions "
                << findIndices(facesToCheck, faceI)
                << abort(FatalError);
        }

        label own = faceOwner[faceI];

        label ownLevel =
        (
            allCellInfo[own].valid(dummyTrackData)
          ? allCellInfo[own].originLevel()
          : cellLevel_[own]
        );

        if (!mesh_.isInternalFace(faceI))
        {
            // Do as if boundary face would have neighbour with one higher
            // refinement level.
            const point& fc = mesh_.faceCentres()[faceI];

            refinementDistanceData neiData
            (
                level0Size,
                2*fc - cc[own],    // est'd cell centre
                ownLevel+1
            );

            allFaceInfo[faceI].updateFace
            (
                mesh_,
                faceI,
                own,        // not used (should be nei)
                neiData,
                FaceCellWave<refinementDistanceData, int>::propagationTol(),
                dummyTrackData
            );
        }
        else
        {
            label nei = faceNeighbour[faceI];

            label neiLevel =
            (
                allCellInfo[nei].valid(dummyTrackData)
              ? allCellInfo[nei].originLevel()
              : cellLevel_[nei]
            );

            if (ownLevel == neiLevel)
            {
                // Fake as if nei>own or own>nei (whichever one 'wins')
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    nei,
                    refinementDistanceData(level0Size, cc[nei], neiLevel+1),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    own,
                    refinementDistanceData(level0Size, cc[own], ownLevel+1),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
            }
            else
            {
                // Difference in level anyway.
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    nei,
                    refinementDistanceData(level0Size, cc[nei], neiLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    own,
                    refinementDistanceData(level0Size, cc[own], ownLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
            }
        }
        seedFaces.append(faceI);
        seedFacesInfo.append(allFaceInfo[faceI]);
    }


    // Create some initial seeds to start walking from. This is only if there
    // are no facesToCheck.
    // Just seed with all faces inbetween different refinement levels for now
    forAll(faceNeighbour, faceI)
    {
        // Check if face already handled in loop above
        if (!allFaceInfo[faceI].valid(dummyTrackData))
        {
            label own = faceOwner[faceI];

            label ownLevel =
            (
                allCellInfo[own].valid(dummyTrackData)
              ? allCellInfo[own].originLevel()
              : cellLevel_[own]
            );

            label nei = faceNeighbour[faceI];

            label neiLevel =
            (
                allCellInfo[nei].valid(dummyTrackData)
              ? allCellInfo[nei].originLevel()
              : cellLevel_[nei]
            );

            if (ownLevel > neiLevel)
            {
                // Set face to owner data. (since face not yet would be copy)
                seedFaces.append(faceI);
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    own,
                    refinementDistanceData(level0Size, cc[own], ownLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFacesInfo.append(allFaceInfo[faceI]);
            }
            else if (neiLevel > ownLevel)
            {
                seedFaces.append(faceI);
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    nei,
                    refinementDistanceData(level0Size, cc[nei], neiLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFacesInfo.append(allFaceInfo[faceI]);
            }
        }
    }

    seedFaces.shrink();
    seedFacesInfo.shrink();

    // face-cell-face transport engine
    FaceCellWave<refinementDistanceData, int> levelCalc
    (
        mesh_,
        seedFaces,
        seedFacesInfo,
        allFaceInfo,
        allCellInfo,
        mesh_.globalData().nTotalCells()+1,
        dummyTrackData
    );


    //if (debug)
    //{
    //    // Dump wanted level
    //    volScalarField wantedLevel
    //    (
    //        IOobject
    //        (
    //            "wantedLevel",
    //            fMesh.time().timeName(),
    //            fMesh,
    //            IOobject::NO_READ,
    //            IOobject::AUTO_WRITE,
    //            false
    //        ),
    //        fMesh,
    //        dimensionedScalar("zero", dimless, 0)
    //    );
    //
    //    forAll(wantedLevel, cellI)
    //    {
    //        wantedLevel[cellI] = allCellInfo[cellI].wantedLevel(cc[cellI]);
    //    }
    //
    //    Pout<< "Writing " << wantedLevel.objectPath() << endl;
    //    wantedLevel.write();
    //}


    // Convert back to labelList of cells to refine.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // 1. Force original refinement cells to be picked up by setting the
    // originLevel of input cells to be a very large level (but within range
    // of 1<< shift inside refinementDistanceData::wantedLevel)
    forAll(cellsToRefine, i)
    {
        label cellI = cellsToRefine[i];

        allCellInfo[cellI].originLevel() = sizeof(label)*8-2;
        allCellInfo[cellI].origin() = cc[cellI];
    }

    // 2. Extend to 2:1. I don't understand yet why this is not done
    // 2. Extend to 2:1. For non-cube cells the scalar distance does not work
    // so make sure it at least provides 2:1.
    PackedBoolList refineCell(mesh_.nCells());
    forAll(allCellInfo, cellI)
    {
        label wanted = allCellInfo[cellI].wantedLevel(cc[cellI]);

        if (wanted > cellLevel_[cellI]+1)
        {
            refineCell.set(cellI);
        }
    }
    faceConsistentRefinement(true, refineCell);

    while (true)
    {
        label nChanged = faceConsistentRefinement(true, refineCell);

        reduce(nChanged, sumOp<label>());

        if (debug)
        {
            Pout<< "hexTetRef8::consistentSlowRefinement2 : Changed " << nChanged
                << " refinement levels due to 2:1 conflicts."
                << endl;
        }

        if (nChanged == 0)
        {
            break;
        }
    }

    // 3. Convert back to labelList.
    label nRefined = 0;

    forAll(refineCell, cellI)
    {
//        if (refineCell.get(cellI))
        if (refineCell[cellI])
        {
            nRefined++;
        }
    }

    labelList newCellsToRefine(nRefined);
    nRefined = 0;

    forAll(refineCell, cellI)
    {
//        if (refineCell.get(cellI))
        if (refineCell[cellI])
        {
            newCellsToRefine[nRefined++] = cellI;
        }
    }

    if (debug)
    {
        Pout<< "hexTetRef8::consistentSlowRefinement2 : From "
            << cellsToRefine.size() << " to " << newCellsToRefine.size()
            << " cells to refine." << endl;

        // Check that newCellsToRefine obeys at least 2:1.

        {
            cellSet cellsIn(mesh_, "cellsToRefineIn", cellsToRefine);
            Pout<< "hexTetRef8::consistentSlowRefinement2 : writing "
                << cellsIn.size() << " to cellSet "
                << cellsIn.objectPath() << endl;
            cellsIn.write();
        }
        {
            cellSet cellsOut(mesh_, "cellsToRefineOut", newCellsToRefine);
            Pout<< "hexTetRef8::consistentSlowRefinement2 : writing "
                << cellsOut.size() << " to cellSet "
                << cellsOut.objectPath() << endl;
            cellsOut.write();
        }

        // Extend to 2:1
        PackedBoolList refineCell(mesh_.nCells());
        forAll(newCellsToRefine, i)
        {
            refineCell.set(newCellsToRefine[i]);
        }
        const PackedBoolList savedRefineCell(refineCell);

        label nChanged = faceConsistentRefinement(true, refineCell);

        {
            cellSet cellsOut2
            (
                mesh_, "cellsToRefineOut2", newCellsToRefine.size()
            );
            forAll(refineCell, cellI)
            {
                if (refineCell.get(cellI))
                {
                    cellsOut2.insert(cellI);
                }
            }
            Pout<< "hexTetRef8::consistentSlowRefinement2 : writing "
                << cellsOut2.size() << " to cellSet "
                << cellsOut2.objectPath() << endl;
            cellsOut2.write();
        }

        if (nChanged > 0)
        {
            forAll(refineCell, cellI)
            {
                if (refineCell.get(cellI) && !savedRefineCell.get(cellI))
                {
                    dumpCell(cellI);
                    FatalErrorInFunction
                        << "Cell:" << cellI << " cc:"
                        << mesh_.cellCentres()[cellI]
                        << " was not marked for refinement but does not obey"
                        << " 2:1 constraints."
                        << abort(FatalError);
                }
            }
        }
    }

    return newCellsToRefine;
}


// Top level driver to insert topo changes to do all refinement.
void Foam::hexTetRef8::setHexRefinement
(
    const labelList& hexCellLabels,
    boolList& faceNeedsSplit,
    const labelList& newEdgeMids,
    const labelListList& cellAddedCells,
    DynamicList<label>& newPointLevel,
    polyTopoChange& meshMod
)
{
    if (debug)
    {
        Pout<< "hexTetRef8::setHexRefinement :"
            << " Allocating " << hexCellLabels.size() << " cell midpoints."
            << endl;
    }


    // Mid point per hex-refined cell.
    // -1  : not hex-refined
    // >=0 : label of mid point.
    labelList cellMidPoint(mesh_.nCells(), -1);

    forAll(hexCellLabels, i)
    {
        label cellI = hexCellLabels[i];

        label anchorPointI = mesh_.faces()[mesh_.cells()[cellI][0]][0];

        cellMidPoint[cellI] = meshMod.setAction
        (
            polyAddPoint
            (
                mesh_.cellCentres()[cellI],     // point
                anchorPointI,                   // master point
                -1,                             // zone for point
                true                            // supports a cell
            )
        );

        newPointLevel(cellMidPoint[cellI]) = cellLevel_[cellI]+1;
    }


    if (debug)
    {
        cellSet splitCells(mesh_, "splitCells", hexCellLabels.size());

        forAll(cellMidPoint, cellI)
        {
            if (cellMidPoint[cellI] >= 0)
            {
                splitCells.insert(cellI);
            }
        }

        Pout<< "hexTetRef8::setHexRefinement : Dumping " << splitCells.size()
            << " cells to split to cellSet " << splitCells.objectPath()
            << endl;

        splitCells.write();
    }



    // get split edges mid points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexTetRef8::setHexRefinement :"
            << " Getting edge midpoints."
            << endl;
    }

    // Face anchor level. There are guaranteed 3 OR 4 points with level
    // <= anchorLevel. These are the corner points.
    labelList faceAnchorLevel(mesh_.nFaces());

    // Number of face anchors (points of level <= faceAnchorLevel)
    labelList nFaceAnchors(mesh_.nFaces());

    for (label faceI = 0; faceI < mesh_.nFaces(); faceI++)
    {
        faceAnchorLevel[faceI] = faceLevel(faceI);
        if (faceAnchorLevel[faceI] == -1)
        {
            faceAnchorLevel[faceI] = faceLevelTet(faceI);
        }

        nFaceAnchors[faceI]
            = countAnchors(mesh_.faces()[faceI], faceAnchorLevel[faceI]);
    }

    // -1  : no need to split edge
    // >=0 : label of introduced mid point (calc in setRefinement)
    labelList edgeMidPoint(mesh_.nEdges(), -1);

    const faceList& cFaces = mesh_.faces();

    forAll(cFaces, faceI)
    {
        DynamicList<label> fEdgesStorage;
        const labelList& fEdges = mesh_.faceEdges
        (
            faceI,
            fEdgesStorage
        );

        forAll(fEdges, j)
        {
            label edgeJ = fEdges[j];

            if (newEdgeMids[edgeJ] >= 0)
            {
                edgeMidPoint[edgeJ] = newEdgeMids[edgeJ];
            }
        }
    }

    if (debug)
    {
        OFstream str(mesh_.time().path()/"edgeMidPoint.obj");

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                const edge& e = mesh_.edges()[edgeI];

                meshTools::writeOBJ(str, e.centre(mesh_.points()));
            }
        }

        Pout<< "hexTetRef8::setHexRefinement :"
            << " Dumping edge centres to split to file "
            << str.name() << endl;
    }


    // Allocate face mids
    // ~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexTetRef8::setHexRefinement :"
            << " Allocating face midpoints."
            << endl;
    }


    // -1  : no need to add mid-point to face
    // >=0 : label of introduced mid point
    labelList faceMidPoint(mesh_.nFaces(), -1);

    for (label faceI = 0; faceI < mesh_.nFaces(); faceI++)
    {
        if (faceNeedsSplit[faceI])
        {
            // we only add a mid-point for square faces
            if (nFaceAnchors[faceI] == 4)
            {
                faceMidPoint[faceI] = 12345;    // mark to add mid-point
            }
        }
    }

    // Synchronize faceMidPoint across coupled patches. (logical or)
    syncTools::syncFaceList
    (
        mesh_,
        faceMidPoint,
        maxEqOp<label>()
    );

    // Introduce face points
    // ~~~~~~~~~~~~~~~~~~~~~

    {
        // Phase 1: determine mid points and sync. See comment for edgeMids
        // above
        pointField bFaceMids
        (
            mesh_.nFaces()-mesh_.nInternalFaces(),
            point(-GREAT, -GREAT, -GREAT)
        );

        forAll(bFaceMids, i)
        {
            label faceI = i+mesh_.nInternalFaces();

            if (faceMidPoint[faceI] >= 0)
            {
                bFaceMids[i] = mesh_.faceCentres()[faceI];
            }
        }
        syncTools::syncBoundaryFacePositions
        (
            mesh_,
            bFaceMids,
            maxEqOp<vector>()
        );

        forAll(faceMidPoint, faceI)
        {
            if (faceMidPoint[faceI] >= 0)
            {
                // Face marked to be split. Replace faceMidPoint with actual
                // point label.

                const face& f = mesh_.faces()[faceI];

                faceMidPoint[faceI] = meshMod.setAction
                (
                    polyAddPoint
                    (
                        (
                            faceI < mesh_.nInternalFaces()
                          ? mesh_.faceCentres()[faceI]
                          : bFaceMids[faceI-mesh_.nInternalFaces()]
                        ),                          // point
                        f[0],                       // master point
                        -1,                         // zone for point
                        true                        // supports a cell
                    )
                );

                // Determine the level of the corner points and midpoint will
                // be one higher.
                newPointLevel(faceMidPoint[faceI]) = faceAnchorLevel[faceI]+1;
            }
        }
    }

    if (debug)
    {
        faceSet splitFaces(mesh_, "splitFaces", hexCellLabels.size());

        forAll(faceNeedsSplit, faceI)
        {
            if (faceNeedsSplit[faceI])
            {
                splitFaces.insert(faceI);
            }
        }

        Pout<< "hexTetRef8::setHexRefinement : Dumping " << splitFaces.size()
            << " faces to split to faceSet " << splitFaces.objectPath() << endl;

        splitFaces.write();
    }


    // Information complete
    // ~~~~~~~~~~~~~~~~~~~~
    // At this point we have all the information we need. We should no
    // longer reference the cellLabels to refine. All the information is:
    // - cells:
    //      -> cellMidPoint >= 0 : cell needs to be split here
    // - faces:
    //      -> isFaceSplit: true : face needs to be split
    //      -> faceMidPoint >= 0 : face needs to be hex-split
    // - edges:
    //      -> edgeMidPoint >= 0 : edge needs to be split



    // Get the corner/anchor points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexTetRef8::setHexRefinement :"
            << " Finding cell anchorPoints (8 per cell)"
            << endl;
    }

    // There will always be 8 points on the hex that have were introduced
    // with the hex and will have the same or lower refinement level.

    // Per cell the 8 corner points.
    labelListList cellAnchorPoints(mesh_.nCells());

    {
        labelList nAnchorPoints(mesh_.nCells(), 0);

        forAll(cellMidPoint, cellI)
        {
            if (cellMidPoint[cellI] >= 0)
            {
                cellAnchorPoints[cellI].setSize(8);

                const labelList &cPoints = mesh_.cellPoints()[cellI];

                forAll(cPoints, j)
                {
                    label pJ = cPoints[j];

                    if
                    (
                        pointLevel_[pJ] <= cellLevel_[cellI]
                    )
                    {
                        if (nAnchorPoints[cellI] == 8)
                        {
                            dumpCell(cellI);
                            FatalErrorInFunction
                                << "cell " << cellI
                                << " of level " << cellLevel_[cellI]
                                << " uses more than 8 points of equal or"
                                << " lower level" << nl
                                << "Points so far:" << cellAnchorPoints[cellI]
                                << abort(FatalError);
                        }

                        cellAnchorPoints[cellI][nAnchorPoints[cellI]++] = pJ;
                    }
                }
            }
        }

        forAll(cellMidPoint, cellI)
        {
            if (cellMidPoint[cellI] >= 0)
            {
                if (nAnchorPoints[cellI] != 8)
                {
                    dumpCell(cellI);

                    const labelList& cPoints = mesh_.cellPoints(cellI);

                    FatalErrorInFunction
                        << "cell " << cellI
                        << " of level " << cellLevel_[cellI]
                        << " does not seem to have 8 points of equal or"
                        << " lower level" << endl
                        << "cellPoints:" << cPoints << endl
                        << "pointLevels:"
                        << UIndirectList<label>(pointLevel_, cPoints)() << endl
                        << abort(FatalError);
                }
            }
        }
    }


    // Faces
    // ~~~~~
    // 1. existing faces that get split (into four always)
    // 2. existing faces that do not get split but get new owner/neighbour
    // 3. new internal faces inside split cells.

    if (debug)
    {
        Pout<< "hexTetRef8::setHexRefinement :"
            << " Marking faces to be handled"
            << endl;
    }

    // Get all affected faces.
    PackedBoolList affectedFace(mesh_.nFaces());

    {
        forAll(cellMidPoint, cellI)
        {
            if (cellMidPoint[cellI] >= 0)
            {
                const cell& cFaces = mesh_.cells()[cellI];

                forAll(cFaces, i)
                {
                    affectedFace.set(cFaces[i]);
                }
            }
        }

        forAll(faceMidPoint, faceI)
        {
            if (faceMidPoint[faceI] >= 0)
            {
                affectedFace.set(faceI);
            }
        }
    }


    // 1. Faces that get split
    // ~~~~~~~~~~~~~~~~~~~~~~~

    // a. normal square faces that get split the old way
    // b. tri-faces that get split the new way

    if (debug)
    {
        Pout<< "hexTetRef8::setRefinement : Splitting faces" << endl;
    }

    forAll(faceMidPoint, faceI)
    {
        if (faceMidPoint[faceI] >= 0 && affectedFace.get(faceI))
        {
            // Face needs to be split and hasn't yet been done in some way
            // (affectedFace - is impossible since this is first change but
            //  just for completeness)

            const face& f = mesh_.faces()[faceI];

            // Has original faceI been used (three faces added, original gets
            // modified)
            bool modifiedFace = false;

            label anchorLevel = faceAnchorLevel[faceI];

            face newFace(4);

            forAll(f, fp)
            {
                label pointI = f[fp];

                if (pointLevel_[pointI] <= anchorLevel)
                {
                    // point is anchor. Start collecting face.

                    DynamicList<label> faceVerts(4);

                    faceVerts.append(pointI);

                    // Walk forward to mid point.
                    // - if next is +2 midpoint is +1
                    // - if next is +1 it is midpoint
                    // - if next is +0 there has to be edgeMidPoint

                    walkFaceToMid
                    (
                        edgeMidPoint,
                        anchorLevel,
                        faceI,
                        fp,
                        faceVerts
                    );

                    faceVerts.append(faceMidPoint[faceI]);

                    walkFaceFromMid
                    (
                        edgeMidPoint,
                        anchorLevel,
                        faceI,
                        fp,
                        faceVerts
                    );

                    // Convert dynamiclist to face.
                    newFace.transfer(faceVerts);

                    //Pout<< "Split face:" << faceI << " verts:" << f
                    //    << " into quad:" << newFace << endl;

                    // Get new owner/neighbour
                    label own, nei;
                    getFaceNeighbours
                    (
                        cellAnchorPoints,
                        cellAddedCells,
                        faceI,
                        pointI,          // Anchor point

                        own,
                        nei
                    );


                    if (debug)
                    {
                        if (mesh_.isInternalFace(faceI))
                        {
                            label oldOwn = mesh_.faceOwner()[faceI];
                            label oldNei = mesh_.faceNeighbour()[faceI];

                            checkInternalOrientation
                            (
                                meshMod,
                                oldOwn,
                                faceI,
                                mesh_.cellCentres()[oldOwn],
                                mesh_.cellCentres()[oldNei],
                                newFace
                            );
                        }
                        else
                        {
                            label oldOwn = mesh_.faceOwner()[faceI];

                            checkBoundaryOrientation
                            (
                                meshMod,
                                oldOwn,
                                faceI,
                                mesh_.cellCentres()[oldOwn],
                                mesh_.faceCentres()[faceI],
                                newFace
                            );
                        }
                    }


                    if (!modifiedFace)
                    {
                        modifiedFace = true;

                        modFace(meshMod, faceI, newFace, own, nei);
                    }
                    else
                    {
                        addFace(meshMod, faceI, newFace, own, nei);
                    }
                }
            }

            // Mark face as having been handled
            affectedFace.unset(faceI);
            faceNeedsSplit[faceI] = false;
        }
        else if
        (
            faceNeedsSplit[faceI] &&
            affectedFace.get(faceI) &&
            (
                (cellMidPoint[mesh_.faceOwner()[faceI]] >= 0)
             || (
                    mesh_.isInternalFace(faceI)
                    ? cellMidPoint[mesh_.faceNeighbour()[faceI]] >= 0
                    : false
                )
            )
        )
        {
            // this is a tri-face on a hex-cell that needs split.

            const face &f = meshMod.faces()[faceI];
            label fLevel = faceLevelTet(faceI);
            label ownOld = mesh_.faceOwner()[faceI];
            label neiOld = -1;
            if (mesh_.isInternalFace(faceI))
            {
                neiOld = mesh_.faceNeighbour()[faceI];
            }

            labelList fAnchors(3, -1);
            labelList fMids(3, -1);

            getFaceAnchorsAndMids
            (
                faceI,
                fAnchors,
                fMids,
                fLevel,
                newPointLevel,
                meshMod
            );

            faceList fList(1);
            fList[0] = f;

            // add anchorFaces
            for (label j = 0; j < 3; j++)
            {
                DynamicList<label> faceVerts;

                labelList nodeLooper(3, -1);
                nodeLooper[0] = fAnchors[j];
                nodeLooper[1] = fMids[j];
                nodeLooper[2] = fMids[fMids.rcIndex(j)];

                forAll(nodeLooper, l)
                {
                    walkBetweenPoints
                    (
                        nodeLooper[l],
                        nodeLooper[nodeLooper.fcIndex(l)],
                        fList,
                        newPointLevel,
                        fLevel+1,
                        faceVerts,
                        false,
                        true
                    );
                }

                face newFace(faceVerts);

                label ownNew = -1;
                label neiNew = -1;

                label nOwnAnchors = countAnchors
                (
                    mesh_.cellPoints()[ownOld],
                    cellLevel_[ownOld]
                );

                if (nOwnAnchors == 8)
                {
                    ownNew = getAnchorCell
                    (
                        cellAnchorPoints,
                        cellAddedCells,
                        ownOld,
                        -1,
                        fAnchors[j]
                    );

                    if (neiOld != -1)
                    {
                        labelList adjAnchors = getAnchors(neiOld);

                        labelListList adjAnchorToes = getAnchorToes
                        (
                            neiOld,
                            adjAnchors,
                            edgeMidPoint
                        );

                        neiNew = getSubCellID
                        (
                            newFace,
                            faceI,
                            cellLevel_[neiOld]+1,
                            newPointLevel,
                            adjAnchorToes,
                            neiOld,
                            adjAnchors,
                            cellAddedCells,
                            meshMod
                        );
                    }
                }
                else if (nOwnAnchors == 4)
                {
                    labelList adjAnchors = getAnchors(ownOld);

                    labelListList adjAnchorToes = getAnchorToes
                    (
                        ownOld,
                        adjAnchors,
                        edgeMidPoint
                    );

                    ownNew = getSubCellID
                    (
                        newFace,
                        faceI,
                        cellLevel_[ownOld]+1,
                        newPointLevel,
                        adjAnchorToes,
                        ownOld,
                        adjAnchors,
                        cellAddedCells,
                        meshMod
                    );

                    if (neiOld != -1)
                    {
                        neiNew = getAnchorCell
                        (
                            cellAnchorPoints,
                            cellAddedCells,
                            neiOld,
                            -1,
                            fAnchors[j]
                        );
                    }
                }

                addFace(meshMod, faceI, newFace, ownNew, neiNew);
            }

            // mod midFace
            DynamicList<label> faceVerts;

            forAll(fMids, j)
            {
                walkBetweenPoints
                (
                    fMids[j],
                    fMids[fMids.fcIndex(j)],
                    fList,
                    newPointLevel,
                    fLevel+1,
                    faceVerts,
                    false,
                    true
                );
            }

            face newFace(faceVerts);

            label ownNew = -1;
            label neiNew = -1;

            label nOwnAnchors = countAnchors
            (
                mesh_.cellPoints()[ownOld],
                cellLevel_[ownOld]
            );

            if (nOwnAnchors == 8)
            {
                label relevantAnchor = -1;
                {
                    label relevantMid;
                    getRelevantAnchorAndMid
                    (
                        relevantAnchor,
                        relevantMid,
                        faceI,
                        mesh_.cells()[ownOld],
                        cellLevel_[ownOld],
                        newPointLevel,
                        edgeMidPoint,
                        meshMod,
                        0.174533
                    );
                }

                ownNew = getAnchorCell
                (
                    cellAnchorPoints,
                    cellAddedCells,
                    ownOld,
                    -1,
                    relevantAnchor
                );

                if (neiOld != -1)
                {
                    labelList adjAnchors = getAnchors(neiOld);
                    labelListList adjAnchorToes = getAnchorToes
                    (
                        neiOld,
                        adjAnchors,
                        edgeMidPoint
                    );

                    neiNew = getSubCellID
                    (
                        newFace,
                        faceI,
                        cellLevel_[neiOld]+1,
                        newPointLevel,
                        adjAnchorToes,
                        neiOld,
                        adjAnchors,
                        cellAddedCells,
                        meshMod
                    );
                }
            }
            else if (nOwnAnchors == 4)
            {
                labelList adjAnchors = getAnchors(ownOld);
                labelListList adjAnchorToes = getAnchorToes
                (
                    ownOld,
                    adjAnchors,
                    edgeMidPoint
                );

                ownNew = getSubCellID
                (
                    newFace,
                    faceI,
                    cellLevel_[ownOld]+1,
                    newPointLevel,
                    adjAnchorToes,
                    ownOld,
                    adjAnchors,
                    cellAddedCells,
                    meshMod
                );

                if (neiOld != -1)
                {
                    label relevantAnchor = -1;
                    {
                        label relevantMid;
                        getRelevantAnchorAndMid
                        (
                            relevantAnchor,
                            relevantMid,
                            faceI,
                            mesh_.cells()[neiOld],
                            cellLevel_[neiOld],
                            newPointLevel,
                            edgeMidPoint,
                            meshMod,
                            0.174533
                        );
                    }
                    neiNew = getAnchorCell
                    (
                        cellAnchorPoints,
                        cellAddedCells,
                        neiOld,
                        -1,
                        relevantAnchor
                    );
                }
            }

            modFace(meshMod, faceI, newFace, ownNew, neiNew);

            faceNeedsSplit[faceI] = false;
            affectedFace.unset(faceI);
        }
    }


    // 2. faces that do not get split but whose owner/neighbour change
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexTetRef8::setHexRefinement :"
            << " Changing owner/neighbour for otherwise unaffected faces"
            << endl;
    }

    forAll(affectedFace, faceI)
    {
        if (affectedFace.get(faceI) && !faceNeedsSplit[faceI])
        {
            const face& f = meshMod.faces()[faceI];

            // The point with the lowest level should be an anchor
            // point of the neighbouring cells.
            //label anchorFp = findMinLevel(f);
            label minLevel = labelMax;
            label anchorFp = -1;

            forAll(f, fp)
            {
                label level = newPointLevel[f[fp]];

                if (level < minLevel)
                {
                    minLevel = level;
                    anchorFp = fp;
                }
            }

            label own, nei;
            getFaceNeighbours
            (
                cellAnchorPoints,
                cellAddedCells,
                faceI,
                f[anchorFp],          // Anchor point

                own,
                nei
            );

            modFace(meshMod, faceI, f, own, nei);

            // Mark face as having been handled
            affectedFace.unset(faceI);
        }
    }


    // 3. new internal faces inside split cells.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    // This is the hard one. We have to find the splitting points between
    // the anchor points. But the edges between the anchor points might have
    // been split (into two,three or four edges).

    if (debug)
    {
        Pout<< "hexTetRef8::setHexRefinement :"
            << " Create new internal faces for split cells"
            << endl;
    }

    forAll(cellMidPoint, cellI)
    {
        if (cellMidPoint[cellI] >= 0)
        {
            createInternalFaces
            (
                cellAnchorPoints,
                cellAddedCells,
                cellMidPoint,
                faceMidPoint,
                faceAnchorLevel,
                edgeMidPoint,
                newPointLevel,
                cellI,
                meshMod
            );
        }
    }

    // Extend pointLevels and cellLevels for the new cells. Could also be done
    // in updateMesh but saves passing cellAddedCells out of this routine.

    // Check
    if (debug)
    {
        label minPointI = labelMax;
        label maxPointI = labelMin;

        forAll(cellMidPoint, cellI)
        {
            if (cellMidPoint[cellI] >= 0)
            {
                minPointI = min(minPointI, cellMidPoint[cellI]);
                maxPointI = max(maxPointI, cellMidPoint[cellI]);
            }
        }
        forAll(faceMidPoint, faceI)
        {
            if (faceMidPoint[faceI] >= 0)
            {
                minPointI = min(minPointI, faceMidPoint[faceI]);
                maxPointI = max(maxPointI, faceMidPoint[faceI]);
            }
        }
        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                minPointI = min(minPointI, edgeMidPoint[edgeI]);
                maxPointI = max(maxPointI, edgeMidPoint[edgeI]);
            }
        }

        if (minPointI != labelMax && minPointI != mesh_.nPoints())
        {
            FatalErrorInFunction
                << "Added point labels not consecutive to existing mesh points."
                << nl
                << "mesh_.nPoints():" << mesh_.nPoints()
                << " minPointI:" << minPointI
                << " maxPointI:" << maxPointI
                << abort(FatalError);
        }
    }
}


//- Insert refinement for tet cells.
//  Mapping:
//  -split cells: 7 new ones get added from original
//  -split faces: original gets modified; 3 new ones get added
//               from original
//  -added internal faces: created out-of-nothing (so will not
//   get mapped!). Note: could make this inflate from point but
//   that will allocate interpolation.
//  -points added to split edge: added from edge start()
void Foam::hexTetRef8::setTetRefinement
(
    const labelList& tetCellLabels,
    const boolList& cellNeedsSplit,
    const boolList& faceNeedsSplit,
    const labelList& newEdgeMids,
    const labelListList& cellAddedCells,
    const DynamicList<label>& newPointLevel,
    polyTopoChange& meshMod
)
{
    // Refinement decision per refined cell.
    boolList isCellTetSplit(mesh_.nCells(), false);

    forAll(tetCellLabels, i)
    {
        label cellI = tetCellLabels[i];

        isCellTetSplit[cellI] = true;
    }

    // get split edges mid points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexTetRef8::setTetRefinement :"
            << " Getting edge midpoints."
            << endl;
    }

    // Face anchor level. There are guaranteed 3 OR 4 points with level
    // <= anchorLevel. These are the corner points.
    labelList faceAnchorLevel(mesh_.nFaces());
    // Number of face anchors (points of level <= faceAnchorLevel)
    labelList nFaceAnchors(mesh_.nFaces());

    for (label faceI = 0; faceI < mesh_.nFaces(); faceI++)
    {
        faceAnchorLevel[faceI] = faceLevel(faceI);
        if (faceAnchorLevel[faceI] == -1)
        {
            faceAnchorLevel[faceI] = faceLevelTet(faceI);
        }

        nFaceAnchors[faceI] = countAnchors
        (
            mesh_.faces()[faceI],
            faceAnchorLevel[faceI]
        );
    }

    // -1  : no need to split edge
    // >=0 : label of introduced mid point (calc in setRefinement)
    labelList edgeMidPoint(mesh_.nEdges(), -1);

    forAll(faceNeedsSplit, faceI)
    {
        if (faceNeedsSplit[faceI])
        {
            if (nFaceAnchors[faceI] == 3)
            {
                const labelList& fEdges = mesh_.faceEdges()[faceI];

                forAll(fEdges, j)
                {
                    label edgeJ = fEdges[j];

                    if (newEdgeMids[edgeJ] >= 0)
                    {
                        edgeMidPoint[edgeJ] = newEdgeMids[edgeJ];
                    }
                }
            }
        }
    }

    if (debug)
    {
        OFstream str(mesh_.time().path()/"edgeMidPoint.obj");

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                const edge& e = mesh_.edges()[edgeI];

                meshTools::writeOBJ(str, e.centre(mesh_.points()));
            }
        }

        Pout<< "hexTetRef8::setTetRefinement :"
            << " Dumping edge centres to split to file " << str.name() << endl;
    }


    // Find if face is Split
    // ~~~~~~~~~~~~~~~~~~~~~

    // -2  : no need to split face
    // -1  : split face
    // >=0 : already split
    labelList isFaceSplit(mesh_.nFaces(), -2);

    for (label faceI = 0; faceI < mesh_.nFaces(); faceI++)
    {
        if (faceNeedsSplit[faceI])
        {
            // we only add a mid-point for square faces
            if (nFaceAnchors[faceI] == 3)
            {
                isFaceSplit[faceI] = -1;    // mark to add mid-point
            }
        }
    }

    // Synchronize isFaceSplit across coupled patches. (logical or)
    syncTools::syncFaceList
    (
        mesh_,
        isFaceSplit,
        maxEqOp<label>()
    );

    // split boundary faces whose owners are not refined
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexTetRef8::setTetRefinement :"
            << " Splitting boundary faces whose owners are not refined."
            << endl;
    }

    label nBoundaryFaces = mesh_.nFaces()-mesh_.nInternalFaces();

    for (label i = 0; i < nBoundaryFaces; i++)
    {
        label faceI = i+mesh_.nInternalFaces();
        label own = mesh_.faceOwner()[faceI];

        if ((isFaceSplit[faceI] == -1) && !cellNeedsSplit[own])
        {
            const face &f = meshMod.faces()[faceI];
            label fLevel = faceLevelTet(faceI);

            labelList fAnchors(3, -1);
            labelList fMids(3, -1);

            getFaceAnchorsAndMids
            (
                faceI,
                fAnchors,
                fMids,
                fLevel,
                newPointLevel,
                meshMod
            );

            faceList fList(1);
            fList[0] = f;

            // add anchorFaces
            for (label j = 0; j < 3; j++)
            {
                DynamicList<label> faceVerts;

                labelList nodeLooper(3, -1);
                nodeLooper[0] = fAnchors[j];
                nodeLooper[1] = fMids[j];
                nodeLooper[2] = fMids[fMids.rcIndex(j)];

                forAll(nodeLooper, l)
                {
                    walkBetweenPoints
                    (
                        nodeLooper[l],
                        nodeLooper[nodeLooper.fcIndex(l)],
                        fList,
                        newPointLevel,
                        fLevel+1,
                        faceVerts,
                        false,
                        true
                    );
                }

                face newFace(faceVerts);

                addFace(meshMod, faceI, newFace, own, -1);
            }

            // mod midFace
            DynamicList<label> faceVerts;

            forAll(fMids, j)
            {
                walkBetweenPoints
                (
                    fMids[j],
                    fMids[fMids.fcIndex(j)],
                    fList,
                    newPointLevel,
                    fLevel+1,
                    faceVerts,
                    false,
                    true
                );
            }

            face newFace(faceVerts);

            modFace(meshMod, faceI, newFace, own, -1);

            isFaceSplit[faceI] = 1; // mark face as split
        }
    }

    // Loop over cells:
    // ~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexTetRef8::setTetRefinement :"
            << " Refining cellFaces of cells marked for refinement"
            << endl;
    }

    forAll (tetCellLabels, i)
    {
        label       cellI           = tetCellLabels[i];
        labelList   cellIFaceIDs    = mesh_.cells()[cellI];
        faceList    cellIFaces      (cellIFaceIDs.size());
        forAll(cellIFaces, j)
        {
//            cellIFaces[j] = meshMod.faces()[cellIFaceIDs[j]];

            const face& f = mesh_.faces()[cellIFaceIDs[j]];
            DynamicList<label> fEdgesStorage;

            const labelList& fEdges = mesh_.faceEdges
            (
                cellIFaceIDs[j],
                fEdgesStorage
            );

            DynamicList<label> newFaceVerts(f.size());

            forAll(f, fp)
            {
                newFaceVerts.append(f[fp]);

                label edgeJ = fEdges[fp];

                if (newEdgeMids[edgeJ] >= 0)
                {
                    newFaceVerts.append(newEdgeMids[edgeJ]);
                }
            }

            cellIFaces[j].transfer(newFaceVerts);
        }

        // Array of anchors at each cell. The anchors are stored
        // such that the first three are in an outward pointing
        // order
        labelList anchors = getAnchors(cellI);

        if (anchors.size() != 4)
        {
            FatalErrorInFunction
                << "Alleged Tet cell " << cellI << " has "
                << anchors.size() << " anchors. "
                << "This must be 4."
                << abort(FatalError);
        }

        //  The sub-cells for cellI.
        //  CellAddedCells contains 8 sub-cells for each cell:
        //      -> the first four are cells adjacent to the
        //         anchors (in order)
        //      -> the next four are the cells adjacent to the
        //         center faces of the cellFaces
        //              -- b/w anchors 0-1-2
        //              -- b/w anchors 0-1-3
        //              -- b/w anchors 0-2-3
        //              -- b/w anchors 1-2-3
        //      -> the last centerFace to be modded recycles
        //         the original cell

        const labelList &cellIAddedCells = cellAddedCells[cellI];

        // Toes for each anchor.
        /*
         *             A1      A2      A3      A4
         *         --                             --
         *    A1 - |   -1    mid12   mid13   mid14 |
         *         |                               |
         *    A2 - | mid12     -1    mid23   mid24 |
         *         |                               |
         *    A3 - | mid13   mid23     -1    mid34 |
         *         |                               |
         *    A4 - | mid14   mid24   mid34     -1  |
         *         --                             --
         */
        labelListList anchorToes = getAnchorToes(cellI, anchors, edgeMidPoint);

        // external subFaces that do not use any anchor point.
        // These faces might or might not be subdivided. The
        // order of points in these faces always points OUT of
        // cellI.

        faceList centerFaces = getCenterFaces
        (
            anchorToes,
            cellIFaces,
            newPointLevel,
            cellLevel_[cellI]+1
        );

        // 1. Add / mod external faces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        forAll(cellIFaceIDs, j)
        {
            label faceJ = cellIFaceIDs[j];

            if (isFaceSplit[faceJ] == -2)
            {
                // doesn't need split: correct them for their own/nei

                face newFace = meshMod.faces()[faceJ];
                label ownNew = -1;
                label neiNew = -1;

                getNeighbours
                (
                    meshMod.faces()[faceJ],
                    faceJ,
                    cellI,
                    anchors,
                    anchorToes,
                    cellAddedCells,
                    ownNew,
                    neiNew,
                    edgeMidPoint,
                    newPointLevel,
                    meshMod
                );

                modFace(meshMod, faceJ, newFace, ownNew, neiNew);

                isFaceSplit[faceJ] = 1;
            }
            else
            {
                if (isFaceSplit[faceJ] == -1)
                {
                    // needs split, but not yet split: split into
                    // sub faces and give them right own/nei

                    labelList faceAnchors(3, -1);
                    labelList mids(3, -1);
                    label referenceLevel;

                    if
                    (
                        countAnchors
                        (
                            mesh_.faces()[faceJ],
                            cellLevel_[cellI]
                        ) != 3
                    )
                    {
                        // face is already refined w.r.t this cell, but needs
                        // refinement for adjacent cell.
                        referenceLevel = cellLevel_[cellI]+2;
                    }
                    else
                    {
                        referenceLevel = cellLevel_[cellI]+1;
                    }

                    // get list of face anchors and their
                    // mid-points in the original order
                    getFaceAnchorsAndMids
                    (
                        faceJ,
                        faceAnchors,
                        mids,
                        referenceLevel-1,
                        newPointLevel,
                        meshMod
                    );

                    label ownNew = -1;
                    label neiNew = -1;

                    // add anchor faces

                    forAll(faceAnchors, k)
                    {
                        // construct the face in the same own/nei orientation
                        // as the parent face. The order will be corrected by
                        // addFace

                        DynamicList<label> faceVerts;

                        labelList nodeLooper(3, -1);
                        nodeLooper[0] = faceAnchors[k];
                        nodeLooper[1] = mids[k];
                        nodeLooper[2] = mids[mids.rcIndex(k)];

                        forAll(nodeLooper, l)
                        {
                            walkBetweenPoints
                            (
                                nodeLooper[l],
                                nodeLooper[nodeLooper.fcIndex(l)],
                                cellIFaces,
                                newPointLevel,
                                referenceLevel,
                                faceVerts,
                                false,
                                true
                            );
                        }

                        face newFace(faceVerts);

                        getNeighbours
                        (
                            newFace,
                            faceJ,
                            cellI,
                            anchors,
                            anchorToes,
                            cellAddedCells,
                            ownNew,
                            neiNew,
                            edgeMidPoint,
                            newPointLevel,
                            meshMod
                        );

                        addFace(meshMod, faceJ, newFace, ownNew, neiNew);
                    }

                    // mod center face

                    DynamicList<label> faceVerts;

                    // get points from anchor to next mid

                    forAll(mids, k)
                    {
                        walkBetweenPoints
                        (
                            mids[k],
                            mids[mids.fcIndex(k)],
                            cellIFaces,
                            newPointLevel,
                            referenceLevel,
                            faceVerts,
                            false,
                            true
                        );
                    }

                    face newFace(faceVerts);

                    getNeighbours
                    (
                        newFace,
                        faceJ,
                        cellI,
                        anchors,
                        anchorToes,
                        cellAddedCells,
                        ownNew,
                        neiNew,
                        edgeMidPoint,
                        newPointLevel,
                        meshMod
                    );

                    modFace(meshMod, faceJ, newFace, ownNew, neiNew);

                    isFaceSplit[faceJ] = 1;
                }
                else if (isFaceSplit[faceJ] >= 0)
                {
                    // needs split, and has already been split in some
                    // other cell.
                }
                else
                {
                    FatalErrorInFunction
                    << "Invalid value " << isFaceSplit[faceJ]
                    << " for face " << faceJ << " in faceSubFacesIndex."
                    << abort(FatalError);
                }
            }
        }


        // Add internal faces
        // ~~~~~~~~~~~~~~~~~~
        //
        // Internal faces are those that do not have any own/nei outside of
        // cellAddedCells for that cell. There are two types of internal faces:
        //      1. internalAnchorFaces: these are those internal faces that lie
        //         next to an anchor. Each internalAnchorFace corresponds to
        //         one unique anchor. Thus, there are 4 of these.
        //      2. internalDiagFaces: these internal faces have the diagonal as
        //         one edge. For each anchor toe (except those on the diagonal)
        //         we have one internalDiagFace. Thus, there are 4 of these.

        // find shortest diagonal
        Pair<label> diag = chooseDiagonal(anchorToes, meshMod);

        // flag for checking if internal anchorFaces have been added
        boolList isInternalAnchorFaceAdded(4, false);

        // check all non-redundant anchor pairs. If the pair's mid point
        // is not one of the diagPoints, work on them

        forAll(anchors, i)
        {
            label as = i; // source anchor index
            label ad = anchors.fcIndex(i);  // destination anchor index

            for (int j = i+1; j < 4; j++)
            {
                // mid point for pair
                label mid = anchorToes[as][ad];

                if
                (
                    diag.first()  == mid ||
                    diag.second() == mid
                )
                {
                    // pair mid is one of the diag points.
                    // ignore.
                }
                else
                {
                    Pair<label> anchorPair(as, ad);

                    // we now need to get those centerFaces that have
                    // "mid" as one of their points. Also get their
                    // face points and the corresponding subCells

                    faceList mCF(2);
                    labelList mCFIndices(2, -1);
                    labelList mCFCells(2, -1);

                    getCenterFacesWithPoint
                    (
                        centerFaces,
                        mid,
                        cellIAddedCells,
                        mCF,
                        mCFIndices,
                        mCFCells
                    );

                    // Add the face that connects mid and the diagonal
                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                    label ownNew = -1;
                    label neiNew = -1;
                    bool slicerOutOfMCFCell0
                    = (mCFCells[0] > mCFCells[1])
                        ? false
                        : true;

                    face newFace = getNewInternalFace
                    (
                        mCF,
                        mCF[0],
                        mCFCells[0],
                        mCFCells[1],
                        mid,
                        diag,
                        slicerOutOfMCFCell0,
                        ownNew,
                        neiNew,
                        0,
                        newPointLevel,
                        cellLevel_[cellI]+1,
                        meshMod
                    );

                    label parentFace = -1;
                    addFace(meshMod, parentFace, newFace, ownNew, neiNew);


                    // Add Internal anchor faces if needed
                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                    forAll(anchorPair, i)
                    {
                        if (!isInternalAnchorFaceAdded[anchorPair[i]])
                        {
                            label anchorMCFIndex = findAnchorMCF
                            (
                                mCF,
                                anchorToes[anchorPair[i]],
                                diag
                            );

                            label anchorCell = cellIAddedCells[anchorPair[i]];

                            bool anchorInternalFaceOutOfMCFCellI
                                = mCFCells[anchorMCFIndex] > anchorCell
                                    ? false
                                    : true;

                            label ownNew = -1;
                            label neiNew = -1;

                            face newFace = getNewInternalFace
                            (
                                centerFaces,
                                mCF[anchorMCFIndex],
                                mCFCells[anchorMCFIndex],
                                anchorCell,
                                mid,
                                diag,
                                anchorInternalFaceOutOfMCFCellI,
                                ownNew,
                                neiNew,
                                1,
                                newPointLevel,
                                cellLevel_[cellI]+1,
                                meshMod
                            );

                            label parentFace = -1;
                            addFace(meshMod, parentFace, newFace, ownNew, neiNew);

                            isInternalAnchorFaceAdded[anchorPair[i]] = true;
                        }
                    }
                }

                // move to the next anchor
                ad = anchors.fcIndex(ad);
            }
        }
    }
}


void Foam::hexTetRef8::storeData
(
    const labelList& pointsToStore,
    const labelList& facesToStore,
    const labelList& cellsToStore
)
{
    savedPointLevel_.resize(2*pointsToStore.size());
    forAll(pointsToStore, i)
    {
        label pointI = pointsToStore[i];
        savedPointLevel_.insert(pointI, pointLevel_[pointI]);
    }

    savedCellLevel_.resize(2*cellsToStore.size());
    forAll(cellsToStore, i)
    {
        label cellI = cellsToStore[i];
        savedCellLevel_.insert(cellI, cellLevel_[cellI]);
    }
}


// Gets called after the mesh change. setRefinement will already have made
// sure the pointLevel_ and cellLevel_ are the size of the new mesh so we
// only need to account for reordering.
void Foam::hexTetRef8::updateMesh(const mapPolyMesh& map)
{
    Map<label> dummyMap(0);

    updateMesh(map, dummyMap, dummyMap, dummyMap);
}


// Gets called after the mesh change. setRefinement will already have made
// sure the pointLevel_ and cellLevel_ are the size of the new mesh so we
// only need to account for reordering.
void Foam::hexTetRef8::updateMesh
(
    const mapPolyMesh& map,
    const Map<label>& pointsToRestore,
    const Map<label>& facesToRestore,
    const Map<label>& cellsToRestore
)
{
    // Update celllevel
    if (debug)
    {
        Pout<< "hexTetRef8::updateMesh :"
            << " Updating various lists"
            << endl;
    }

    {
        const labelList& reverseCellMap = map.reverseCellMap();

        if (debug)
        {
            Pout<< "hexTetRef8::updateMesh :"
                << " reverseCellMap:" << map.reverseCellMap().size()
                << " cellMap:" << map.cellMap().size()
                << " nCells:" << mesh_.nCells()
                << " nOldCells:" << map.nOldCells()
                << " cellLevel_:" << cellLevel_.size()
                << " reversePointMap:" << map.reversePointMap().size()
                << " pointMap:" << map.pointMap().size()
                << " nPoints:" << mesh_.nPoints()
                << " nOldPoints:" << map.nOldPoints()
                << " pointLevel_:" << pointLevel_.size()
                << endl;
        }

        if (reverseCellMap.size() == cellLevel_.size())
        {
            // Assume it is after hexTetRef8 that this routine is called.
            // Just account for reordering. We cannot use cellMap since
            // then cells created from cells would get cellLevel_ of
            // cell they were created from.
            reorder(reverseCellMap, mesh_.nCells(), -1, cellLevel_);
        }
        else
        {
            // Map data
            const labelList& cellMap = map.cellMap();

            labelList newCellLevel(cellMap.size());
            forAll(cellMap, newCellI)
            {
                label oldCellI = cellMap[newCellI];

                if (oldCellI == -1)
                {
                    newCellLevel[newCellI] = -1;
                }
                else
                {
                    newCellLevel[newCellI] = cellLevel_[oldCellI];
                }
            }
            cellLevel_.transfer(newCellLevel);
        }

        // See if any cells to restore. This will be for some new cells
        // the corresponding old cell.
        forAllConstIter(Map<label>, cellsToRestore, iter)
        {
            label newCellI = iter.key();
            label storedCellI = iter();

            Map<label>::iterator fnd = savedCellLevel_.find(storedCellI);

            if (fnd == savedCellLevel_.end())
            {
                FatalErrorInFunction
                    << "Problem : trying to restore old value for new cell "
                    << newCellI << " but cannot find old cell " << storedCellI
                    << " in map of stored values " << savedCellLevel_
                    << abort(FatalError);
            }
            cellLevel_[newCellI] = fnd();
        }

        //if (findIndex(cellLevel_, -1) != -1)
        //{
        //    WarningInFunction
        //        << "Problem : "
        //        << "cellLevel_ contains illegal value -1 after mapping
        //        << " at cell " << findIndex(cellLevel_, -1) << endl
        //        << "This means that another program has inflated cells"
        //        << " (created cells out-of-nothing) and hence we don't know"
        //        << " their cell level. Continuing with illegal value."
        //        << abort(FatalError);
        //}
    }


    // Update pointlevel
    {
        const labelList& reversePointMap = map.reversePointMap();

        if (reversePointMap.size() == pointLevel_.size())
        {
            // Assume it is after hexTetRef8 that this routine is called.
            reorder(reversePointMap, mesh_.nPoints(), -1,  pointLevel_);
        }
        else
        {
            // Map data
            const labelList& pointMap = map.pointMap();

            labelList newPointLevel(pointMap.size());

            forAll(pointMap, newPointI)
            {
                label oldPointI = pointMap[newPointI];

                if (oldPointI == -1)
                {
                    //FatalErrorInFunction
                    //    << "Problem : point " << newPointI
                    //    << " at " << mesh_.points()[newPointI]
                    //    << " does not originate from another point"
                    //    << " (i.e. is inflated)." << nl
                    //    << "Hence we cannot determine the new pointLevel"
                    //    << " for it." << abort(FatalError);
                    newPointLevel[newPointI] = -1;
                }
                else
                {
                    newPointLevel[newPointI] = pointLevel_[oldPointI];
                }
            }
            pointLevel_.transfer(newPointLevel);
        }

        // See if any points to restore. This will be for some new points
        // the corresponding old point (the one from the call to storeData)
        forAllConstIter(Map<label>, pointsToRestore, iter)
        {
            label newPointI = iter.key();
            label storedPointI = iter();

            Map<label>::iterator fnd = savedPointLevel_.find(storedPointI);

            if (fnd == savedPointLevel_.end())
            {
                FatalErrorInFunction
                    << "Problem : trying to restore old value for new point "
                    << newPointI << " but cannot find old point "
                    << storedPointI
                    << " in map of stored values " << savedPointLevel_
                    << abort(FatalError);
            }
            pointLevel_[newPointI] = fnd();
        }

        //if (findIndex(pointLevel_, -1) != -1)
        //{
        //    WarningInFunction
        //        << "Problem : "
        //        << "pointLevel_ contains illegal value -1 after mapping"
        //        << " at point" << findIndex(pointLevel_, -1) << endl
        //        << "This means that another program has inflated points"
        //        << " (created points out-of-nothing) and hence we don't know"
        //        << " their point level. Continuing with illegal value."
        //        //<< abort(FatalError);
        //}
    }

    // Update refinement tree
    if (history_.active())
    {
        history_.updateMesh(map);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // Update face removal engine
    faceRemover_.updateMesh(map);

    // Clear cell shapes
    cellShapesPtr_.clear();
}


// Gets called after mesh subsetting. Maps are from new back to old.
void Foam::hexTetRef8::subset
(
    const labelList& pointMap,
    const labelList& faceMap,
    const labelList& cellMap
)
{
    // Update celllevel
    if (debug)
    {
        Pout<< "hexTetRef8::subset :"
            << " Updating various lists"
            << endl;
    }

    if (history_.active())
    {
        WarningInFunction
            << "Subsetting will not work in combination with unrefinement."
            << nl
            << "Proceed at your own risk." << endl;
    }


    // Update celllevel
    {
        labelList newCellLevel(cellMap.size());

        forAll(cellMap, newCellI)
        {
            newCellLevel[newCellI] = cellLevel_[cellMap[newCellI]];
        }

        cellLevel_.transfer(newCellLevel);

        if (findIndex(cellLevel_, -1) != -1)
        {
            FatalErrorInFunction
                << "Problem : "
                << "cellLevel_ contains illegal value -1 after mapping:"
                << cellLevel_
                << abort(FatalError);
        }
    }

    // Update pointlevel
    {
        labelList newPointLevel(pointMap.size());

        forAll(pointMap, newPointI)
        {
            newPointLevel[newPointI] = pointLevel_[pointMap[newPointI]];
        }

        pointLevel_.transfer(newPointLevel);

        if (findIndex(pointLevel_, -1) != -1)
        {
            FatalErrorInFunction
                << "Problem : "
                << "pointLevel_ contains illegal value -1 after mapping:"
                << pointLevel_
                << abort(FatalError);
        }
    }

    // Update refinement tree
    if (history_.active())
    {
        history_.subset(pointMap, faceMap, cellMap);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // Nothing needs doing to faceRemover.
    //faceRemover_.subset(pointMap, faceMap, cellMap);

    // Clear cell shapes
    cellShapesPtr_.clear();
}


// Gets called after the mesh distribution
void Foam::hexTetRef8::distribute(const mapDistributePolyMesh& map)
{
    if (debug)
    {
        Pout<< "hexTetRef8::distribute :"
            << " Updating various lists"
            << endl;
    }

    // Update celllevel
    map.distributeCellData(cellLevel_);
    // Update pointlevel
    map.distributePointData(pointLevel_);

    // Update refinement tree
    if (history_.active())
    {
        history_.distribute(map);
    }

    // Update face removal engine
    faceRemover_.distribute(map);

    // Clear cell shapes
    cellShapesPtr_.clear();
}


void Foam::hexTetRef8::checkMesh() const
{
    const scalar smallDim = 1e-6 * mesh_.bounds().mag();

    if (debug)
    {
        Pout<< "hexTetRef8::checkMesh : Using matching tolerance smallDim:"
            << smallDim << endl;
    }

    // Check owner on coupled faces.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // There should be only one coupled face between two cells. Why? Since
    // otherwise mesh redistribution might cause multiple faces between two
    // cells
    {
        labelList nei(mesh_.nFaces()-mesh_.nInternalFaces());
        forAll(nei, i)
        {
            nei[i] = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
        }

        // Replace data on coupled patches with their neighbour ones.
        syncTools::swapBoundaryFaceList(mesh_, nei);

        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (pp.coupled())
            {
                // Check how many faces between owner and neighbour. Should
                // be only one.
                HashTable<label, labelPair, labelPair::Hash<> >
                    cellToFace(2*pp.size());

                label faceI = pp.start();

                forAll(pp, i)
                {
                    label own = mesh_.faceOwner()[faceI];
                    label bFaceI = faceI-mesh_.nInternalFaces();

                    if (!cellToFace.insert(labelPair(own, nei[bFaceI]), faceI))
                    {
                        dumpCell(own);
                        FatalErrorInFunction
                            << "Faces do not seem to be correct across coupled"
                            << " boundaries" << endl
                            << "Coupled face " << faceI
                            << " between owner " << own
                            << " on patch " << pp.name()
                            << " and coupled neighbour " << nei[bFaceI]
                            << " has two faces connected to it:"
                            << faceI << " and "
                            << cellToFace[labelPair(own, nei[bFaceI])]
                            << abort(FatalError);
                    }

                    faceI++;
                }
            }
        }
    }

    // Check face areas.
    // ~~~~~~~~~~~~~~~~~

    {
        scalarField neiFaceAreas(mesh_.nFaces()-mesh_.nInternalFaces());
        forAll(neiFaceAreas, i)
        {
            neiFaceAreas[i] = mag(mesh_.faceAreas()[i+mesh_.nInternalFaces()]);
        }

        // Replace data on coupled patches with their neighbour ones.
        syncTools::swapBoundaryFaceList(mesh_, neiFaceAreas);

        forAll(neiFaceAreas, i)
        {
            label faceI = i+mesh_.nInternalFaces();

            const scalar magArea = mag(mesh_.faceAreas()[faceI]);

            if (mag(magArea - neiFaceAreas[i]) > smallDim)
            {
                const face& f = mesh_.faces()[faceI];
                label patchI = mesh_.boundaryMesh().whichPatch(faceI);

                dumpCell(mesh_.faceOwner()[faceI]);

                FatalErrorInFunction
                    << "Faces do not seem to be correct across coupled"
                    << " boundaries" << endl
                    << "Coupled face " << faceI
                    << " on patch " << patchI
                    << " " << mesh_.boundaryMesh()[patchI].name()
                    << " coords:" << UIndirectList<point>(mesh_.points(), f)()
                    << " has face area:" << magArea
                    << " (coupled) neighbour face area differs:"
                    << neiFaceAreas[i]
                    << " to within tolerance " << smallDim
                    << abort(FatalError);
            }
        }
    }


    // Check number of points on faces.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    {
        labelList nVerts(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(nVerts, i)
        {
            nVerts[i] = mesh_.faces()[i+mesh_.nInternalFaces()].size();
        }

        // Replace data on coupled patches with their neighbour ones.
        syncTools::swapBoundaryFaceList(mesh_, nVerts);

        forAll(nVerts, i)
        {
            label faceI = i+mesh_.nInternalFaces();

            const face& f = mesh_.faces()[faceI];

            if (f.size() != nVerts[i])
            {
                dumpCell(mesh_.faceOwner()[faceI]);

                label patchI = mesh_.boundaryMesh().whichPatch(faceI);

                FatalErrorInFunction
                    << "Faces do not seem to be correct across coupled"
                    << " boundaries" << endl
                    << "Coupled face " << faceI
                    << " on patch " << patchI
                    << " " << mesh_.boundaryMesh()[patchI].name()
                    << " coords:" << UIndirectList<point>(mesh_.points(), f)()
                    << " has size:" << f.size()
                    << " (coupled) neighbour face has size:"
                    << nVerts[i]
                    << abort(FatalError);
            }
        }
    }


    // Check points of face
    // ~~~~~~~~~~~~~~~~~~~~
    {
        // Anchor points.
        pointField anchorPoints(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(anchorPoints, i)
        {
            label faceI = i+mesh_.nInternalFaces();
            const point& fc = mesh_.faceCentres()[faceI];
            const face& f = mesh_.faces()[faceI];
            const vector anchorVec(mesh_.points()[f[0]] - fc);

            anchorPoints[i] = anchorVec;
        }

        // Replace data on coupled patches with their neighbour ones. Apply
        // rotation transformation (but not separation since is relative vector
        // to point on same face.
        syncTools::swapBoundaryFaceList(mesh_, anchorPoints);

        forAll(anchorPoints, i)
        {
            label faceI = i+mesh_.nInternalFaces();
            const point& fc = mesh_.faceCentres()[faceI];
            const face& f = mesh_.faces()[faceI];
            const vector anchorVec(mesh_.points()[f[0]] - fc);

            if (mag(anchorVec - anchorPoints[i]) > smallDim)
            {
                dumpCell(mesh_.faceOwner()[faceI]);

                label patchI = mesh_.boundaryMesh().whichPatch(faceI);

                FatalErrorInFunction
                    << "Faces do not seem to be correct across coupled"
                    << " boundaries" << endl
                    << "Coupled face " << faceI
                    << " on patch " << patchI
                    << " " << mesh_.boundaryMesh()[patchI].name()
                    << " coords:" << UIndirectList<point>(mesh_.points(), f)()
                    << " has anchor vector:" << anchorVec
                    << " (coupled) neighbour face anchor vector differs:"
                    << anchorPoints[i]
                    << " to within tolerance " << smallDim
                    << abort(FatalError);
            }
        }
    }

    if (debug)
    {
        Pout<< "hexTetRef8::checkMesh : Returning" << endl;
    }
}


void Foam::hexTetRef8::checkRefinementLevels
(
    const label maxPointDiff,
    const labelList& pointsToCheck
) const
{
    if (debug)
    {
        Pout<< "hexTetRef8::checkRefinementLevels :"
            << " Checking 2:1 refinement level" << endl;
    }

    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorInFunction
            << "cellLevel size should be number of cells"
            << " and pointLevel size should be number of points."<< nl
            << "cellLevel:" << cellLevel_.size()
            << " mesh.nCells():" << mesh_.nCells() << nl
            << "pointLevel:" << pointLevel_.size()
            << " mesh.nPoints():" << mesh_.nPoints()
            << abort(FatalError);
    }


    // Check 2:1 consistency.
    // ~~~~~~~~~~~~~~~~~~~~~~

    {
        // Internal faces.
        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            label own = mesh_.faceOwner()[faceI];
            label nei = mesh_.faceNeighbour()[faceI];

            if (mag(cellLevel_[own] - cellLevel_[nei]) > 1)
            {
                dumpCell(own);
                dumpCell(nei);

                FatalErrorInFunction
                    << "Celllevel does not satisfy 2:1 constraint." << nl
                    << "On face " << faceI << " owner cell " << own
                    << " has refinement " << cellLevel_[own]
                    << " neighbour cell " << nei << " has refinement "
                    << cellLevel_[nei]
                    << abort(FatalError);
            }
        }

        // Coupled faces. Get neighbouring value
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(neiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

            neiLevel[i] = cellLevel_[own];
        }

        // No separation
        syncTools::swapBoundaryFaceList(mesh_, neiLevel);

        forAll(neiLevel, i)
        {
            label faceI = i+mesh_.nInternalFaces();

            label own = mesh_.faceOwner()[faceI];

            if (mag(cellLevel_[own] - neiLevel[i]) > 1)
            {
                dumpCell(own);

                label patchI = mesh_.boundaryMesh().whichPatch(faceI);

                FatalErrorInFunction
                    << "Celllevel does not satisfy 2:1 constraint."
                    << " On coupled face " << faceI
                    << " on patch " << patchI << " "
                    << mesh_.boundaryMesh()[patchI].name()
                    << " owner cell " << own << " has refinement "
                    << cellLevel_[own]
                    << " (coupled) neighbour cell has refinement "
                    << neiLevel[i]
                    << abort(FatalError);
            }
        }
    }


    // Check pointLevel is synchronized
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    {
        labelList syncPointLevel(pointLevel_);

        // Get min level
        syncTools::syncPointList
        (
            mesh_,
            syncPointLevel,
            minEqOp<label>(),
            labelMax
        );


        forAll(syncPointLevel, pointI)
        {
            if (pointLevel_[pointI] != syncPointLevel[pointI])
            {
                FatalErrorInFunction
                    << "PointLevel is not consistent across coupled patches."
                    << endl
                    << "point:" << pointI << " coord:" << mesh_.points()[pointI]
                    << " has level " << pointLevel_[pointI]
                    << " whereas the coupled point has level "
                    << syncPointLevel[pointI]
                    << abort(FatalError);
            }
        }
    }


    // Check 2:1 across points (instead of faces)
    if (maxPointDiff != -1)
    {
        // Determine per point the max cell level.
        labelList maxPointLevel(mesh_.nPoints(), 0);

        forAll(maxPointLevel, pointI)
        {
            const labelList& pCells = mesh_.pointCells(pointI);

            label& pLevel = maxPointLevel[pointI];

            forAll(pCells, i)
            {
                pLevel = max(pLevel, cellLevel_[pCells[i]]);
            }
        }

        // Sync maxPointLevel to neighbour
        syncTools::syncPointList
        (
            mesh_,
            maxPointLevel,
            maxEqOp<label>(),
            labelMin            // null value
        );

        // Check 2:1 across boundary points
        forAll(pointsToCheck, i)
        {
            label pointI = pointsToCheck[i];

            const labelList& pCells = mesh_.pointCells(pointI);

            forAll(pCells, i)
            {
                label cellI = pCells[i];

                if
                (
                    mag(cellLevel_[cellI]-maxPointLevel[pointI])
                  > maxPointDiff
                )
                {
                    dumpCell(cellI);

                    FatalErrorInFunction
                        << "Too big a difference between"
                        << " point-connected cells." << nl
                        << "cell:" << cellI
                        << " cellLevel:" << cellLevel_[cellI]
                        << " uses point:" << pointI
                        << " coord:" << mesh_.points()[pointI]
                        << " which is also used by a cell with level:"
                        << maxPointLevel[pointI]
                        << abort(FatalError);
                }
            }
        }
    }


    //- Gives problems after first splitting off inside mesher.
    //// Hanging points
    //{
    //    // Any patches with points having only two edges.
    //
    //    boolList isHangingPoint(mesh_.nPoints(), false);
    //
    //    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    //
    //    forAll(patches, patchI)
    //    {
    //        const polyPatch& pp = patches[patchI];
    //
    //        const labelList& meshPoints = pp.meshPoints();
    //
    //        forAll(meshPoints, i)
    //        {
    //            label pointI = meshPoints[i];
    //
    //            const labelList& pEdges = mesh_.pointEdges()[pointI];
    //
    //            if (pEdges.size() == 2)
    //            {
    //                isHangingPoint[pointI] = true;
    //            }
    //        }
    //    }
    //
    //    syncTools::syncPointList
    //    (
    //        mesh_,
    //        isHangingPoint,
    //        andEqOp<bool>(),        // only if all decide it is hanging point
    //        true,                   // null
    //        false                   // no separation
    //    );
    //
    //    //OFstream str(mesh_.time().path()/"hangingPoints.obj");
    //
    //    label nHanging = 0;
    //
    //    forAll(isHangingPoint, pointI)
    //    {
    //        if (isHangingPoint[pointI])
    //        {
    //            nHanging++;
    //
    //            Pout<< "Hanging boundary point " << pointI
    //                << " at " << mesh_.points()[pointI]
    //                << endl;
    //            //meshTools::writeOBJ(str, mesh_.points()[pointI]);
    //        }
    //    }
    //
    //    if (returnReduce(nHanging, sumOp<label>()) > 0)
    //    {
    //        FatalErrorInFunction
    //            << "Detected a point used by two edges only (hanging point)"
    //            << nl << "This is not allowed"
    //            << abort(FatalError);
    //    }
    //}
}


const Foam::cellShapeList& Foam::hexTetRef8::cellShapes() const
{
    if (cellShapesPtr_.empty())
    {
        if (debug)
        {
            Pout<< "hexTetRef8::cellShapes() : calculating splitHex cellShapes."
                << " cellLevel:" << cellLevel_.size()
                << " pointLevel:" << pointLevel_.size()
                << endl;
        }

        const cellShapeList& meshShapes = mesh_.cellShapes();
        cellShapesPtr_.reset(new cellShapeList(meshShapes));

        label nSplitHex = 0;
        label nUnrecognised = 0;

        forAll(cellLevel_, cellI)
        {
            if (meshShapes[cellI].model().index() == 0)
            {
                label level = cellLevel_[cellI];

                // Return true if we've found 6 quads
                DynamicList<face> quads;
                bool haveQuads = matchHexShape
                (
                    cellI,
                    level,
                    quads
                );

                if (haveQuads)
                {
                    faceList faces(quads.xfer());
                    cellShapesPtr_()[cellI] = degenerateMatcher::match(faces);
                    nSplitHex++;
                }
                else
                {
                    nUnrecognised++;
                }
            }
        }
        if (debug)
        {
            Pout<< "hexTetRef8::cellShapes() :"
                << " nCells:" << mesh_.nCells() << " of which" << nl
                << "    primitive:" << (mesh_.nCells()-nSplitHex-nUnrecognised)
                << nl
                << "    split-hex:" << nSplitHex << nl
                << "    poly     :" << nUnrecognised << nl
                << endl;
        }
    }
    return cellShapesPtr_();
}



//
// Unrefinement
// ~~~~~~~~~~~~
//


Foam::labelList Foam::hexTetRef8::getSplitPoints() const
{
    if (debug)
    {
        checkRefinementLevels(-1, labelList(0));
    }

    if (debug)
    {
        Pout<< "hexTetRef8::getSplitPoints :"
            << " Calculating unrefineable points" << endl;
    }


    if (!history_.active())
    {
        FatalErrorInFunction
            << "Only call if constructed with history capability"
            << abort(FatalError);
    }

    // Master cell
    // -1 undetermined
    // -2 certainly not split point
    // >= label of master cell
    labelList splitMaster(mesh_.nPoints(), -1);
    labelList splitMasterLevel(mesh_.nPoints(), 0);

    // Unmark all with not 8 cells
    //const labelListList& pointCells = mesh_.pointCells();

    for (label pointI = 0; pointI < mesh_.nPoints(); pointI++)
    {
        const labelList& pCells = mesh_.pointCells(pointI);

        if (pCells.size() != 8)
        {
            splitMaster[pointI] = -2;
        }
    }

    // Unmark all with different master cells
    const labelList& visibleCells = history_.visibleCells();

    forAll(visibleCells, cellI)
    {
        const labelList& cPoints = mesh_.cellPoints(cellI);

        if (visibleCells[cellI] != -1 && history_.parentIndex(cellI) >= 0)
        {
            label parentIndex = history_.parentIndex(cellI);

            // Check same master.
            forAll(cPoints, i)
            {
                label pointI = cPoints[i];

                label masterCellI = splitMaster[pointI];

                if (masterCellI == -1)
                {
                    // First time visit of point. Store parent cell and
                    // level of the parent cell (with respect to cellI). This
                    // is additional guarantee that we're referring to the
                    // same master at the same refinement level.

                    splitMaster[pointI] = parentIndex;
                    splitMasterLevel[pointI] = cellLevel_[cellI] - 1;
                }
                else if (masterCellI == -2)
                {
                    // Already decided that point is not splitPoint
                }
                else if
                (
                    (masterCellI != parentIndex)
                 || (splitMasterLevel[pointI] != cellLevel_[cellI] - 1)
                )
                {
                    // Different masters so point is on two refinement
                    // patterns
                    splitMaster[pointI] = -2;
                }
            }
        }
        else
        {
            // Either not visible or is unrefined cell
            forAll(cPoints, i)
            {
                label pointI = cPoints[i];

                splitMaster[pointI] = -2;
            }
        }
    }

    // Unmark boundary faces
    for
    (
        label faceI = mesh_.nInternalFaces();
        faceI < mesh_.nFaces();
        faceI++
    )
    {
        const face& f = mesh_.faces()[faceI];

        forAll(f, fp)
        {
            splitMaster[f[fp]] = -2;
        }
    }


    // Collect into labelList

    label nSplitPoints = 0;

    forAll(splitMaster, pointI)
    {
        if (splitMaster[pointI] >= 0)
        {
            nSplitPoints++;
        }
    }

    labelList splitPoints(nSplitPoints);
    nSplitPoints = 0;

    forAll(splitMaster, pointI)
    {
        if (splitMaster[pointI] >= 0)
        {
            splitPoints[nSplitPoints++] = pointI;
        }
    }

    return splitPoints;
}


//void Foam::hexTetRef8::markIndex
//(
//    const label maxLevel,
//    const label level,
//    const label index,
//    const label markValue,
//    labelList& indexValues
//) const
//{
//    if (level < maxLevel && indexValues[index] == -1)
//    {
//        // Mark
//        indexValues[index] = markValue;
//
//        // Mark parent
//        const splitCell8& split = history_.splitCells()[index];
//
//        if (split.parent_ >= 0)
//        {
//            markIndex
//            (
//              maxLevel, level+1, split.parent_, markValue, indexValues);
//            )
//        }
//    }
//}
//
//
//// Get all cells which (down to level) originate from the same cell.
//// level=0 returns cell only, level=1 returns the 8 cells this cell
//// originates from, level=2 returns 64 cells etc.
//// If the cell does not originate from refinement returns just itself.
//void Foam::hexTetRef8::markCellClusters
//(
//    const label maxLevel,
//    labelList& cluster
//) const
//{
//    cluster.setSize(mesh_.nCells());
//    cluster = -1;
//
//    const DynamicList<splitCell8>& splitCells = history_.splitCells();
//
//    // Mark all splitCells down to level maxLevel with a cell originating from
//    // it.
//
//    labelList indexLevel(splitCells.size(), -1);
//
//    forAll(visibleCells, cellI)
//    {
//        label index = visibleCells[cellI];
//
//        if (index >= 0)
//        {
//            markIndex(maxLevel, 0, index, cellI, indexLevel);
//        }
//    }
//
//    // Mark cells with splitCell
//}


void Foam::hexTetRef8::consistentUnrefinement
(
    const labelList& pointsToUnrefine,
    const labelList& tetCellsToUnrefine,
    labelList& newSplitPoints,
    labelList& newTetsToUnrefine,
    const bool maxSet
) const
{
    if (debug)
    {
        Pout<< "hexTetRef8::consistentUnrefinement :"
            << " Determining 2:1 consistent unrefinement" << endl;
    }

    if (maxSet)
    {
        FatalErrorInFunction
            << "maxSet not implemented yet."
            << abort(FatalError);
    }

    const labelList &visibleCells
        = history_.visibleCells();

    const DynamicList<refinementMRAHistory::splitCell8> &splitCells
        = history_.splitCells();

    // Loop, modifying pointsToUnrefine, until no more changes to due to 2:1
    // conflicts.
    // maxSet = false : unselect points to refine
    // maxSet = true: select points to refine

    // Maintain boolList for pointsToUnrefine and cellsToUnrefine
    PackedBoolList unrefinePoint(mesh_.nPoints());

    forAll(pointsToUnrefine, i)
    {
        label pointI = pointsToUnrefine[i];

        unrefinePoint.set(pointI);
    }

    // Maintain boolList for tetCellsToUnrefine
    PackedBoolList unrefineTetCell(mesh_.nCells());

    forAll(tetCellsToUnrefine, i)
    {
        unrefineTetCell.set(tetCellsToUnrefine[i]);
    }


    while (true)
    {
        // Construct cells to unrefine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        PackedBoolList unrefineCell(mesh_.nCells());

        forAll(unrefinePoint, pointI)
        {
            if (unrefinePoint.get(pointI))
            {
                const labelList& pCells = mesh_.pointCells(pointI);

                forAll(pCells, j)
                {
                    unrefineCell.set(pCells[j]);
                }
            }
        }

        forAll(unrefineTetCell, cellI)
        {
            if (unrefineTetCell.get(cellI))
            {
                unrefineCell.set(cellI);
            }
        }


        label nChanged = 0;


        // Check 2:1 consistency taking refinement into account
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Internal faces.
        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            label own = mesh_.faceOwner()[faceI];
            label ownLevel = cellLevel_[own] - unrefineCell.get(own);

            label nei = mesh_.faceNeighbour()[faceI];
            label neiLevel = cellLevel_[nei] - unrefineCell.get(nei);

            label nOwnAnchors
                = countAnchors(mesh_.cellPoints()[own], cellLevel_[own]);

            label nNeiAnchors
                = countAnchors(mesh_.cellPoints()[nei], cellLevel_[nei]);

            bool isHexTetInterface = !(nNeiAnchors == nOwnAnchors);

            if
            (
//                ownLevel < (neiLevel-1)
                ((ownLevel < (neiLevel-1)) ||
                    ((ownLevel < neiLevel) && isHexTetInterface))
            )
            {
                // Since was 2:1 this can only occur if own is marked for
                // unrefinement.

                if (maxSet)
                {
                    unrefineCell.set(nei);
                }
                else
                {
                    // could also combine with unset:
                    // if (!unrefineCell.unset(own))
                    // {
                    //     FatalErrorInFunction
                    //         << "problem cell already unset"
                    //         << abort(FatalError);
                    // }
                    if (unrefineCell.get(own) == 0)
                    {
                        FatalErrorInFunction
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(own);
                }
                nChanged++;
            }
            else if
            (
//                neiLevel < (ownLevel-1)
                ((neiLevel < (ownLevel-1)) ||
                    ((neiLevel < ownLevel) && isHexTetInterface))
            )
            {
                if (maxSet)
                {
                    unrefineCell.set(own);
                }
                else
                {
                    if (unrefineCell.get(nei) == 0)
                    {
                        FatalErrorInFunction
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(nei);
                }
                nChanged++;
            }
        }


        // Coupled faces. Swap owner level to get neighbouring cell level.
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
        labelList nNeiAnchors(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(neiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

            neiLevel[i] = cellLevel_[own] - unrefineCell.get(own);
            nNeiAnchors[i]
                = countAnchors(mesh_.cellPoints()[own], cellLevel_[own]);
        }

        // Swap to neighbour
        syncTools::swapBoundaryFaceList(mesh_, neiLevel);
        syncTools::swapBoundaryFaceList(mesh_, nNeiAnchors);

        forAll(neiLevel, i)
        {
            label faceI = i+mesh_.nInternalFaces();
            label own = mesh_.faceOwner()[faceI];
            label ownLevel = cellLevel_[own] - unrefineCell.get(own);

            label nOwnAnchors
                = countAnchors(mesh_.cellPoints()[own], cellLevel_[own]);

            bool isHexTetInterface = !(nNeiAnchors[i] == nOwnAnchors);

            if
            (
                (ownLevel < (neiLevel[i]-1))
             || ((ownLevel < neiLevel[i]) && isHexTetInterface)
            )
            {
                if (!maxSet)
                {
                    if (unrefineCell.get(own) == 0)
                    {
                        FatalErrorInFunction
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(own);
                    nChanged++;
                }
            }
            else if
            (
                (neiLevel[i] < (ownLevel-1))
             || ((neiLevel[i] < ownLevel) && isHexTetInterface)
            )
            {
                if (maxSet)
                {
                    if (unrefineCell.get(own) == 1)
                    {
                        FatalErrorInFunction
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.set(own);
                    nChanged++;
                }
            }
        }

        // Convert cellsToUnrefine back into points/tetCell to unrefine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Knock out any point whose cell neighbour cannot be unrefined.
        forAll(unrefinePoint, pointI)
        {
            if (unrefinePoint.get(pointI))
            {
                const labelList& pCells = mesh_.pointCells(pointI);

                forAll(pCells, j)
                {
                    if (!unrefineCell.get(pCells[j]))
                    {
                        unrefinePoint.unset(pointI);
                        break;
                    }
                }
            }
        }

        // Knock out any cells whose brothers cannot be unrefined.
        forAll(unrefineTetCell, cellI)
        {
            if (unrefineTetCell.get(cellI))
            {
                label parentSplitCellID
                    = splitCells[visibleCells[cellI]].parent_;

                const FixedList<label, 8> brotherCells
                    = splitCells[parentSplitCellID].addedCellIDsPtr_();

                bool isUnrefinementCancelled = false;

                forAll(brotherCells, j)
                {
                    if (!unrefineCell.get(brotherCells[j]))
                    {
                        isUnrefinementCancelled = true;
                        break;
                    }
                }

                if (isUnrefinementCancelled)
                {
                    forAll(brotherCells, j)
                    {
                        unrefineTetCell.unset(brotherCells[j]);
                        nChanged++;
                    }
                }
            }
        }

        reduce(nChanged, sumOp<label>());

        if (debug)
        {
            Pout<< "hexTetRef8::consistentUnrefinement :"
                << " Changed " << nChanged
                << " refinement levels due to 2:1 conflicts."
                << endl;
        }

        if (nChanged == 0)
        {
            break;
        }
    }


    // Convert back to labelList.
    label nSet = 0;

    forAll(unrefinePoint, pointI)
    {
        if (unrefinePoint.get(pointI))
        {
            nSet++;
        }
    }

    newSplitPoints.resize(nSet);
    nSet = 0;

    forAll(unrefinePoint, pointI)
    {
        if (unrefinePoint.get(pointI))
        {
            newSplitPoints[nSet++] = pointI;
        }
    }

    // pick out those tet-cells from whom the others have
    // been created
    forAll(unrefineTetCell, cellI)
    {
        if (unrefineTetCell.get(cellI))
        {
            label parentSplitCellID
                = splitCells[visibleCells[cellI]].parent_;

            FixedList<label, 8> brotherCells
                = splitCells[parentSplitCellID].addedCellIDsPtr_();

            label parent = -1;
            {
                labelList broCellsForMin(brotherCells);
                parent = min(broCellsForMin);
            }

            forAll(brotherCells, j)
            {
                if (brotherCells[j] != parent)
                {
                    unrefineTetCell.unset(brotherCells[j]);
                }
            }
        }
    }

    nSet = 0;
    forAll(unrefineTetCell, cellI)
    {
        if (unrefineTetCell.get(cellI))
        {
            nSet++;
        }
    }

    newTetsToUnrefine.resize(nSet);
    nSet = 0;

    forAll(unrefineTetCell, cellI)
    {
        if (unrefineTetCell.get(cellI))
        {
            newTetsToUnrefine[nSet++] = cellI;
        }
    }
}


Foam::Map<Foam::label> Foam::hexTetRef8::setUnrefinement
(
    const labelList& splitPointLabels,
    const labelList& unrefineTetCellsParent,
    polyTopoChange& meshMod
)
{
    if (!history_.active())
    {
        FatalErrorInFunction
            << "Only call if constructed with history capability"
            << abort(FatalError);
    }

    const labelList &visibleCells
        = history_.visibleCells();

    const DynamicList<refinementMRAHistory::splitCell8> &splitCells
        = history_.splitCells();

    if (debug)
    {
        Pout<< "hexTetRef8::setUnrefinement :"
            << " Checking initial mesh just to make sure" << endl;

        checkMesh();

        forAll(cellLevel_, cellI)
        {
            if (cellLevel_[cellI] < 0)
            {
                FatalErrorInFunction
                    << "Illegal cell level " << cellLevel_[cellI]
                    << " for cell " << cellI
                    << abort(FatalError);
            }
        }


        // Write to sets.
        pointSet pSet(mesh_, "splitPoints", splitPointLabels);
        pSet.write();

        cellSet cSet(mesh_, "splitPointCells", splitPointLabels.size());

        forAll(splitPointLabels, i)
        {
            const labelList& pCells = mesh_.pointCells(splitPointLabels[i]);

            forAll(pCells, j)
            {
                cSet.insert(pCells[j]);
            }
        }
        forAll(unrefineTetCellsParent, i)
        {
            label parentSplitCellID
                = splitCells[visibleCells[unrefineTetCellsParent[i]]].parent_;

            const FixedList<label, 8>  brotherCells
                = splitCells[parentSplitCellID].addedCellIDsPtr_();

            forAll(brotherCells, j)
            {
                cSet.insert(brotherCells[j]);
            }
        }

        cSet.write();

        Pout<< "hexTetRef8::setRefinement : Dumping " << pSet.size()
            << " points and "
            << cSet.size() << " cells for unrefinement to" << nl
            << "    pointSet " << pSet.objectPath() << nl
            << "    cellSet " << cSet.objectPath()
            << endl;
    }

    // construct a map from faces to anchorFace
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Map<label> faceToCenterFace(3*unrefineTetCellsParent.size());

    forAll(unrefineTetCellsParent, i)
    {
        label parentSplitCellID
            = splitCells[visibleCells[unrefineTetCellsParent[i]]].parent_;

        const FixedList<label, 8> brotherCells
            = splitCells[parentSplitCellID].addedCellIDsPtr_();

        labelList anchorCells(4, -1);
        labelListList anchorCellAnchors(4, labelList(4, -1));

        forAll (anchorCells, j)
        {
            anchorCells[j] = brotherCells[j];
            anchorCellAnchors[j] = getAnchors(anchorCells[j]);
        }

        labelList masterAnchors(4, -1);
        labelListList masterAnchorToes(4, labelList(4, -1));

        getMasterAnchorsAndToes
        (
            anchorCellAnchors,
            cellLevel_[unrefineTetCellsParent[i]]-1,
            masterAnchors,
            masterAnchorToes
        );

        labelList subCellFaces;
        {
            std::set<label> subCellFaceSet;
            forAll(brotherCells, j)
            {
                forAll(mesh_.cells()[brotherCells[j]], k)
                {
                    subCellFaceSet.insert(mesh_.cells()[brotherCells[j]][k]);
                }
            }

            subCellFaces.resize(subCellFaceSet.size());

            std::set<label>::iterator it;
            label j = 0;
            for (it = subCellFaceSet.begin(); it != subCellFaceSet.end(); ++it)
            {
                subCellFaces[j++] = *it;
            }
        }

        labelListList centerFaceAndAnchorFaces(4, labelList(4, -1));

        getCenterFaceAndAnchorFaces
        (
            subCellFaces,
            masterAnchors,
            masterAnchorToes,
            centerFaceAndAnchorFaces
        );

        forAll(centerFaceAndAnchorFaces, j)
        {
            if
            (
                centerFaceAndAnchorFaces[j][0] != -1 &&
                centerFaceAndAnchorFaces[j][1] != -1 &&
                centerFaceAndAnchorFaces[j][2] != -1 &&
                centerFaceAndAnchorFaces[j][3] != -1
            )
            {
                faceToCenterFace.insert
                (
                    centerFaceAndAnchorFaces[j][1],
                    centerFaceAndAnchorFaces[j][0]
                );
                faceToCenterFace.insert
                (
                    centerFaceAndAnchorFaces[j][2],
                    centerFaceAndAnchorFaces[j][0]
                );
                faceToCenterFace.insert
                (
                    centerFaceAndAnchorFaces[j][3],
                    centerFaceAndAnchorFaces[j][0]
                );
            }
        }
    }

    // Find split faces
    // ~~~~~~~~~~~~~~~~

    labelList cellRegion;
    labelList cellRegionMaster;
    labelList facesToRemove;

    {
        labelHashSet splitFaces
        (
            12*splitPointLabels.size()
          + 8*unrefineTetCellsParent.size()
        );

        forAll(splitPointLabels, i)
        {
            const labelList& pFaces = mesh_.pointFaces()[splitPointLabels[i]];

            forAll(pFaces, j)
            {
                splitFaces.insert(pFaces[j]);
            }
        }
        forAll(unrefineTetCellsParent, i)
        {
            labelList pFaces = getSplitFaces(unrefineTetCellsParent[i]);

            forAll(pFaces, j)
            {
                splitFaces.insert(pFaces[j]);
            }
        }

        // Check with faceRemover what faces will get removed. Note that this
        // can be more (but never less) than splitFaces provided.
        faceRemover_.compatibleRemoves
        (
            splitFaces.toc(),   // pierced faces
            cellRegion,         // per cell -1 or region it is merged into
            cellRegionMaster,   // per region the master cell
            facesToRemove       // new faces to be removed.
        );

        if (facesToRemove.size() != splitFaces.size())
        {
            FatalErrorInFunction
                << "Ininitial set of split points to unrefine does not"
                << " seem to be consistent or not mid points of refined cells"
                << abort(FatalError);
        }
    }

    // Redo the region master so it is consistent with our master.
    // This will guarantee that the new cell (for which faceRemover uses
    // the region master) is already compatible with our refinement structure.

    forAll(splitPointLabels, i)
    {
        label pointI = splitPointLabels[i];

        // Get original cell label

        const labelList& pCells = mesh_.pointCells(pointI);

        // Check
        if (pCells.size() != 8)
        {
            FatalErrorInFunction
                << "splitPoint " << pointI
                << " should have 8 cells using it. It has " << pCells
                << abort(FatalError);
        }


        // Check that the lowest numbered pCells is the master of the region
        // (should be guaranteed by directRemoveFaces)
        //if (debug)
        {
            label masterCellI = min(pCells);

            forAll(pCells, j)
            {
                label cellI = pCells[j];

                label region = cellRegion[cellI];

                if (region == -1)
                {
                    FatalErrorInFunction
                        << "Ininitial set of split points to unrefine does not"
                        << " seem to be consistent or not mid points"
                        << " of refined cells" << nl
                        << "cell:" << cellI << " on splitPoint " << pointI
                        << " has no region to be merged into"
                        << abort(FatalError);
                }

                if (masterCellI != cellRegionMaster[region])
                {
                    FatalErrorInFunction
                        << "cell:" << cellI << " on splitPoint:" << pointI
                        << " in region " << region
                        << " has master:" << cellRegionMaster[region]
                        << " which is not the lowest numbered cell"
                        << " among the pointCells:" << pCells
                        << abort(FatalError);
                }
            }
        }
    }

    forAll(unrefineTetCellsParent, i)
    {
        // Check that the lowest numbered pCells is the master of the region
        // (should be guaranteed by directRemoveFaces)
        //if (debug)
        {
            label masterCellI = unrefineTetCellsParent[i];

            label parentSplitCellID
                = splitCells[visibleCells[unrefineTetCellsParent[i]]].parent_;

            const FixedList<label, 8> pCells
                = splitCells[parentSplitCellID].addedCellIDsPtr_();

            forAll(pCells, j)
            {
                label cellI = pCells[j];

                label region = cellRegion[cellI];

                if (region == -1)
                {
                    FatalErrorInFunction
                        << "Ininitial set of split points to unrefine does not"
                        << " seem to be consistent or not mid points"
                        << " of refined cells"
                        << abort(FatalError);
                }

                if (masterCellI != cellRegionMaster[region])
                {
                    FatalErrorInFunction
                        << "cell:" << cellI
                        << " in region " << region
                        << " has master:" << cellRegionMaster[region]
                        << " which is not the lowest numbered cell"
                        << " among the pointCells:" << pCells
                        << abort(FatalError);
                }
            }
        }
    }

    // Insert all commands to combine cells. Never fails so don't have to
    // test for success.
    faceRemover_.setRefinement
    (
        facesToRemove,
        cellRegion,
        cellRegionMaster,
        meshMod
    );

    // Remove the 8 cells that originated from merging around the split point
    // and adapt cell levels (not that pointLevels stay the same since points
    // either get removed or stay at the same position.
    forAll(splitPointLabels, i)
    {
        label pointI = splitPointLabels[i];

        const labelList& pCells = mesh_.pointCells(pointI);

        label masterCellI = min(pCells);

        forAll(pCells, j)
        {
            cellLevel_[pCells[j]]--;
        }

        history_.combineCells(masterCellI, pCells);
    }

    forAll(unrefineTetCellsParent, i)
    {
        label parentSplitCellID
            = splitCells[visibleCells[unrefineTetCellsParent[i]]].parent_;

        const FixedList<label, 8> pCells
            = splitCells[parentSplitCellID].addedCellIDsPtr_();

        labelList pCellsTransfer(pCells);

        label masterCellI = unrefineTetCellsParent[i];

        forAll(pCells, j)
        {
            cellLevel_[pCells[j]]--;
        }

        history_.combineCells(masterCellI, pCellsTransfer);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // history_.updateMesh will take care of truncating.

    return faceToCenterFace;
}


// Write refinement to polyMesh directory.
bool Foam::hexTetRef8::write() const
{
    bool writeOk =
        cellLevel_.write()
     && pointLevel_.write()
     && level0Edge_.write();

    if (history_.active())
    {
        writeOk = writeOk && history_.write();
    }

    return writeOk;
}


// ==========================
// TET Unrefinement Functions
// ==========================

//- For given parent TET cell ID, get the split faces.
//  These are faces that are INSIDE the TET (4 diag faces
//  and 4 anchor internal faces)
Foam::labelList Foam::hexTetRef8::getSplitFaces
(
    label minCell
)
{
    const labelList &visibleCells
        = history_.visibleCells();

    const DynamicList<refinementMRAHistory::splitCell8> &splitCells
        = history_.splitCells();


    // get brothers
    label parentSplitCellID
        = splitCells[visibleCells[minCell]].parent_;

    const FixedList<label, 8>  brotherCells
        = splitCells[parentSplitCellID].addedCellIDsPtr_();



    // the inner cells are the last 4 cells. Loop over
    // them, and gather only those faces which have both
    // own and nei as one of the brotherCells

    HashSet<label, Hash<label> > foundFaces;
    HashSet<label, Hash<label> > brotherCellsHash;

    forAll(brotherCells, i)
    {
        brotherCellsHash.insert(brotherCells[i]);
    }

    labelList splitFaces(8, -1);

    label nFoundFaces = 0;

    for (label i = 4; i < 8; i++)
    {
        label cellI = brotherCells[i];

        labelList faces = mesh_.cells()[cellI];

        forAll(faces, j)
        {
            label faceJ = faces[j];

            if (!foundFaces[faceJ])
            {
                label own = mesh_.faceOwner()[faceJ];
                label nei = -1;
                if (mesh_.isInternalFace(faceJ))
                {
                    nei = mesh_.faceNeighbour()[faceJ];
                }

                if (brotherCellsHash[own] && brotherCellsHash[nei])
                {
                    foundFaces.insert(faceJ);
                    splitFaces[nFoundFaces] = faceJ;
                    nFoundFaces++;
                }
            }

            if (nFoundFaces == 8)
            {
                break;
            }
        }

        if (nFoundFaces == 8)
        {
            break;
        }
    }

    if (nFoundFaces < 8)
    {
        FatalErrorInFunction
            << "Cannot find 8 splitFaces for Tet Cell " << minCell
            << abort(FatalError);
    }

    return splitFaces;
}


void Foam::hexTetRef8::getMasterAnchorsAndToes
(
    const labelListList &anchorCellAnchors,
    label masterCellLevel,
    labelList &masterAnchors,
    labelListList &masterAnchorToes
)
{
    // find the master anchors
    forAll(anchorCellAnchors, i)
    {
        forAll(anchorCellAnchors[i], j)
        {
            if (pointLevel_[anchorCellAnchors[i][j]] <= masterCellLevel)
            {
                masterAnchors[i] = anchorCellAnchors[i][j];
                break;
            }
        }
    }

    // find the master anchor-toes
    forAll(masterAnchorToes, i)
    {
        forAll(anchorCellAnchors[i], j)
        {
            if (anchorCellAnchors[i][j] != masterAnchors[i])
            {
                label mid = anchorCellAnchors[i][j];
                label matchingAnchorIndex = -1;
                bool mAIFound = false;

                forAll(anchorCellAnchors, k)
                {
                    if (k != i)
                    {
                        forAll(anchorCellAnchors[k], l)
                        {
                            if (anchorCellAnchors[k][l] == mid)
                            {
                                mAIFound = true;
                                matchingAnchorIndex = k;
                                break;
                            }
                        }

                        if (mAIFound)
                        {
                            break;
                        }
                    }
                }

                if (masterAnchorToes[i][matchingAnchorIndex] != -1)
                {
                    if (masterAnchorToes[i][matchingAnchorIndex] != mid)
                    {
                        FatalErrorInFunction
                            << "Non-matching results."
                            << abort(FatalError);
                    }
                }
                masterAnchorToes[i][matchingAnchorIndex] = mid;

                if (masterAnchorToes[matchingAnchorIndex][i] != -1)
                {
                    if (masterAnchorToes[matchingAnchorIndex][i] != mid)
                    {
                        FatalErrorInFunction
                            << "Non-matching results."
                            << abort(FatalError);
                    }
                }
                masterAnchorToes[matchingAnchorIndex][i] = mid;
            }
        }

    }
}


Foam::label Foam::hexTetRef8::findFaceID
(
    const labelHashSet &pointsToFind,
    const labelList &cellFaces
)
{
    forAll(cellFaces, i)
    {
        label nPointsMatched = 0;

        const face &f = mesh_.faces()[cellFaces[i]];

        forAll(f, j)
        {
            if (pointsToFind[f[j]])
            {
                nPointsMatched++;
            }
            if (nPointsMatched == 3)
            {
                return cellFaces[i];
            }
        }
    }

    return -1;

    FatalErrorInFunction
        << "Couldn't find the right face"
        << abort(FatalError);
}


void Foam::hexTetRef8::getCenterFaceAndAnchorFaces
(
    const labelList &cellFaces,
    const labelList &masterAnchors,
    const labelListList &masterAnchorToes,
    labelListList &centerFaceAndAnchorFaces
)
{
    // for each combination of three masterAnchors
    // find the center and anchor faces

    // 0-1-2 : store in index 0
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // find centerFace
        {
            labelHashSet pointsToFind;
            pointsToFind.insert(masterAnchorToes[0][1]);
            pointsToFind.insert(masterAnchorToes[1][2]);
            pointsToFind.insert(masterAnchorToes[2][0]);

            centerFaceAndAnchorFaces[0][0]
                = findFaceID(pointsToFind, cellFaces);
        }

        // find anchorFace for 0
        {
            labelHashSet pointsToFind;
            pointsToFind.insert(masterAnchors[0]);
            pointsToFind.insert(masterAnchorToes[0][1]);
            pointsToFind.insert(masterAnchorToes[0][2]);

            centerFaceAndAnchorFaces[0][1]
                = findFaceID(pointsToFind, cellFaces);
        }

        // find anchorFace for 1
        {
            labelHashSet pointsToFind;
            pointsToFind.insert(masterAnchors[1]);
            pointsToFind.insert(masterAnchorToes[1][0]);
            pointsToFind.insert(masterAnchorToes[1][2]);

            centerFaceAndAnchorFaces[0][2]
                = findFaceID(pointsToFind, cellFaces);
        }

        // find anchorFace for 2
        {
            labelHashSet pointsToFind;
            pointsToFind.insert(masterAnchors[2]);
            pointsToFind.insert(masterAnchorToes[2][0]);
            pointsToFind.insert(masterAnchorToes[2][1]);

            centerFaceAndAnchorFaces[0][3]
                = findFaceID(pointsToFind, cellFaces);
        }
    }

    // 0-1-3 : store in index 1
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // find centerFace
        {
            labelHashSet pointsToFind;
            pointsToFind.insert(masterAnchorToes[0][1]);
            pointsToFind.insert(masterAnchorToes[1][3]);
            pointsToFind.insert(masterAnchorToes[3][0]);

            centerFaceAndAnchorFaces[1][0]
                = findFaceID(pointsToFind, cellFaces);
        }

        // find anchorFace for 0
        {
            labelHashSet pointsToFind;
            pointsToFind.insert(masterAnchors[0]);
            pointsToFind.insert(masterAnchorToes[0][1]);
            pointsToFind.insert(masterAnchorToes[0][3]);

            centerFaceAndAnchorFaces[1][1]
                = findFaceID(pointsToFind, cellFaces);
        }

        // find anchorFace for 1
        {
            labelHashSet pointsToFind;
            pointsToFind.insert(masterAnchors[1]);
            pointsToFind.insert(masterAnchorToes[1][0]);
            pointsToFind.insert(masterAnchorToes[1][3]);

            centerFaceAndAnchorFaces[1][2]
                = findFaceID(pointsToFind, cellFaces);
        }

        // find anchorFace for 2
        {
            labelHashSet pointsToFind;
            pointsToFind.insert(masterAnchors[3]);
            pointsToFind.insert(masterAnchorToes[3][0]);
            pointsToFind.insert(masterAnchorToes[3][1]);

            centerFaceAndAnchorFaces[1][3]
                = findFaceID(pointsToFind, cellFaces);
        }
    }

    // 0-2-3 : store in index 2
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // find centerFace
        {
            labelHashSet pointsToFind;
            pointsToFind.insert(masterAnchorToes[0][2]);
            pointsToFind.insert(masterAnchorToes[2][3]);
            pointsToFind.insert(masterAnchorToes[3][0]);

            centerFaceAndAnchorFaces[2][0]
                = findFaceID(pointsToFind, cellFaces);
        }

        // find anchorFace for 0
        {
            labelHashSet pointsToFind;
            pointsToFind.insert(masterAnchors[0]);
            pointsToFind.insert(masterAnchorToes[0][2]);
            pointsToFind.insert(masterAnchorToes[0][3]);

            centerFaceAndAnchorFaces[2][1]
                = findFaceID(pointsToFind, cellFaces);
        }

        // find anchorFace for 1
        {
            labelHashSet pointsToFind;
            pointsToFind.insert(masterAnchors[2]);
            pointsToFind.insert(masterAnchorToes[2][0]);
            pointsToFind.insert(masterAnchorToes[2][3]);

            centerFaceAndAnchorFaces[2][2]
                = findFaceID(pointsToFind, cellFaces);
        }

        // find anchorFace for 2
        {
            labelHashSet pointsToFind;
            pointsToFind.insert(masterAnchors[3]);
            pointsToFind.insert(masterAnchorToes[3][0]);
            pointsToFind.insert(masterAnchorToes[3][2]);

            centerFaceAndAnchorFaces[2][3]
                = findFaceID(pointsToFind, cellFaces);
        }
    }

    // 1-2-3 : store in index 3
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // find centerFace
        {
            labelHashSet pointsToFind;
            pointsToFind.insert(masterAnchorToes[1][2]);
            pointsToFind.insert(masterAnchorToes[2][3]);
            pointsToFind.insert(masterAnchorToes[3][1]);

            centerFaceAndAnchorFaces[3][0]
                = findFaceID(pointsToFind, cellFaces);
        }

        // find anchorFace for 0
        {
            labelHashSet pointsToFind;
            pointsToFind.insert(masterAnchors[1]);
            pointsToFind.insert(masterAnchorToes[1][2]);
            pointsToFind.insert(masterAnchorToes[1][3]);

            centerFaceAndAnchorFaces[3][1]
                = findFaceID(pointsToFind, cellFaces);
        }

        // find anchorFace for 1
        {
            labelHashSet pointsToFind;
            pointsToFind.insert(masterAnchors[2]);
            pointsToFind.insert(masterAnchorToes[2][1]);
            pointsToFind.insert(masterAnchorToes[2][3]);

            centerFaceAndAnchorFaces[3][2]
                = findFaceID(pointsToFind, cellFaces);
        }

        // find anchorFace for 2
        {
            labelHashSet pointsToFind;
            pointsToFind.insert(masterAnchors[3]);
            pointsToFind.insert(masterAnchorToes[3][1]);
            pointsToFind.insert(masterAnchorToes[3][2]);

            centerFaceAndAnchorFaces[3][3]
                = findFaceID(pointsToFind, cellFaces);
        }
    }
}
// ************************************************************************* //
