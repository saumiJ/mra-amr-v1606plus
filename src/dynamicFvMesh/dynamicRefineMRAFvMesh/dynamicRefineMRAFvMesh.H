/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

Class
    Foam::dynamicRefineMRAFvMesh

Description
    A fvMesh with built-in refinement.

    Determines which cells to refine/unrefine and does all in update().


        // How often to refine
        refineInterval  1;
        // Field to be refinement on
        field           alpha.water;
        // Refine field inbetween lower..upper
        lowerRefineLevel 0.001;
        upperRefineLevel 0.999;
        // If value < unrefineLevel (default=GREAT) unrefine
        //unrefineLevel   10;
        // Have slower than 2:1 refinement
        nBufferLayers   1;
        // Refine cells only up to maxRefinement levels
        maxRefinement   2;
        // Stop refinement if maxCells reached
        maxCells        200000;
        // Flux field and corresponding velocity field. Fluxes on changed
        // faces get recalculated by interpolating the velocity. Use 'none'
        // on surfaceScalarFields that do not need to be reinterpolated, use
        // NaN to detect use of mapped variable
        correctFluxes
        (
            (phi none)  //NaN)   //none)
            (nHatf none)   //none)
            (rho*phi none)   //none)
            (ghf none)  //NaN)   //none)
        );
        // Write the refinement level as a volScalarField
        dumpLevel       true;


SourceFiles
    dynamicRefineMRAFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicRefineMRAFvMesh_H
#define dynamicRefineMRAFvMesh_H

#include "dynamicFvMesh.H"
#include "hexTetRef8.H"
#include "PackedBoolList.H"
#include "Switch.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "mapDistributePolyMesh.H"

// #include "WaveletDecomp.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class dynamicRefineMRAFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicRefineMRAFvMesh
:
    public dynamicFvMesh
{
protected:

        //- Mesh cutting engine
        hexTetRef8 meshCutter_;

        //- Dump cellLevel for postprocessing
        Switch dumpLevel_;

        //- Dump cellID for postprocessing
        Switch dumpID_;

        //- Fluxes to map
        HashTable<word> correctFluxes_;

        //- Number of refinement/unrefinement steps done so far.
        label nRefinementIterations_;

        //- Protected cells (usually since not hexes)
        PackedBoolList protectedCell_;

        //- cell treatment flag field
        PackedBoolList isCellTreated_;

        //- did cell need refinement due to MRA?
        PackedBoolList isCellMRADetailed_;

        //- max number of cells allowed
        label maxCells_;

        //- max level of refinement allowed
        label maxRefinement_;

        //- is mesh intial refinement done?
        bool isInitialized_;

        //- global accuracy threshold for unrefinement
        double globalAccuracyThreshold_;

        //- normalizing parameter for detail
        double normalizingParameter_;

        //- fraction of local threshold at which unrefinement is triggered
        double unrefTriggerFrac_;

        //- level upto whose FV accuracy we need refinement
        label relativeAccuracyLevel_;


    // Private Member Functions

        //- find the splitCellID of root of a cell
        label topParentID
        (
            label p
        );

        //- Find a list of cells that share with "cell"
        //  the ancestor "levelDiff" levels above
        void findRelatedCells
        (
            label cell,
            label levelDiff,
            DynamicList<label> &relatedCells,
            PackedBoolList &untouchableCells
        );

        //- Compute the cell average at the requested difference
        //  in level from the input cell
        double computeAverage
        (
            label cellI,
            label refLevel
        );

        //- For the given cell, determine its fate
        //  (refine? unrefine? leave untouched?) based on the vFld
        void treatCell
        (
            label cellI,
            const volScalarField &vFld,
            PackedBoolList &refineCell,
            PackedBoolList &unrefineCell
        );

        //- Conduct refinement
        bool conductRefinement
        (
            PackedBoolList &refineCell,
            PackedBoolList &unrefineCell
        );

        //- Conduct unrefinement
        bool conductUnrefinement
        (
            PackedBoolList &unrefineCell
        );

        //- Count set/unset elements in packedlist.
        static label count
        (
            const PackedBoolList&,
            const unsigned int
        );

        //- Read the projection parameters from dictionary
        void readDict();

        //- Decompose polyhedra into tetrahedra
        autoPtr<mapPolyMesh> decomposePolyhedra();

        //- Refine HEX and TET cells that are passed in via the
        //  lists. Update mesh and fields.
        autoPtr<mapPolyMesh> refine
        (
            const labelList&
        );

        //- Unrefine cells. Gets passed in centre points of HEX cells
        //  to combine and IDs of tetCells to combine
        autoPtr<mapPolyMesh> unrefine
        (
            const labelList&,
            const labelList&
        );

        //- load balance the mesh
        void balanceMesh
        (
            const Foam::floatScalar allowableImbalance
        );

        // Selection of cells to un/refine

            //- Calculates approximate value for refinement level so
            //  we don't go above maxCell
            scalar getRefineLevel
            (
                const label maxCells,
                const label maxRefinement,
                const scalar refineLevel,
                const scalarField&
            ) const;

            //- Get per cell max of connected point
            scalarField maxPointField
            (
                const scalarField&
            ) const;

            //- Get point max of connected cell
            scalarField maxCellField
            (
                const volScalarField&
            ) const;

            scalarField cellToPoint
            (
                const scalarField& vFld
            ) const;

            scalarField error
            (
                const scalarField& fld,
                const scalar minLevel,
                const scalar maxLevel
            ) const;

            //- Select candidate cells for refinement
            virtual void selectRefineCandidates
            (
                const scalar lowerRefineLevel,
                const scalar upperRefineLevel,
                const scalarField& vFld,
                PackedBoolList& candidateCell
            ) const;

            //- Subset candidate cells for refinement
            virtual labelList selectRefineCells
            (
                const label maxCells,
                const label maxRefinement,
                const PackedBoolList& candidateCell,
                const bool maxSet = true
            ) const;

            //- Extend markedCell with cell-face-cell.
            void extendMarkedCells
            (
                PackedBoolList& markedCell
            ) const;

            //- Check all cells have 8 anchor points
            void checkEightAnchorPoints
            (
                PackedBoolList& protectedCell,
                label& nProtected
            ) const;

            //- determine list of protected cells
            label identifyProtectedCells();

            //- Check if cell has a parent
            bool isCellTreatable
            (
                label cellI
            );

        // Multiphase initialization

            //- Initialize phase by refining to capture interface
            void initializeSetFields();

private:

        //- Disallow default bitwise copy construct
        dynamicRefineMRAFvMesh(const dynamicRefineMRAFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicRefineMRAFvMesh&);

public:

    //- Runtime type information
    TypeName("dynamicRefineMRAFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit dynamicRefineMRAFvMesh(const IOobject& io);


    //- Destructor
    virtual ~dynamicRefineMRAFvMesh();


    // Member Functions

        //- Direct access to the refinement engine
        const hexTetRef8& meshCutter() const
        {
            return meshCutter_;
        }

        //- Cells which should not be refined/unrefined
        const PackedBoolList& protectedCell() const
        {
            return protectedCell_;
        }

        //- Cells which should not be refined/unrefined
        PackedBoolList& protectedCell()
        {
            return protectedCell_;
        }

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Template to update all volField boundaries
        template<class Type> void correctBoundaries();


    // Writing

        //- Write using given format, version and compression
        virtual bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dynamicRefineMRAFvMeshTemplates.C"
#endif

#endif

// ************************************************************************* //
