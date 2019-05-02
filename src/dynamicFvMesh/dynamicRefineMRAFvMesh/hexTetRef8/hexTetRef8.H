/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::hexTetRef8

Description

    Breaks up polyhedral cells into tetrahedra, and then refines/unrefines
    both hexes and tets using polyTopoChange

SourceFiles
    hexTetRef8.C

\*---------------------------------------------------------------------------*/

#ifndef hexTetRef8_H
#define hexTetRef8_H

#include "labelIOList.H"
#include "face.H"
#include "HashSet.H"
#include "DynamicList.H"
#include "primitivePatch.H"
#include "removeFaces.H"
#include "refinementMRAHistory.H"
#include "PackedBoolList.H"
#include "uniformDimensionedFields.H"
#include "cellShapeList.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyPatch;
class polyTopoChange;
class mapPolyMesh;
class mapDistributePolyMesh;

/*---------------------------------------------------------------------------*\
                           Class hexTetRef8 Declaration
\*---------------------------------------------------------------------------*/

class hexTetRef8
{
    // Private data

        //- Reference to underlying mesh.
        const polyMesh& mesh_;

        //- Per cell the refinement level
        labelIOList cellLevel_;

        //- Per point the refinement level
        labelIOList pointLevel_;

        //- Typical edge length between unrefined points
        uniformDimensionedScalarField level0Edge_;

        //- Refinement history
        refinementMRAHistory history_;

        //- Face remover engine
        removeFaces faceRemover_;

        //- Level of saved points
        Map<label> savedPointLevel_;

        //- Level of saved cells
        Map<label> savedCellLevel_;

        //- Cell shapes when seen as split hexes
        mutable autoPtr<cellShapeList> cellShapesPtr_;

        //- Flag to turn on tetRefTesting
        bool tetTesting_;


    // Private Member Functions

        //- Reorder according to map.
        static void reorder
        (
            const labelList& map,
            const label len,
            const label null,
            labelList& elems
        );

        //- Get patch and zone info
        void getFaceInfo
        (
            const label faceI,
            label& patchID,
            label& zoneID,
            label& zoneFlip
        ) const;

//        //- Get zone info for new faces
//        void getNewFaceInfo
//        (
//            polyTopoChange& meshMod,
//            const label faceI,
//            label& patchID,
//            label& zoneID,
//            label& zoneFlip
//        ) const;

        //- Adds a face on top of existing faceI. Reverses if necessary.
        label addFace
        (
            polyTopoChange& meshMod,
            const label faceI,
            const face& newFace,
            const label own,
            const label nei
        ) const;

        //- Adds internal face from point. No checks on reversal.
        label addInternalFace
        (
            polyTopoChange& meshMod,
            const label meshFaceI,
            const label meshPointI,
            const face& newFace,
            const label own,
            const label nei
        ) const;

        //- Modifies existing faceI for either new owner/neighbour or new face
        //  points. Reverses if necessary.
        void modFace
        (
            polyTopoChange& meshMod,
            const label faceI,
            const face& newFace,
            const label own,
            const label nei
        ) const;

        // Modifies new faceI for either new owner/neighbour or new face points.
        void modNewFace
        (
            polyTopoChange& meshMod,
            const label faceI,
            const label parentFace,
            const face& newFace,
            const label own,
            const label nei
        ) const;

        scalar getLevel0EdgeLength() const;

        //- Get cell added to point of cellI (if any)
        label getAnchorCell
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const label cellI,
            const label faceI,
            const label pointI
        ) const;

        //- Get new owner and neighbour (in unspecified order) of pointI
        //  on faceI.
        void getFaceNeighbours
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const label faceI,
            const label pointI,

            label& own,
            label& nei
        ) const;


        //- Get index of minimum pointlevel.
        label findMinLevel(const labelList& f) const;
        //- Get maximum pointlevel.
        label findMaxLevel(const labelList& f) const;
        //- Count number of vertices <= anchorLevel
        label countAnchors(const labelList&, const label) const;
        label countAnchors(const labelList&, const label, const DynamicList<label>&) const;
        //- Debugging: dump cell as .obj file
        void dumpCell(const label cellI) const;
        //- Find index of point with wantedLevel, starting from fp.
        label findLevel
        (
            const label faceI,
            const face& f,
            const label startFp,
            const bool searchForward,
            const label wantedLevel
        ) const;

        ////- Print levels of list of points.
        //void printLevels(Ostream&, const labelList&) const;

        //- debug:check orientation of added internal face
        static void checkInternalOrientation
        (
            polyTopoChange& meshMod,
            const label cellI,
            const label faceI,
            const point& ownPt,
            const point& neiPt,
            const face& newFace
        );

        //- debug:check orientation of new boundary face
        static void checkBoundaryOrientation
        (
            polyTopoChange& meshMod,
            const label cellI,
            const label faceI,
            const point& ownPt,
            const point& boundaryPt,
            const face& newFace
        );

        //- If p0 and p1 are existing vertices check if edge is split and insert
        // splitPoint.
        void insertEdgeSplit
        (
            const labelList& edgeMidPoint,
            const label p0,
            const label p1,
            DynamicList<label>& verts
        ) const;

        //- Store in maps correspondence from midpoint to anchors and faces.
        label storeMidPointInfo
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const labelList& cellMidPoint,
            const labelList& edgeMidPoint,
            const label cellI,
            const label faceI,
            const bool faceOrder,
            const label midPointI,
            const label anchorPointI,
            const label faceMidPointI,

            Map<edge>& midPointToAnchors,
            Map<edge>& midPointToFaceMids,
            polyTopoChange& meshMod
        ) const;

        //- Create all internal faces from an unsplit face.
        void createInternalFromSplitFace
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const labelList& cellMidPoint,
            const labelList& faceMidPoint,
            const labelList& edgeMidPoint,
            const label cellI,
            const label faceI,

            Map<edge>& midPointToAnchors,
            Map<edge>& midPointToFaceMids,
            polyTopoChange& meshMod,
            label& nFacesAdded
        ) const;

        //- Get the relevant anchor for a tet-face on a hex-cell
        void getRelevantAnchorAndMid
        (
            label &relevantAnchor,
            label &relevantMid,
            const label &faceI,
            const cell &faceList,
            const label &referenceLevel,
            const DynamicList<label> &newPointLevel,
            const labelList &edgeMidPoint,
            const polyTopoChange &meshMod,
            const floatScalar &maxAngleBetweenFaces = 0.0174533 // 1 degree
        ) const;

        //- Create all internal faces to split cellI into 8.
        void createInternalFaces
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
        ) const;

        //- Store vertices from startFp upto face split point.
        //  Used when splitting face into 4.
        void walkFaceToMid
        (
            const labelList& edgeMidPoint,
            const label cLevel,
            const label faceI,
            const label startFp,
            DynamicList<label>& faceVerts
        ) const;

        //- Same as walkFaceToMid but now walk back.
        void walkFaceFromMid
        (
            const labelList& edgeMidPoint,
            const label cLevel,
            const label faceI,
            const label startFp,
            DynamicList<label>& faceVerts
        ) const;

        //- Updates refineCell so consistent 2:1 refinement. Returns local
        //  number of cells changed.
        label faceConsistentRefinement
        (
            const bool maxSet,
            PackedBoolList& refineCell
        ) const;

        //- Check wanted refinement for 2:1 consistency
        void checkWantedRefinementLevels(const labelList&) const;


        // Cellshape recognition

            //- Collect all points on face of certain level
            void collectLevelPoints
            (
                const labelList& f,
                const label level,
                DynamicList<label>& points
            ) const;

            //- Collect all points on face (in local numbering) of certain level
            void collectLevelPoints
            (
                const labelList& meshPoints,
                const labelList& f,
                const label level,
                DynamicList<label>& points
            ) const;

            //- Collect all faces with four corner points and return true if
            //  hex was matched (6 faces of each four corner points)
            bool matchHexShape
            (
                const label cellI,
                const label cellLevel,
                DynamicList<face>& quads
            ) const;


        //- Disallow default bitwise copy construct
        hexTetRef8(const hexTetRef8&);

        //- Disallow default bitwise assignment
        void operator=(const hexTetRef8&);


        // Decompose Polyhedra into Tetrahedra Helper functions

            //- For the triangulated cell surface, find the faceID
            //  of the face that shares edge e with face ownFace
            //  and return the cellID of that face's cell
            label edgeFaceNeighbour
            (
                edge e,
                label ownFace,
                const DynamicList<DynamicList<label> > &triListList,
                const DynamicList<DynamicList<label> > &newCellListList,
                const polyTopoChange &meshMod
            );

            //- Break up cell faces into triangles and add them to
            //  the refinement engine. Also assigns cells to these
            //  faces.
            void triangulateCellFaces
            (
                label                            cellI,
                PackedBoolList                   &isFaceSplit,
                DynamicList<DynamicList<label> > &triListListLocal,
                DynamicList<DynamicList<label> > &triListListGlobal,
                labelList                        &indexToSplitFaces,
                DynamicList<DynamicList<label> > &newCellListListLocal,
                DynamicList<label>               &newCellLevel,
                Foam::polyTopoChange             &meshMod
            );

            //- For the decomposed cellI faces, add internal faces
            //  between edges and cell center.
            void addInternalFacesForTets
            (
                label                            cellI,
                label                            cellICenter,
                DynamicList<DynamicList<label> > &triListListLocal,
                DynamicList<DynamicList<label> > &newCellListListLocal,
                DynamicList<DynamicList<label> > &triListListGlobal,
                labelList                        &indexToSplitFaces,
                polyTopoChange                   &meshMod
            );


        // TET Refinement Helper Functions (see Public functions for main func)

            //- get anchor points for a cell
            labelList getAnchors
            (
                label cell
            );

            //- Find the mid point between given anchors
            label pointBetweenAnchors
            (
                Pair<label> anchorPair,
                const labelList &anchors,
                label cell,
                const labelList &edgeMidPoint
            );

            //- Find anchor toes: these are those points (new or old) that,
            //  for each anchor, are the midPoints of the edges that connect
            //  that anchor to other anchors.
            labelListList getAnchorToes
            (
                label cellI,
                const labelList& anchors,
                const labelList& edgeMidPoint
            );

            // Gets cell level such that the face has four points <= level.
            label faceLevelTet
            (
                const label faceI
            ) const;

            //- Find if point is anchor. If yes, get anchorIndex.
            //  If not, return -1
            label whichAnchor
            (
                label pointID,
                const labelList &anchors
            ) const;

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
            label whichCenterFace
            (
                const face &f,
                const labelListList &anchorToes,
                label referenceLevel,
                const DynamicList<label> &newPointLevel
            );

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
            label getSubCellID
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
            );

            //- for given face, find its anchors and the mid points
            //  between them
            void getFaceAnchorsAndMids
            (
                label faceID,
                labelList &fAnchors,
                labelList &fMids,
                label fLevel,
                const DynamicList<label> &newPointLevel,
                const polyTopoChange &meshMod,
                const bool midsWanted = true
            ) const;

            //- get the right own/nei
            void getNeighbours
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
            );

            //- Compute the cartesian distance between
            //  two points
            doubleScalar cartDistance
            (
                point p1,
                point p2
            );

            //- Choose which two opposite mids to connect
            //  such that tets are of best quality. Ref
            //  Bey 1995
            Pair<label> chooseDiagonal
            (
                const labelListList &anchorToes,
                polyTopoChange &meshMod
            );

            //- From the list of centerFaces, return the
            //  indices of those that contain mid
            void getCenterFacesWithPoint
            (
                const faceList &centerFaces,
                label mid,
                const labelList &cellIAddedCells,
                faceList &mCF,
                labelList &mCFIndices,
                labelList &mCFCells
            );

            //- Given two anchorToes, find if there is a point between
            //  them in the mesh. If yes, return them in the right order
            label walkBetweenPoints
            (
                label p1,
                label p2,
                const faceList &searchFaceList,
                const DynamicList<label> &newPointLevel,
                label referenceLevel,
                DynamicList<label> &newFace,
                bool isOrderImportantForSearch,
                bool isOrderImportantForAddition
            );

            face getNewInternalFace
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
            );

            //- Return diagonal point that is NOT
            //  present in face
            label findOtherDiagPoint
            (
                Pair<label> diag,
                face f
            );

            //- Find the right midCenterFace for adding internal
            //  anchor face for given anchor
            label findAnchorMCF
            (
                const faceList& mCF,
                const labelList& anchorToes,
                Pair<label> diag
            );

            //- Find the cell adjacent to the given face of the
            //  given cell in the mesh. Returns -1 if boundary cell
            label getAdjacentCell
            (
                label face,
                label cell
            );

            //- For given parent TET cell ID, get the split faces.
            //  These are faces that are INSIDE the TET (4 diag faces
            //  and 4 anchor internal faces)
            labelList getSplitFaces
            (
                label minCell
            );

            //- For a given cell, get its centerFaces in the right
            //  order.
            //      -> b/w anchors 0-1-2
            //      -> b/w anchors 0-1-3
            //      -> b/w anchors 0-2-3
            //      -> b/w anchors 1-2-3
            faceList getCenterFaces
            (
                const labelListList &anchorToes,
                const faceList &cellFaces,
                const DynamicList<label> &newPointLevel,
                const label referenceLevel
            );

public:

    //- Runtime type information
    ClassName("hexTetRef8");


    // Constructors

        //- Construct from mesh, read_if_present refinement data
        //  (from write below). If readHistory is true does read_if_present
        //  of refinement history. If false clears all history
        hexTetRef8(const polyMesh& mesh, const bool readHistory = true);

        //- Construct from mesh and un/refinement data and optional size of
        //  starting cells
        hexTetRef8
        (
            const polyMesh& mesh,
            const labelList& cellLevel,
            const labelList& pointLevel,
            const refinementMRAHistory& history,
            const scalar level0Edge = -GREAT
        );

        //- Construct from mesh and refinement data and optional size of
        //  starting cells
        hexTetRef8
        (
            const polyMesh& mesh,
            const labelList& cellLevel,
            const labelList& pointLevel,
            const scalar level0Edge = -GREAT
        );


    // Member Functions

        // Access

            const polyMesh& mesh() const
            {
                return mesh_;
            }

            const labelIOList& cellLevel() const
            {
                return cellLevel_;
            }

            const labelIOList& pointLevel() const
            {
                return pointLevel_;
            }

            const refinementMRAHistory& history() const
            {
                return history_;
            }

            //- Typical edge length between unrefined points
            scalar level0EdgeLength() const
            {
                return level0Edge_.value();
            }

        // Refinement history

            //- Fill cellAvg of all history cells
            void fillHistory
            (
                const volScalarField& vFld
            );

        // Decompose Polyhedra into Tetrahedra

            //- Break up polyhedra into tetrahedra
            void decomposePolyhedra
            (
                polyTopoChange &meshMod,
                PackedBoolList &isNonHex

            );

        // Refinement

            //- Gets level such that the face has four points <= level.
            label faceLevel(const label faceI) const;

            //- Given valid mesh and current cell level and proposed
            //  cells to refine calculate any clashes (due to 2:1) and return
            //  ok list of cells to refine.
            //  Either adds cells to refine to set (maxSet = true) or
            //  removes cells to refine (maxSet = false)
            labelList consistentRefinement
            (
                const labelList& cellsToRefine,
                const bool maxSet
            ) const;

            //- Like consistentRefinement but slower:
            //
            //  - specify number of cells between consecutive refinement levels
            //    (consistentRefinement equivalent to 1)
            //  - specify max level difference between point-connected cells.
            //    (-1 to disable) Note that with normal 2:1 limitation
            //    (maxFaceDiff=1) there can be 8:1 size difference across point
            //    connected cells so maxPointDiff allows you to make that less.
            //  cellsToRefine : cells we're thinking about refining. It will
            //                  extend this set. All refinement levels will be
            //                  at least maxFaceDiff layers thick.
            //  facesToCheck  : additional faces where to implement the
            //                  maxFaceDiff thickness (usually only boundary
            //                  faces)
            labelList consistentSlowRefinement
            (
                const label maxFaceDiff,
                const labelList& cellsToRefine,
                const labelList& facesToCheck,
                const label maxPointDiff,
                const labelList& pointsToCheck
            ) const;

            //- Like consistentSlowRefinement but uses different meshWave
            //  (proper distance instead of toplogical count). No point checks
            //  yet.
            labelList consistentSlowRefinement2
            (
                const label maxFaceDiff,
                const labelList& cellsToRefine,
                const labelList& facesToCheck
            ) const;

            //- creates new cells and calls the setHexRefinement and setTetRefinement
            //  functions to set the mesh connectivity changes
            //  All selected cells will be split into 8.
            //  Returns per element in cells the 8 cells they were split into.
            //  Guarantees that the 0th element is the original cell label.
            labelListList setRefinement
            (
                const labelList& cellsToRefine,
                polyTopoChange& meshMod
            );

            //- Insert refinement for hex cells.
            //  Mapping:
            //  -split cells: 7 new ones get added from original
            //  -split faces: original gets modified; 3 new ones get added
            //               from original
            //  -added internal faces: added from original cell face(if
            //   that was internal) or created out-of-nothing (so will not
            //   get mapped!). Note: could make this inflate from point but
            //   that will allocate interpolation.
            //  -points added to split edge: added from edge start()
            //  -midpoints added: added from cellPoints[0].
            void setHexRefinement
            (
                const labelList& hexCellLabels,
                boolList& faceNeedsSplit,
                const labelList& newEdgeMids,
                const labelListList& cellAddedCells,
                DynamicList<label>& newPointLevel,
                polyTopoChange& meshMod
            );

            //- Insert refinement for tet cells.
            //  Mapping:
            //  -split cells: 7 new ones get added from original
            //  -split faces: original gets modified; 3 new ones get added
            //               from original
            //  -added internal faces: created out-of-nothing (so will not
            //   get mapped!). Note: could make this inflate from point but
            //   that will allocate interpolation.
            //  -points added to split edge: added from edge start()
            void setTetRefinement
            (
                const labelList& tetCellLabels,
                const boolList& cellNeedsSplit,
                const boolList& faceNeedsSplit,
                const labelList& newEdgeMids,
                const labelListList& cellAddedCells,
                const DynamicList<label>& newPointLevel,
                polyTopoChange& meshMod
            );

            //- Update local numbering for changed mesh.
            void updateMesh(const mapPolyMesh&);


            // Restoring : is where other processes delete and reinsert data.
            //             These callbacks allow this to restore the cellLevel
            //             and pointLevel for reintroduced points.
            //             Is not related to undoing my refinement

                //- Signal points/face/cells for which to store data
                void storeData
                (
                    const labelList& pointsToStore,
                    const labelList& facesToStore,
                    const labelList& cellsToStore
                );

                //- Update local numbering + undo
                //  Data to restore given as new pointlabel + stored pointlabel
                //  (i.e. what was in pointsToStore)
                void updateMesh
                (
                    const mapPolyMesh&,
                    const Map<label>& pointsToRestore,
                    const Map<label>& facesToRestore,
                    const Map<label>& cellsToRestore
                );


            //- Update local numbering for subsetted mesh.
            //  Gets new-to-old maps. Not compatible with unrefinement.
            void subset
            (
                const labelList& pointMap,
                const labelList& faceMap,
                const labelList& cellMap
            );

            //- Update local numbering for mesh redistribution
            void distribute(const mapDistributePolyMesh&);

            //- Debug: Check coupled mesh for correctness
            void checkMesh() const;

            //- Debug: Check 2:1 consistency across faces.
            //  maxPointDiff==-1 : only check 2:1 across faces
            //  maxPointDiff!=-1 : check point-connected cells.
            void checkRefinementLevels
            (
                const label maxPointDiff,
                const labelList& pointsToCheck
            ) const;

            //- Utility: get hexes as cell shapes
            const cellShapeList& cellShapes() const;


        // Unrefinement (undoing refinement, not arbitrary coarsening)

            //- Return the points at the centre of top-level split cells
            //  that can be unsplit.
            labelList getSplitPoints() const;

            //- Given proposed
            //  splitPoints to unrefine according to calculate any clashes
            //  (due to 2:1) and return ok list of points to unrefine.
            //  Either adds points to refine to set (maxSet = true) or
            //  removes points to refine (maxSet = false)
            void consistentUnrefinement
            (
                const labelList& pointsToUnrefine,
                const labelList& tetCellsToUnrefine,
                labelList& newSplitPoints,
                labelList& newTetsToUnrefine,
                const bool maxSet
            ) const;

            //- Remove some refinement. Needs to be supplied output of
            //  consistentUnrefinement. Only call if undoable set.
            //  All 8 pointCells of a split point will be combined into
            //  the lowest numbered cell of those 8.
            Map<label> setUnrefinement
            (
                const labelList& splitPointLabels,
                const labelList& unrefineTetCellsParent,
                polyTopoChange&
            );


        // Tet Unrefinement

            //- When unrefining a Tet cell, this function returns the
            //  anchors and toes for the master tet cell into which
            //  the 8 subcells will be unrefined - used in constructing
            //  a face-to-face map (see setUnrefinement)
            void getMasterAnchorsAndToes
            (
                const labelListList &anchorCellAnchors,
                label masterCellLevel,
                labelList &masterAnchors,
                labelListList &masterAnchorToes
            );

            //- Find the FaceID of a given set of points from a list of
            //  face labels
            label findFaceID
            (
                const labelHashSet &pointsToFind,
                const labelList &cellFaces
            );

            //- Return the center and anchor faces from a list of cell
            //  faces. These are then used in tet-unrefinement
            void getCenterFaceAndAnchorFaces
            (
                const labelList &cellFaces,
                const labelList &masterAnchors,
                const labelListList &masterAnchorToes,
                labelListList &centerFaceAndAnchorFaces
            );


        // Write

            // Set instance for mesh files
            void setInstance(const fileName& inst);

            //- Force writing refinement+history to polyMesh directory.
            bool write() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
