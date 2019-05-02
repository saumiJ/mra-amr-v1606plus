/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

#include "refinementMRAHistory.H"
#include "mapPolyMesh.H"
#include "mapDistributePolyMesh.H"
#include "polyMesh.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(refinementMRAHistory, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct null
Foam::refinementMRAHistory::splitCell8::splitCell8()
:
    parent_(-1),
    cellAvg_(0.0),
    addedCellIDsPtr_(NULL)
{}


//- Construct as child element of parent
Foam::refinementMRAHistory::splitCell8::splitCell8(const label parent)
:
    parent_(parent),
    cellAvg_(0.0),
    addedCellIDsPtr_(NULL)
{}


//- Construct from Istream
Foam::refinementMRAHistory::splitCell8::splitCell8(Istream& is)
{
    is >> *this;
}


//- Construct as (deep) copy.
Foam::refinementMRAHistory::splitCell8::splitCell8(const splitCell8& sc)
:
    parent_(sc.parent_),
    cellAvg_(sc.cellAvg_),
    addedCellIDsPtr_
    (
        sc.addedCellIDsPtr_.valid()
      ? new FixedList<label, 8>(sc.addedCellIDsPtr_())
      : NULL
    )
{}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

//- Copy operator since autoPtr otherwise 'steals' storage.
void Foam::refinementMRAHistory::splitCell8::operator=(const splitCell8& s)
{
    // Check for assignment to self
    if (this == &s)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    parent_ = s.parent_;

    cellAvg_ = s.cellAvg_;

    addedCellIDsPtr_.reset
    (
        s.addedCellIDsPtr_.valid()
      ? new FixedList<label, 8>(s.addedCellIDsPtr_())
      : NULL
    );
}


bool Foam::refinementMRAHistory::splitCell8::operator==(const splitCell8& s) const
{
    if (addedCellIDsPtr_.valid() != s.addedCellIDsPtr_.valid())
    {
        return false;
    }
    else if (parent_ != s.parent_)
    {
        return false;
    }
    else if (cellAvg_ != s.cellAvg_)
    {
        return false;
    }
    else if (addedCellIDsPtr_.valid())
    {
        return addedCellIDsPtr_() == s.addedCellIDsPtr_();
    }
    else
    {
        return true;
    }
}


bool Foam::refinementMRAHistory::splitCell8::operator!=(const splitCell8& s) const
{
    return !operator==(s);
}


// * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, refinementMRAHistory::splitCell8& sc)
{
    labelList addedCellIDs;

    is >> sc.parent_ >> sc.cellAvg_ >> addedCellIDs;

    if (addedCellIDs.size())
    {
        sc.addedCellIDsPtr_.reset(new FixedList<label, 8>(addedCellIDs));
    }
    else
    {
        sc.addedCellIDsPtr_.reset(NULL);
    }

    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const refinementMRAHistory::splitCell8& sc
)
{
    // Output as labelList so we can have 0 sized lists. Alternative is to
    // output as fixedlist with e.g. -1 elements and check for this upon
    // reading. However would cause much more data to be transferred.

    if (sc.addedCellIDsPtr_.valid())
    {
        return os
            << sc.parent_
            << token::SPACE
            << sc.cellAvg_
            << token::SPACE
            << labelList(sc.addedCellIDsPtr_());
    }
    else
    {
        return os << sc.parent_ << token::SPACE << sc.cellAvg_ << token::SPACE << labelList(0);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::refinementMRAHistory::writeEntry
(
    const List<splitCell8>& splitCells,
    const splitCell8& split
)
{
    // Write me:
    if (split.addedCellIDsPtr_.valid())
    {
        Pout<< "parent:" << split.parent_ << " cellAvg: " << split.cellAvg_
            << " subCellIDs:" << split.addedCellIDsPtr_()
            << endl;
    }
    else
    {
        Pout<< "parent:" << split.parent_ << " cellAvg: " << split.cellAvg_
            << " no subcells"
            << endl;
    }

    if (split.parent_ >= 0)
    {
        Pout<< "parent data:" << endl;
        // Write my parent
        string oldPrefix = Pout.prefix();
        Pout.prefix() = "  " + oldPrefix;
        writeEntry(splitCells, splitCells[split.parent_]);
        Pout.prefix() = oldPrefix;
    }
}


void Foam::refinementMRAHistory::writeDebug
(
    const labelList& visibleCells,
    const List<splitCell8>& splitCells
)
{
    string oldPrefix = Pout.prefix();
    Pout.prefix() = "";

    forAll(visibleCells, cellI)
    {
        label index = visibleCells[cellI];

        if (index >= 0)
        {
            Pout<< "Cell from refinement:" << cellI << " index:" << index
                << endl;

            string oldPrefix = Pout.prefix();
            Pout.prefix() = "  " + oldPrefix;
            writeEntry(splitCells, splitCells[index]);
            Pout.prefix() = oldPrefix;
        }
        else
        {
            Pout<< "Unrefined cell:" << cellI << " index:" << index << endl;
        }
    }
    Pout.prefix() = oldPrefix;
}


void Foam::refinementMRAHistory::checkIndices() const
{
    // Check indices.
    forAll(visibleCells_, i)
    {
        // the change here to -1 is because there might be cases where the cell is not
        // split, and hence will not have an entry in splitCells
        if (visibleCells_[i] < -1 || visibleCells_[i] >= splitCells_.size())
        {
            FatalErrorInFunction
                << "Illegal entry " << visibleCells_[i]
                << " in visibleCells at location" << i << nl
                << "It points outside the range of splitCells : 0.."
                << splitCells_.size()-1
                << abort(FatalError);
        }
    }
}


Foam::label Foam::refinementMRAHistory::allocateSplitCell
(
    const label parent,
    const label myID,
    const label i
)
{
    label index = -1;

    if (freeSplitCells_.size())
    {
        index = freeSplitCells_.remove();

        splitCells_[index] = splitCell8(parent);
    }
    else
    {
        index = splitCells_.size();

        splitCells_.append(splitCell8(parent));
    }


    // Update the parent field
    if (parent >= 0)
    {
        splitCell8& parentSplit = splitCells_[parent];

        if (parentSplit.addedCellIDsPtr_.empty())
        {
            // Allocate storage on parent for the 8 subcells.
            parentSplit.addedCellIDsPtr_.reset(new FixedList<label, 8>(-1));
        }


        // Store me on my parent
        FixedList<label, 8>& parentSplitIDs = parentSplit.addedCellIDsPtr_();

        parentSplitIDs[i] = myID;
    }

    return index;
}


void Foam::refinementMRAHistory::freeSplitCell(const label index)
{
    splitCell8& split = splitCells_[index];

    // Make sure parent does not point to me anymore.
    if (split.parent_ >= 0)
    {
        autoPtr<FixedList<label, 8> >& subCellIDsPtr =
            splitCells_[split.parent_].addedCellIDsPtr_;

        if (subCellIDsPtr.valid())
        {
            FixedList<label, 8> subCells = getAddedCellsIndices(split.parent_);

            label myPos = findIndex(subCells, index);

            if (myPos == -1)
            {
                FatalErrorInFunction
                    << "Problem: cannot find myself in"
                    << " parents' children" << abort(FatalError);
            }
            else
            {
                FixedList<label, 8>& subCellIDs = subCellIDsPtr();

                subCellIDs[myPos] = -1;
            }
        }
    }

    // Mark splitCell as free
    split.parent_ = -2;
    split.addedCellIDsPtr_.clear();

    split.cellAvg_ = 0.0;

    // Add to cache of free splitCells
    freeSplitCells_.append(index);
}


// Mark entry in splitCells. Recursively mark its parent and subs.
void Foam::refinementMRAHistory::markSplit
(
    const label index,
    labelList& oldToNew,
    DynamicList<splitCell8>& newSplitCells
) const
{
    if (oldToNew[index] == -1)
    {
        // Not yet compacted.

        const splitCell8& split = splitCells_[index];

        oldToNew[index] = newSplitCells.size();
        newSplitCells.append(split);

        if (split.parent_ >= 0)
        {
            markSplit(split.parent_, oldToNew, newSplitCells);
        }
        if (split.addedCellIDsPtr_.valid())
        {
//            FixedList<label, 8> splits = getAddedCellsIndices(index);

//            forAll(splits, i)
            for (label i = 0; i < 8; i++)
            {
//                if (splits[i] >= 0)
                if (visibleCells_[split.addedCellIDsPtr_()[i]])
                {
                    markSplit(visibleCells_[split.addedCellIDsPtr_()[i]], oldToNew, newSplitCells);
                }
            }
        }
    }
}


// Mark index and all its descendants
void Foam::refinementMRAHistory::mark
(
    const label val,
    const label index,
    labelList& splitToVal
) const
{
    splitToVal[index] = val;

    const splitCell8& split = splitCells_[index];

    if (split.addedCellIDsPtr_.valid())
    {
//        FixedList<label, 8> splits = getAddedCellsIndices(index);

//        forAll(splits, i)
        for (label i = 0; i < 8; i++)
        {
//            if (splits[i] >= 0)
            if (visibleCells_[split.addedCellIDsPtr_()[i]])
            {
                mark(val, visibleCells_[split.addedCellIDsPtr_()[i]], splitToVal);
            }
        }
    }
}


Foam::label Foam::refinementMRAHistory::markCommonCells
(
    labelList& cellToCluster
) const
{
    label clusterI = 0;

    labelList splitToCluster(splitCells_.size(), -1);

    // Pass1: find top of all clusters
    forAll(visibleCells_, cellI)
    {
        label index = visibleCells_[cellI];

        if (index >= 0)
        {
            // Find highest ancestor
            while (splitCells_[index].parent_ != -1)
            {
                index = splitCells_[index].parent_;
            }

            // Mark tree with clusterI
            if (splitToCluster[index] == -1)
            {
                mark(clusterI, index, splitToCluster);
                clusterI++;
            }
        }
    }

    // Pass2: mark all cells with cluster
    cellToCluster.setSize(visibleCells_.size(), -1);

    forAll(visibleCells_, cellI)
    {
        label index = visibleCells_[cellI];

        if (index >= 0)
        {
            cellToCluster[cellI] = splitToCluster[index];
        }
    }

    return clusterI;
}


void Foam::refinementMRAHistory::add
(
    boolList& blockedFace,
    PtrList<labelList>& specifiedProcessorFaces,
    labelList& specifiedProcessor,
    List<labelPair>& explicitConnections
) const
{
    const polyMesh& mesh = dynamic_cast<const polyMesh&>(db());

    blockedFace.setSize(mesh.nFaces(), true);

    // Find common parent for all cells
    labelList cellToCluster;
    markCommonCells(cellToCluster);


    // Unblock all faces inbetween same cluster

    label nUnblocked = 0;

    forAll(mesh.faceNeighbour(), faceI)
    {
        label ownCluster = cellToCluster[mesh.faceOwner()[faceI]];
        label neiCluster = cellToCluster[mesh.faceNeighbour()[faceI]];

        if (ownCluster != -1 && ownCluster == neiCluster)
        {
            if (blockedFace[faceI])
            {
                blockedFace[faceI] = false;
                nUnblocked++;
            }
        }
    }

    if (refinementMRAHistory::debug)
    {
        reduce(nUnblocked, sumOp<label>());
        Info<< type() << " : unblocked " << nUnblocked << " faces" << endl;
    }

    syncTools::syncFaceList(mesh, blockedFace, andEqOp<bool>());
}


void Foam::refinementMRAHistory::apply
(
    const boolList& blockedFace,
    const PtrList<labelList>& specifiedProcessorFaces,
    const labelList& specifiedProcessor,
    const List<labelPair>& explicitConnections,
    labelList& decomposition
) const
{
    const polyMesh& mesh = dynamic_cast<const polyMesh&>(db());

    // Find common parent for all cells
    labelList cellToCluster;
    label nClusters = markCommonCells(cellToCluster);

    // Unblock all faces inbetween same cluster


    labelList clusterToProc(nClusters, -1);

    label nChanged = 0;

    forAll(mesh.faceNeighbour(), faceI)
    {
        label own = mesh.faceOwner()[faceI];
        label nei = mesh.faceNeighbour()[faceI];

        label ownCluster = cellToCluster[own];
        label neiCluster = cellToCluster[nei];

        if (ownCluster != -1 && ownCluster == neiCluster)
        {
            if (clusterToProc[ownCluster] == -1)
            {
                clusterToProc[ownCluster] = decomposition[own];
            }

            if (decomposition[own] != clusterToProc[ownCluster])
            {
                decomposition[own] = clusterToProc[ownCluster];
                nChanged++;
            }
            if (decomposition[nei] != clusterToProc[ownCluster])
            {
                decomposition[nei] = clusterToProc[ownCluster];
                nChanged++;
            }
        }
    }

    if (refinementMRAHistory::debug)
    {
        reduce(nChanged, sumOp<label>());
        Info<< type() << " : changed decomposition on " << nChanged
            << " cells" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementMRAHistory::refinementMRAHistory(const IOobject& io)
:
    regIOobject(io),
    refCount(),
    active_(false)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningInFunction
            << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
            << " does not support automatic rereading."
            << endl;
    }

    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }

    // When running in redistributPar + READ_IF_PRESENT it can happen
    // that some processors do have refinementMRAHistory and some don't so
    // test for active has to be outside of above condition.
    active_ = (returnReduce(visibleCells_.size(), sumOp<label>()) > 0);

    if (debug)
    {
        Pout<< "refinementMRAHistory::refinementMRAHistory :"
            << " constructed history from IOobject :"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << endl;
    }
}


//- Read or construct
Foam::refinementMRAHistory::refinementMRAHistory
(
    const IOobject& io,
    const List<splitCell8>& splitCells,
    const labelList& visibleCells,
    const bool active
)
:
    regIOobject(io),
    refCount(),
    active_(active),
    splitCells_(splitCells),
    freeSplitCells_(0),
    visibleCells_(visibleCells)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningInFunction
            << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
            << " does not support automatic rereading."
            << endl;
    }

    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }

    // Check indices.
    checkIndices();

    if (debug)
    {
        Pout<< "refinementMRAHistory::refinementMRAHistory :"
            << " constructed history from IOobject or components :"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << endl;
    }
}


// Construct from initial number of cells (all visible)
Foam::refinementMRAHistory::refinementMRAHistory
(
    const IOobject& io,
    const label nCells
)
:
    regIOobject(io),
    refCount(),
    active_(false),
    freeSplitCells_(0)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningInFunction
            << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
            << " does not support automatic rereading."
            << endl;
    }

    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }
    else
    {
        visibleCells_.setSize(nCells);
        splitCells_.setCapacity(nCells);

        for (label cellI = 0; cellI < nCells; cellI++)
        {
            visibleCells_[cellI] = cellI;
            splitCells_.append(splitCell8());
        }
    }

    active_ = (returnReduce(visibleCells_.size(), sumOp<label>()) > 0);


    // Check indices.
    checkIndices();

    if (debug)
    {
        Pout<< "refinementMRAHistory::refinementMRAHistory :"
            << " constructed history from IOobject or initial size :"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << endl;
    }
}


// Construct from initial number of cells (all visible)
Foam::refinementMRAHistory::refinementMRAHistory
(
    const IOobject& io,
    const label nCells,
    const bool active
)
:
    regIOobject(io),
    refCount(),
    active_(active),
    freeSplitCells_(0)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningInFunction
            << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
            << " does not support automatic rereading."
            << endl;
    }

    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }
    else
    {
        visibleCells_.setSize(nCells);
        splitCells_.setCapacity(nCells);

        for (label cellI = 0; cellI < nCells; cellI++)
        {
            visibleCells_[cellI] = cellI;
            splitCells_.append(splitCell8());
        }
    }

    // Check indices.
    checkIndices();

    if (debug)
    {
        Pout<< "refinementMRAHistory::refinementMRAHistory :"
            << " constructed history from IOobject or initial size :"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << endl;
    }
}


// Construct as copy
Foam::refinementMRAHistory::refinementMRAHistory
(
    const IOobject& io,
    const refinementMRAHistory& rh
)
:
    regIOobject(io),
    refCount(),
    active_(rh.active_),
    splitCells_(rh.splitCells()),
    freeSplitCells_(rh.freeSplitCells()),
    visibleCells_(rh.visibleCells())
{
    if (debug)
    {
        Pout<< "refinementMRAHistory::refinementMRAHistory : constructed initial"
            << " history." << endl;
    }
}


// Construct from multiple
Foam::refinementMRAHistory::refinementMRAHistory
(
    const IOobject& io,
    const UPtrList<const labelList>& cellMaps,
    const UPtrList<const refinementMRAHistory>& refs
)
:
    regIOobject(io),
    refCount(),
    active_(false)
{
    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        WarningInFunction
            << "read option IOobject::MUST_READ, READ_IF_PRESENT or "
            << "MUST_READ_IF_MODIFIED"
            << " suggests that a read constructor would be more appropriate."
            << endl;
    }

    const polyMesh& mesh = dynamic_cast<const polyMesh&>(db());


    // Determine offsets into splitCells
    labelList offsets(refs.size()+1);
    offsets[0] = 0;
    forAll(refs, refI)
    {
        const DynamicList<splitCell8>& subSplits = refs[refI].splitCells();
        offsets[refI+1] = offsets[refI]+subSplits.size();
    }

    // Construct merged splitCells
    splitCells_.setSize(offsets.last());
    forAll(refs, refI)
    {
        const labelList& cellMap = cellMaps[refI];

        const DynamicList<splitCell8>& subSplits = refs[refI].splitCells();
        forAll(subSplits, i)
        {
            splitCell8& newSplit = splitCells_[offsets[refI]+i];

            // Copy
            newSplit = subSplits[i];

            // Offset indices
            if (newSplit.parent_ >= 0)
            {
                newSplit.parent_ += offsets[refI];
            }

            if (newSplit.addedCellIDsPtr_.valid())
            {
                FixedList<label, 8>& splitIDs = newSplit.addedCellIDsPtr_();

                forAll(splitIDs, i)
                {
                    if (splitIDs[i] >= 0)
                    {
                        splitIDs[i] = cellMap[splitIDs[i]];
                    }
                }
            }
        }
    }


    // Construct merged visibleCells
    visibleCells_.setSize(mesh.nCells(), -1);
    forAll(refs, refI)
    {
        const labelList& cellMap = cellMaps[refI];
        const labelList& subVis = refs[refI].visibleCells();

        forAll(subVis, i)
        {
            label& newVis = visibleCells_[cellMap[i]];

            newVis = subVis[i];
            if (newVis >= 0)
            {
                newVis += offsets[refI];
            }
        }
    }


    // Is active if any of the refinementHistories is active (assumes active
    // flag parallel synchronised)
    active_ = false;
    forAll(refs, refI)
    {
        if (refs[refI].active())
        {
            active_ = true;
            break;
        }
    }

    // Check indices.
    checkIndices();

    if (debug)
    {
        Pout<< "refinementMRAHistory::refinementMRAHistory :"
            << " constructed history from multiple refinementHistories :"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << endl;
    }
}


// Construct from Istream
Foam::refinementMRAHistory::refinementMRAHistory(const IOobject& io, Istream& is)
:
    regIOobject(io),
    refCount(),
    splitCells_(is),
    freeSplitCells_(0),
    visibleCells_(is)
{
    active_ = (returnReduce(visibleCells_.size(), sumOp<label>()) > 0);

    // Check indices.
    checkIndices();

    if (debug)
    {
        Pout<< "refinementMRAHistory::refinementMRAHistory :"
            << " constructed history from Istream"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::refinementMRAHistory> Foam::refinementMRAHistory::clone
(
    const IOobject& io,
    // Per visible cell the processor it is going to
    const labelList& decomposition,
    // Per splitCell entry the processor it moves to
    const labelList& splitCellProc,
    // Per splitCell entry the number of live cells that move to that processor
    const labelList& splitCellNum,

    const label procI,

    // From old to new cells
    const labelList& cellMap,

    // From old to new splitCells
    labelList& oldToNewSplit
) const
{
    oldToNewSplit.setSize(splitCells_.size());
    oldToNewSplit = -1;

    // Compacted splitCells
    DynamicList<splitCell8> newSplitCells(splitCells_.size());

    // Loop over all entries. Note: could recurse like countProc so only
    // visit used entries but is probably not worth it.

    forAll(splitCells_, index)
    {
        if (splitCellProc[index] == procI && splitCellNum[index] == 8)
        {
            // Entry moves in its whole to procI
            oldToNewSplit[index] = newSplitCells.size();
            newSplitCells.append(splitCells_[index]);
        }
    }

    // Add live cells that are subsetted.
    forAll(visibleCells_, cellI)
    {
        label index = visibleCells_[cellI];

        if (index >= 0 && decomposition[cellI] == procI)
        {
            label parent = splitCells_[index].parent_;

            // Create new splitCell with parent
            oldToNewSplit[index] = newSplitCells.size();
            newSplitCells.append(splitCell8(parent));
        }
    }

    //forAll(oldToNewSplit, index)
    //{
    //    Pout<< "old:" << index << " new:" << oldToNewSplit[index]
    //        << endl;
    //}

    newSplitCells.shrink();

    // Renumber contents of newSplitCells
    forAll(newSplitCells, index)
    {
        splitCell8& split = newSplitCells[index];

        if (split.parent_ >= 0)
        {
            split.parent_ = oldToNewSplit[split.parent_];
        }
        if (split.addedCellIDsPtr_.valid())
        {
            FixedList<label, 8>& splitIDs = split.addedCellIDsPtr_();

            forAll(splitIDs, i)
            {
                if (splitIDs[i] >= 0)
                {
                    splitIDs[i] = cellMap[splitIDs[i]];
                }
            }
        }
    }


    // Count number of cells
    label nSub = 0;
    forAll(decomposition, cellI)
    {
        if (decomposition[cellI] == procI)
        {
            nSub++;
        }
    }

    labelList newVisibleCells(nSub);
    nSub = 0;

    forAll(visibleCells_, cellI)
    {
        if (decomposition[cellI] == procI)
        {
            label index = visibleCells_[cellI];
            if (index >= 0)
            {
                index = oldToNewSplit[index];
            }
            newVisibleCells[nSub++] = index;
        }
    }

    return autoPtr<refinementMRAHistory>
    (
        new refinementMRAHistory
        (
            io,
            newSplitCells,
            newVisibleCells,
            active_
        )
    );
}


Foam::autoPtr<Foam::refinementMRAHistory> Foam::refinementMRAHistory::clone
(
    const IOobject& io,
    const labelList& cellMap
) const
{
    if (active_)
    {
        // Mark selected cells with '1'
        labelList decomposition(visibleCells_.size(), 0);
        forAll(cellMap, i)
        {
            decomposition[cellMap[i]] = 1;
        }


        // Per splitCell entry the processor it moves to
        labelList splitCellProc(splitCells_.size(), -1);
        // Per splitCell entry the number of live cells that move to that
        // processor
        labelList splitCellNum(splitCells_.size(), 0);

        forAll(visibleCells_, cellI)
        {
            label index = visibleCells_[cellI];

            if (index >= 0)
            {
                countProc
                (
                    splitCells_[index].parent_,
                    decomposition[cellI],
                    splitCellProc,
                    splitCellNum
                );
            }
        }

        labelList oldToNewSplit;
        return clone
        (
            io,
            decomposition,
            splitCellProc,
            splitCellNum,
            1,      //procI,
            cellMap,
            oldToNewSplit
        );
    }
    else
    {
        return autoPtr<refinementMRAHistory>
        (
            new refinementMRAHistory
            (
                io,
                DynamicList<splitCell8>(0),
                labelList(0),
                false
            )
        );
    }
}


void Foam::refinementMRAHistory::resize(const label size)
{
    label oldSize = visibleCells_.size();

    if (debug)
    {
        Pout<< "refinementMRAHistory::resize from " << oldSize << " to " << size
            << " cells" << endl;
    }

    visibleCells_.setSize(size);

    // Set additional elements to -1.
    for (label i = oldSize; i < visibleCells_.size(); i++)
    {
        visibleCells_[i] = -1;
    }
}


void Foam::refinementMRAHistory::updateMesh(const mapPolyMesh& map)
{
    if (active())
    {
        const labelList& reverseCellMap = map.reverseCellMap();

        // Note that only the live cells need to be renumbered.

        labelList newVisibleCells(map.cellMap().size(), -1);

        forAll(visibleCells_, cellI)
        {
            if (visibleCells_[cellI] != -1)
            {
                label index = visibleCells_[cellI];

                // Check not already set
                if (splitCells_[index].addedCellIDsPtr_.valid())
                {
                    FatalErrorInFunction
                        << "Problem" << abort(FatalError);
                }

                label newCellI = reverseCellMap[cellI];

                if (newCellI >= 0)
                {
                    newVisibleCells[newCellI] = index;
                }
            }
        }

        forAll(splitCells_, sCI)
        {
            if (!splitCells_[sCI].addedCellIDsPtr_.empty())
            {
                FixedList<label, 8>& splitIDs = splitCells_[sCI].addedCellIDsPtr_();

                forAll(splitIDs, j)
                {
                    splitIDs[j] = reverseCellMap[splitIDs[j]];
                }
            }
        }

        if (debug)
        {
            Pout<< "refinementMRAHistory::updateMesh : from "
                << visibleCells_.size()
                << " to " << newVisibleCells.size()
                << " cells" << endl;
        }

        visibleCells_.transfer(newVisibleCells);
    }
}


// Update numbering for subsetting
void Foam::refinementMRAHistory::subset
(
    const labelList& pointMap,
    const labelList& faceMap,
    const labelList& cellMap
)
{
    if (active())
    {
        labelList newVisibleCells(cellMap.size(), -1);

        forAll(newVisibleCells, cellI)
        {
            label oldCellI = cellMap[cellI];

            label index = visibleCells_[oldCellI];

            // Check that cell is live (so its parent has no refinement)
            if (index >= 0 && splitCells_[index].addedCellIDsPtr_.valid())
            {
                FatalErrorInFunction
                    << "Problem" << abort(FatalError);
            }

            newVisibleCells[cellI] = index;
        }

        labelList reverseCellMap(visibleCells_.size(), -1);

        forAll(cellMap, cellI)
        {
            label oldCellI = cellMap[cellI];

            reverseCellMap[oldCellI] = cellI;
        }

        forAll(splitCells_, sCI)
        {
            if (!splitCells_[sCI].addedCellIDsPtr_.empty())
            {
                FixedList<label, 8>& splitIDs = splitCells_[sCI].addedCellIDsPtr_();

                forAll(splitIDs, j)
                {
                    splitIDs[j] = reverseCellMap[splitIDs[j]];
                }
            }
        }

        if (debug)
        {
            Pout<< "refinementMRAHistory::updateMesh : from "
                << visibleCells_.size()
                << " to " << newVisibleCells.size()
                << " cells" << endl;
        }

        visibleCells_.transfer(newVisibleCells);
    }
}


void Foam::refinementMRAHistory::countProc
(
    const label index,
    const label newProcNo,
    labelList& splitCellProc,
    labelList& splitCellNum
) const
{
    if (splitCellProc[index] != newProcNo)
    {
        // Different destination processor from other cells using this
        // parent. Reset count.
        splitCellProc[index] = newProcNo;
        splitCellNum[index] = 1;
    }
    else
    {
        splitCellNum[index]++;

        // Increment parent if whole splitCell moves to same processor
        if (splitCellNum[index] == 8)
        {
            if (debug)
            {
                Pout<< "Moving " << splitCellNum[index]
                    << " cells originating from cell " << index
                    << " from processor " << Pstream::myProcNo()
                    << " to processor " << splitCellProc[index]
                    << endl;
            }

            label parent = splitCells_[index].parent_;

            if (parent >= 0)
            {
                countProc(parent, newProcNo, splitCellProc, splitCellNum);
            }
        }
    }
}


void Foam::refinementMRAHistory::distribute(const mapDistributePolyMesh& map)
{
    if (!active())
    {
        FatalErrorInFunction
            << "Calling distribute on inactive history" << abort(FatalError);
    }


    if (!Pstream::parRun())
    {
        return;
    }

    // Remove unreferenced history.
    compact();

    //Pout<< nl << "--BEFORE:" << endl;
    //writeDebug();
    //Pout<< "---------" << nl << endl;


    // Distribution is only partially functional.
    // If all 8 cells resulting from a single parent are sent across in one
    // go it will also send across that part of the refinement history.
    // If however e.g. first 1 and then the other 7 are sent across the
    // history will not be reconstructed.

    // Determine clusters. This is per every entry in splitCells_ (that is
    // a parent of some refinement) a label giving the processor it goes to
    // if all its children are going to the same processor.

    // Per visible cell the processor it goes to.
    labelList destination(visibleCells_.size());

    const labelListList& subCellMap = map.cellMap().subMap();

    forAll(subCellMap, procI)
    {
        const labelList& newToOld = subCellMap[procI];

        forAll(newToOld, i)
        {
            label oldCellI = newToOld[i];

            destination[oldCellI] = procI;
        }
    }

    // Per splitCell entry the processor it moves to
    labelList splitCellProc(splitCells_.size(), -1);
    // Per splitCell entry the number of live cells that move to that processor
    labelList splitCellNum(splitCells_.size(), 0);

    forAll(visibleCells_, cellI)
    {
        label index = visibleCells_[cellI];

        if (index >= 0)
        {
            countProc
            (
                splitCells_[index].parent_,
                destination[cellI],
                splitCellProc,
                splitCellNum
            );
        }
    }

    //Pout<< "refinementMRAHistory::distribute :"
    //    << " splitCellProc:" << splitCellProc << endl;
    //
    //Pout<< "refinementMRAHistory::distribute :"
    //    << " splitCellNum:" << splitCellNum << endl;


    // Create subsetted refinement tree consisting of all parents that
    // move in their whole to other processor.
    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        //Pout<< "-- Subetting for processor " << procI << endl;

        // subsetted cellIDs of old cells
        labelList oldToNewCells(visibleCells_.size(), -1);

        const labelList& newToOld = subCellMap[procI];

        forAll(newToOld, i)
        {
            label oldCellI = newToOld[i];

            oldToNewCells[oldCellI] = i;
        }

        // From uncompacted to compacted splitCells.
        labelList oldToNew(splitCells_.size(), -1);

        // Compacted splitCells. Similar to subset routine below.
        DynamicList<splitCell8> newSplitCells(splitCells_.size());

        // Loop over all entries. Note: could recurse like countProc so only
        // visit used entries but is probably not worth it.

        forAll(splitCells_, index)
        {
            if (splitCellProc[index] == procI && splitCellNum[index] == 8)
            {
                // Entry moves in its whole to procI
                oldToNew[index] = newSplitCells.size();
                newSplitCells.append(splitCells_[index]);
            }
        }

        const labelList& subMap = subCellMap[procI];

        // Add live cells that are subsetted.
        forAll(visibleCells_, cellI)
        {
            label index = visibleCells_[cellI];

            if (index >= 0 && destination[cellI] == procI)
            {
                label parent = splitCells_[index].parent_;

                // Create new splitCell with parent
                oldToNew[index] = newSplitCells.size();
                newSplitCells.append(splitCell8(parent));
            }
        }

        //forAll(oldToNew, index)
        //{
        //    Pout<< "old:" << index << " new:" << oldToNew[index]
        //        << endl;
        //}

        newSplitCells.shrink();

        // Renumber contents of newSplitCells
        forAll(newSplitCells, index)
        {
            splitCell8& split = newSplitCells[index];

            if (split.parent_ >= 0)
            {
                split.parent_ = oldToNew[split.parent_];
            }
            if (split.addedCellIDsPtr_.valid())
            {
                FixedList<label, 8>& splitIDs = split.addedCellIDsPtr_();

                forAll(splitIDs, i)
                {
                    if (splitIDs[i] >= 0)
                    {
                        splitIDs[i] = oldToNewCells[splitIDs[i]];
                    }
                }
            }
        }

        // New visible cells.
        labelList newVisibleCells(subMap.size(), -1);

        forAll(subMap, newCellI)
        {
            label oldCellI = subMap[newCellI];

            label oldIndex = visibleCells_[oldCellI];

            if (oldIndex >= 0)
            {
                newVisibleCells[newCellI] = oldToNew[oldIndex];
            }
        }

        //Pout<< nl << "--Subset for domain:" << procI << endl;
        //writeDebug(newVisibleCells, newSplitCells);
        //Pout<< "---------" << nl << endl;


        // Send to neighbours
        OPstream toNbr(Pstream::blocking, procI);
        toNbr << newSplitCells << newVisibleCells;
    }


    // Receive from neighbours and merge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Remove all entries. Leave storage intact.
    splitCells_.clear();

    const polyMesh& mesh = dynamic_cast<const polyMesh&>(db());

    visibleCells_.setSize(mesh.nCells());
    visibleCells_ = -1;

    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        IPstream fromNbr(Pstream::blocking, procI);
        List<splitCell8> newSplitCells(fromNbr);
        labelList newVisibleCells(fromNbr);

        //Pout<< nl << "--Received from domain:" << procI << endl;
        //writeDebug(newVisibleCells, newSplitCells);
        //Pout<< "---------" << nl << endl;


        // newSplitCells contain indices only into newSplitCells so
        // renumbering can be done here.
        label offset = splitCells_.size();

        //Pout<< "**Renumbering data from proc " << procI << " with offset "
        //    << offset << endl;

        const labelList& constructMap = map.cellMap().constructMap()[procI];

        forAll(newSplitCells, index)
        {
            splitCell8& split = newSplitCells[index];

            if (split.parent_ >= 0)
            {
                split.parent_ += offset;
            }
            if (split.addedCellIDsPtr_.valid())
            {
                FixedList<label, 8>& splitIDs = split.addedCellIDsPtr_();

                forAll(splitIDs, i)
                {
                    if (splitIDs[i] >= 0)
                    {
                        splitIDs[i] = constructMap[splitIDs[i]];
                    }
                }
            }

            splitCells_.append(split);
        }


        // Combine visibleCell.

        forAll(newVisibleCells, i)
        {
            if (newVisibleCells[i] >= 0)
            {
                visibleCells_[constructMap[i]] = newVisibleCells[i] + offset;
            }
        }
    }
    splitCells_.shrink();

    //Pout<< nl << "--AFTER:" << endl;
    //writeDebug();
    //Pout<< "---------" << nl << endl;
}


void Foam::refinementMRAHistory::compact()
{
    if (debug)
    {
        Pout<< "refinementMRAHistory::compact() Entering with:"
            << " freeSplitCells_:" << freeSplitCells_.size()
            << " splitCells_:" << splitCells_.size()
            << " visibleCells_:" << visibleCells_.size()
            << endl;

        // Check all free splitCells are marked as such
        forAll(freeSplitCells_, i)
        {
            label index = freeSplitCells_[i];

            if (splitCells_[index].parent_ != -2)
            {
                FatalErrorInFunction
                    << "Problem index:" << index
                    << abort(FatalError);
            }
        }

        // Check none of the visible cells are marked as free
        forAll(visibleCells_, cellI)
        {
            if
            (
                visibleCells_[cellI] >= 0
             && splitCells_[visibleCells_[cellI]].parent_ == -2
            )
            {
                FatalErrorInFunction
                    << "Problem : visible cell:" << cellI
                    << " is marked as being free." << abort(FatalError);
            }
        }
    }

    DynamicList<splitCell8> newSplitCells(splitCells_.size());

    // From uncompacted to compacted splitCells.
    labelList oldToNew(splitCells_.size(), -1);

    // Mark all used splitCell entries. These are either indexed by visibleCells
    // or indexed from other splitCell entries.

    // Mark from visibleCells
    forAll(visibleCells_, cellI)
    {
        label index = visibleCells_[cellI];

        if (index >= 0)
        {
            // Make sure we only mark visible indices if they either have a
            // parent or subsplits.
            if
            (
                splitCells_[index].parent_ != -1
             || splitCells_[index].addedCellIDsPtr_.valid()
            )
            {
                markSplit(index, oldToNew, newSplitCells);
            }
        }
    }

    // Mark from splitCells
    forAll(splitCells_, index)
    {
        if (splitCells_[index].parent_ == -2)
        {
            // freed cell.
        }
        else if
        (
            splitCells_[index].parent_ == -1
         && splitCells_[index].addedCellIDsPtr_.empty()
        )
        {
            // recombined cell. No need to keep since no parent and no subsplits
            // Note that gets marked if reachable from other index!
        }
        else
        {
            // Is used element.
            markSplit(index, oldToNew, newSplitCells);
        }
    }


    // Now oldToNew is fully complete and compacted elements are in
    // newSplitCells.
    // Renumber contents of newSplitCells and visibleCells.
    forAll(newSplitCells, index)
    {
        splitCell8& split = newSplitCells[index];

        if (split.parent_ >= 0)
        {
            split.parent_ = oldToNew[split.parent_];
        }
    }


    if (debug)
    {
        Pout<< "refinementMRAHistory::compact : compacted splitCells from "
            << splitCells_.size() << " to " << newSplitCells.size() << endl;
    }

    splitCells_.transfer(newSplitCells);
    freeSplitCells_.clearStorage();


    if (debug)
    {
        Pout<< "refinementMRAHistory::compact() NOW:"
            << " freeSplitCells_:" << freeSplitCells_.size()
            << " splitCells_:" << splitCells_.size()
            << " newSplitCells:" << newSplitCells.size()
            << " visibleCells_:" << visibleCells_.size()
            << endl;
    }


    // Adapt indices in visibleCells_
    forAll(visibleCells_, cellI)
    {
        label index = visibleCells_[cellI];

        if (index >= 0)
        {
            // Note that oldToNew can be -1 so it resets newVisibleCells.
            visibleCells_[cellI] = oldToNew[index];
        }
        else
        {
            // Keep -1 value.
        }
    }
}


void Foam::refinementMRAHistory::writeDebug() const
{
    writeDebug(visibleCells_, splitCells_);
}


const Foam::FixedList<Foam::label, 8> Foam::refinementMRAHistory::getAddedCellsIndices
(
    const label parent
) const
{
    FixedList<Foam::label, 8> addedCells(-1);

    if (parent >= 0)
    {
        const splitCell8& parentSplit = splitCells_[parent];

        if (!parentSplit.addedCellIDsPtr_.empty())
        {
            // Allocate storage on parent for the 8 subcells.
            const FixedList<label, 8> &childrenIDs = parentSplit.addedCellIDsPtr_();

            forAll(addedCells, i)
            {
                label childI = childrenIDs[i];
                if (childI >= 0)
                {
                    addedCells[i] = visibleCells_[childI];
                }
            }
        }
        else
        {
            FatalErrorInFunction
                << "Cell with splitCellID " << parent
                << " has not been split."
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Invalid parentIndex " << parent
            << abort(FatalError);
    }

    return addedCells;
}


void Foam::refinementMRAHistory::storeSplit
(
    const label cellI,
    const labelList& addedCells
)
{
    label parentIndex = -1;

    if (visibleCells_[cellI] != -1)
    {
        // Was already live. The current live cell becomes the
        // parent of the cells split off from it.

        parentIndex = visibleCells_[cellI];

        // It is no longer live (note that actually cellI gets alive
        // again below since is addedCells[0])
        visibleCells_[cellI] = -1;
    }
    else
    {
        // Create 0th level. -1 parent to denote this.
        parentIndex = allocateSplitCell(-1, -1, -1);
    }

    // Create live entries for added cells that point to the
    // cell they were created from (parentIndex)
    forAll(addedCells, i)
    {
        label addedCellI = addedCells[i];

        // Create entries for the split off cells. All of them
        // are visible.
        visibleCells_[addedCellI] = allocateSplitCell(parentIndex, addedCellI, i);
    }
}


void Foam::refinementMRAHistory::combineCells
(
    const label masterCellI,
    const labelList& combinedCells
)
{
    // Save the parent structure
    label parentIndex = splitCells_[visibleCells_[masterCellI]].parent_;

    double avgVal = 0.0;

    // Remove the information for the combined cells
    forAll(combinedCells, i)
    {
        label cellI = combinedCells[i];

        avgVal += splitCells_[visibleCells_[cellI]].cellAvg_;

        freeSplitCell(visibleCells_[cellI]);
        visibleCells_[cellI] = -1;
    }

    splitCell8& parentSplit = splitCells_[parentIndex];
    parentSplit.cellAvg_ = avgVal * 0.125;
    parentSplit.addedCellIDsPtr_.reset(NULL);
    visibleCells_[masterCellI] = parentIndex;
}



//- Fill the cell averages of splitCells
void Foam::refinementMRAHistory::fillCellAvgs
(
    const volScalarField& vFld
)
{
    // TODO: add check for same size of vFld and visibleCells

    // update cellAvg for visible cells

    DynamicList<label> sameLevelSplits;
    PackedBoolList isFilled(splitCells_.size(), 0);

    forAll (visibleCells_, cellI)
    {
        label splitID = visibleCells_[cellI];

        if (splitID >= 0)
        {
            splitCell8& split = splitCells_[splitID];
            split.cellAvg_ = vFld[cellI];

            if (split.parent_ >= 0)
            {
                sameLevelSplits.append(splitID);
            }

            isFilled.set(splitID, 1);
        }
    }

    bool isChanging = true;

    while (isChanging)
    {
        isChanging = false;

        DynamicList<label> sameLevelSplitsNew;

        forAll (sameLevelSplits, i)
        {
            label splitID = sameLevelSplits[i];
            label parentSplitID = splitCells_[splitID].parent_;

            splitCell8& parentSplit = splitCells_[parentSplitID];

            if (!isFilled[parentSplitID])
            {
                FixedList<label, 8> splits = getAddedCellsIndices(parentSplitID);

                double cellAvg = 0.0;

                forAll (splits, j)
                {
                    cellAvg += splitCells_[splits[j]].cellAvg_;
                }

                parentSplit.cellAvg_ = cellAvg * 0.125;

                isFilled.set(parentSplitID, 1);

                if (parentSplit.parent_ >= 0)
                {
                    sameLevelSplitsNew.append(parentSplitID);
                    isChanging = true;
                }
            }
        }

        sameLevelSplits.transfer(sameLevelSplitsNew);
    }
}

bool Foam::refinementMRAHistory::read()
{
    bool ok = readData(readStream(typeName));
    close();

    active_ = (returnReduce(visibleCells_.size(), sumOp<label>()) > 0);

    return ok;
}


bool Foam::refinementMRAHistory::readData(Istream& is)
{
    is >> *this;
    return !is.bad();
}


bool Foam::refinementMRAHistory::writeData(Ostream& os) const
{
    os << *this;

    return os.good();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, refinementMRAHistory& rh)
{
    rh.freeSplitCells_.clearStorage();

    is >> rh.splitCells_ >> rh.visibleCells_;

    // Check indices.
    rh.checkIndices();

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const refinementMRAHistory& rh)
{
    const_cast<refinementMRAHistory&>(rh).compact();

    return os   << "// splitCells" << nl
                << rh.splitCells_ << nl
                << "// visibleCells" << nl
                << rh.visibleCells_;
}


// ************************************************************************* //
