/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "ensightSurfaceReader.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ensightSurfaceReader, 0);
    addToRunTimeSelectionTable(surfaceReader, ensightSurfaceReader, fileName);
}


void Foam::ensightSurfaceReader::skip(const label n, IFstream& is) const
{
    label i = 0;
    token t;
    while (is.good() && (i < n))
    {
        is >> t;
        i++;

        if (debug)
        {
            Info<< "Skipping token " << t << endl;
        }
    }

    if (i != n)
    {
        WarningInFunction
            << "Requested to skip " << n << "tokens, but stream exited after "
            << i << " tokens. Last token read: " << t
            << endl;
    }
}


void Foam::ensightSurfaceReader::readGeometryHeader(ensightReadFile& is) const
{
    // Binary flag string if applicable
    is.readBinaryHeader();

    string buffer;

    // Ensight Geometry File
    is.read(buffer);
    if (debug) Info<< "buffer: " << buffer << endl;

    // Description - 1
    is.read(buffer);
    if (debug) Info<< "buffer: " << buffer << endl;

    // Node info
    is.read(buffer);
    if (debug) Info<< "buffer: " << buffer << endl;

    // Element info
    is.read(buffer);
    if (debug) Info<< "buffer: " << buffer << endl;

    // Part
    is.read(buffer);
    if (debug) Info<< "buffer: " << buffer << endl;

    // Part number
    label ibuffer;
    is.read(ibuffer);
    if (debug) Info<< "ibuffer: " << ibuffer << endl;

    // Description - 2
    is.read(buffer);
    if (debug) Info<< "buffer: " << buffer << endl;

    // Co-ordinates
    is.read(buffer);
    if (debug) Info<< "buffer: " << buffer << endl;
}


void Foam::ensightSurfaceReader::debugSection
(
    const word& expected,
    IFstream& is
) const
{
    word actual(is);

    if (expected != actual)
    {
        FatalIOErrorInFunction(is)
            << "Expected section header '" << expected
            << "' but read the word '" << actual << "'"
            << exit(FatalIOError);
    }

    if (debug)
    {
        Info<< "Read section header: " << expected << endl;
    }
}


void Foam::ensightSurfaceReader::readCase(IFstream& is)
{
    if (debug)
    {
        InfoInFunction<< endl;
    }

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << is.name()
            << exit(FatalError);
    }

    string buffer;

    // Read the file
    debugSection("FORMAT", is);
    skip(3, is); // type: ensight gold

    debugSection("GEOMETRY", is);
    readSkip(is, 2, meshFileName_);

    debugSection("VARIABLE", is);

    DynamicList<word> fieldNames(10);
    DynamicList<word> fieldFileNames(10);
    word fieldName;
    word fieldFileName;
    while (is.good())
    {
        word primitiveType(is); // scalar, vector

        if (primitiveType == "TIME")
        {
            break;
        }

        readSkip(is, 3, fieldName); // p, U etc
        fieldNames.append(fieldName);

        is >> fieldFileName; // surfaceName.****.fieldName
        fieldFileNames.append(fieldFileName);
    }
    fieldNames_.transfer(fieldNames);
    fieldFileNames_.transfer(fieldFileNames);

    if (debug)
    {
        Info<< "fieldNames: " << fieldNames_ << nl
            << "fieldFileNames: " << fieldFileNames_ << endl;
    }

    // Start reading time information
    skip(3, is); // time set: 1
    readSkip(is, 3, nTimeSteps_);
    readSkip(is, 3, timeStartIndex_);
    readSkip(is, 2, timeIncrement_);

    if (debug)
    {
        Info<< "nTimeSteps: " << nTimeSteps_ << nl
            << "timeStartIndex: " << timeStartIndex_ << nl
            << "timeIncrement: " << timeIncrement_ << endl;
    }

    // Read the time values
    skip(2, is);
    timeValues_.setSize(nTimeSteps_);
    for (label i = 0; i < nTimeSteps_; i++)
    {
        scalar t(readScalar(is));

        timeValues_[i].value() = t;
        // TODO: use character representation of t directly instead of
        // regenerating from scalar value
        timeValues_[i].name() = Foam::name(t);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightSurfaceReader::ensightSurfaceReader(const fileName& fName)
:
    surfaceReader(fName),
    streamFormat_(IOstream::ASCII),
    baseDir_(fName.path()),
    meshFileName_(),
    fieldNames_(),
    fieldFileNames_(),
    nTimeSteps_(0),
    timeStartIndex_(0),
    timeIncrement_(1),
    timeValues_(),
    surfPtr_(NULL)
{
    IFstream is(fName);
    readCase(is);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightSurfaceReader::~ensightSurfaceReader()
{}


// * * * * * * * * * * * * * Public Member Functions   * * * * * * * * * * * //

const Foam::meshedSurface& Foam::ensightSurfaceReader::geometry()
{
    if (debug)
    {
        InfoInFunction<< endl;
    }

    if (!surfPtr_.valid())
    {
        IFstream isBinary(baseDir_/meshFileName_, IOstream::BINARY);

        if (!isBinary.good())
        {
            FatalErrorInFunction
                << "Cannot read file " << isBinary.name()
                << exit(FatalError);
        }

        streamFormat_ = IOstream::BINARY;
        {
            istream& is = isBinary.stdStream();

            char buffer[80];
            is.read(buffer, 80);

            char test[80];
            label nChar = 0;
            for (label i = 0; i < 80; ++i)
            {
                if (buffer[i] == '\0')
                {
                    break;
                }
                test[i] = buffer[i];
                nChar++;
            }

            string testStr(test, nChar);

            if
            (
                (testStr.find("binary", 0) == string::npos)
             && (testStr.find("Binary", 0) == string::npos)
            )
            {
                streamFormat_ = IOstream::ASCII;
            }
        }

        if (debug)
        {
            Info<< "stream format: ";
            if (streamFormat_ == IOstream::ASCII)
            {
                Info<< "ascii" << endl;
            }
            else
            {
                Info<< "binary" << endl;
            }
        }
            

        ensightReadFile is(baseDir_/meshFileName_, streamFormat_);

        if (debug)
        {
            Info<< "File: " << is.name() << endl;
        }

        readGeometryHeader(is);
        
        label nPoints;
        is.read(nPoints);

        if (debug)
        {
            Info<< "nPoints: " << nPoints << endl;
        }

        pointField points(nPoints);
        {
            scalarField x(nPoints);
            for (label dir = 0; dir < 3; dir++)
            {
                forAll(points, pointI)
                {
                    is.read(x[pointI]);
                }

                points.replace(dir, x);
            }
        }


        // Read faces - may be a mix of tris, quads and polys
        DynamicList<face> faces(ceil(nPoints/3));
        DynamicList<Tuple2<string, label> > schema(faces.size());
        string faceType = "";
        label nFace = 0;
        while (is.good()) // (is.peek() != EOF)
        {
            is.read(faceType);

            if (!is.good())
            {
                break;
            }

            if (debug)
            {
                Info<< "faceType: " << faceType << endl;
            }

            if (faceType == "tria3")
            {
                is.read(nFace);

                label np = 3;
                for (label faceI = 0; faceI < nFace; ++faceI)
                {
                    face f(np);
                    for (label fpI = 0; fpI < np; fpI++)
                    {
                        is.read(f[fpI]);
                    }

                    faces.append(f);
                }
            }
            else if (faceType == "quad4")
            {
                is.read(nFace);

                label np = 4;
                for (label faceI = 0; faceI < nFace; ++faceI)
                {
                    face f(np);
                    for (label fpI = 0; fpI < np; fpI++)
                    {
                        is.read(f[fpI]);
                    }

                    faces.append(f);
                }
            }
            else if (faceType == "nsided")
            {
                is.read(nFace);

                labelList np(nFace);
                for (label faceI = 0; faceI < nFace; ++faceI)
                {
                    is.read(np[faceI]);
                }
                for (label faceI = 0; faceI < nFace; ++faceI)
                {
                    face f(np[faceI]);
                    for (label fpI = 0; fpI < f.size(); ++fpI)
                    {
                        is.read(f[fpI]);
                    }

                    faces.append(f);
                }
            }
            else
            {
                if (debug)
                {
                    WarningInFunction
                        << "Unknown face type: " << faceType
                        << ".  Aborting read and continuing with current "
                        << "elements only" << endl;
                }

                break;
            }
            schema.append(Tuple2<string, label>(faceType, nFace));
        }

        schema_.transfer(schema);

        if (debug)
        {
            Info<< "read nFaces: " << faces.size() << nl
                << "file schema: " << schema_ << endl;
        }

        // Convert from 1-based Ensight addressing to 0-based OF addressing
        forAll(faces, faceI)
        {
            face& f = faces[faceI];

            forAll(f, fpI)
            {
                f[fpI]--;
            }
        }

        surfPtr_.reset(new meshedSurface(xferMove(points), faces.xfer()));
    }

    return surfPtr_();
}


Foam::instantList Foam::ensightSurfaceReader::times() const
{
    return timeValues_;
}


Foam::wordList Foam::ensightSurfaceReader::fieldNames
(
    const label timeIndex
) const
{
    return fieldNames_;
}


Foam::tmp<Foam::Field<Foam::scalar> > Foam::ensightSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const scalar& refValue
) const
{
    return readField<scalar>(timeIndex, fieldIndex);
}


Foam::tmp<Foam::Field<Foam::vector> > Foam::ensightSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const vector& refValue
) const
{
    return readField<vector>(timeIndex, fieldIndex);
}


Foam::tmp<Foam::Field<Foam::sphericalTensor> >
Foam::ensightSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const sphericalTensor& refValue
) const
{
    return readField<sphericalTensor>(timeIndex, fieldIndex);
}


Foam::tmp<Foam::Field<Foam::symmTensor> > Foam::ensightSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const symmTensor& refValue
) const
{
    return readField<symmTensor>(timeIndex, fieldIndex);
}


Foam::tmp<Foam::Field<Foam::tensor> > Foam::ensightSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const tensor& refValue
) const
{
    return readField<tensor>(timeIndex, fieldIndex);
}


// ************************************************************************* //
