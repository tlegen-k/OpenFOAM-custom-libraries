/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 OpenFOAM Foundation
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
    Foam::my_thermalImmiscibleIncompressibleTwoPhaseMixture

Description
    A thermal immiscible incompressible two-phase mixture transport model

SourceFiles
    my_thermalImmiscibleIncompressibleTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef my_thermalImmiscibleIncompressibleTwoPhaseMixture_H
#define my_thermalImmiscibleIncompressibleTwoPhaseMixture_H

#include "my_thermalIncompressibleTwoPhaseMixture.H"
#include "interfaceProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class my_thermalImmiscibleIncompressibleTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class my_thermalImmiscibleIncompressibleTwoPhaseMixture
:
    public my_thermalIncompressibleTwoPhaseMixture,
    public interfaceProperties
{

public:

    // Constructors

        //- Construct from components
        my_thermalImmiscibleIncompressibleTwoPhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );

    //- Destructor
    virtual ~my_thermalImmiscibleIncompressibleTwoPhaseMixture()
    {}

    // Member Functions

        //- Correct the transport and interface properties
        virtual void correct()
        {
            my_thermalIncompressibleTwoPhaseMixture::correct();
            interfaceProperties::correct();
        }

        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
