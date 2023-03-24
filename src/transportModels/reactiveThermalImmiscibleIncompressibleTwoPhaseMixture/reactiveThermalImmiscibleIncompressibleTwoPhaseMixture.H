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
    Foam::reactiveThermalImmiscibleIncompressibleTwoPhaseMixture

Description
    A reactive thermal immiscible incompressible two-phase mixture transport model

SourceFiles
    reactiveThermalImmiscibleIncompressibleTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef reactiveThermalImmiscibleIncompressibleTwoPhaseMixture_H
#define reactiveThermalImmiscibleIncompressibleTwoPhaseMixture_H

#include "my_thermalIncompressibleTwoPhaseMixture.H"
#include "interfaceProperties.H"

#include "speciesTable.H"
#include "PtrList.H"
#include "volFields.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class reactiveThermalImmiscibleIncompressibleTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class reactiveThermalImmiscibleIncompressibleTwoPhaseMixture
:
    public my_thermalIncompressibleTwoPhaseMixture,
    public interfaceProperties //,
{
protected:
    // Protected data

    // modified_start from basicMultiComponentMixture //
        //- Table of specie names
        speciesTable species_;

        //- List of specie active flags
        List<bool> active_;

        //- Species mass fractions
        PtrList<volScalarField> Y_;

        //- Number of species
        label nSpecie_;

        //- Number of reactions
        label nReaction_;

        //- List of reaction rate per species [kg/m3/s]
        PtrList<volScalarField::Internal> RR_;

public:
    //- Run time type information
    TypeName("reactiveThermalImmiscibleIncompressibleTwoPhaseMixture");

    // Constructors

        //- Construct from components
        reactiveThermalImmiscibleIncompressibleTwoPhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            // //--modified_start from basicMultiComponentMixture--//
            const dictionary&,
            const wordList& specieNames,
            const fvMesh&,
            const word&
            // //--modified_end--//
        );


    //- Destructor
    virtual ~reactiveThermalImmiscibleIncompressibleTwoPhaseMixture()
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

        // modified_start from basicMultiComponentMixture //
    // Member functions

        //- Return the table of species
        inline const speciesTable& species() const;

        //- Does the mixture include this specie?
        inline bool contains(const word& specieName) const;

        //- Return true for active species
        inline bool active(label speciei) const;

        //- Return the bool list of active species
        inline const List<bool>& active() const;

        //- Set speciei active
        inline void setActive(label speciei);

        //- Set speciei inactive
        inline void setInactive(label speciei);

        //- Return the mass-fraction fields
        inline PtrList<volScalarField>& Y();

        //- Return the const mass-fraction fields
        inline const PtrList<volScalarField>& Y() const;

        //- Return the mass-fraction field for a specie given by index
        inline volScalarField& Y(const label i);

        //- Return the const mass-fraction field for a specie given by index
        inline const volScalarField& Y(const label i) const;

        //- Return the mass-fraction field for a specie given by name
        inline volScalarField& Y(const word& specieName);

        //- Return the const mass-fraction field for a specie given by name
        inline const volScalarField& Y(const word& specieName) const;
        // modified_end from basicMultiComponentMixture

        // modified_start from StandardChemistryModel //

        //- The number of species
        virtual inline label nSpecie() const;

        //- The number of reactios
        virtual inline label nReaction() const;


        // ODE functions (overriding abstract functions in ODE.H)

        //-Number of ODEs to solve
        inline virtual label nEqns() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reactiveThermalImmiscibleIncompressibleTwoPhaseMixtureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
