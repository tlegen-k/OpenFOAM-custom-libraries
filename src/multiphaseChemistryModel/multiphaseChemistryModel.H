/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 AUTHOR,AFFILIATION
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
    Foam::multiphaseChemistryModel

Description
    Library to solve chemical kinetics in multiphase flow

SourceFiles
    multiphaseChemistryModelI.H
    multiphaseChemistryModel.C
    multiphaseChemistryModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseChemistryModel_H
#define multiphaseChemistryModel_H

#include "reactionMultiphase.H" //defines class ReactionThermo
#include "ODESystem.H"
#include "volFields.H"
#include "simpleMatrix.H"
#include "fvScalarMatrix.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class multiphaseChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class multiphaseChemistryModel
// :
    // public ODESystem
{
    // Private Member Functions

        //- No copy construct
        multiphaseChemistryModel
        (
            const multiphaseChemistryModel&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const multiphaseChemistryModel&
        ) = delete;


protected:

    // typedef ThermoType thermoType;


    // Protected data
        //- from basicChemistryModel
        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the field of specie mass fractions
        // PtrList<volScalarField>& Y_;

        //- Reference to the species concentrations
        PtrList<volScalarField>& C_;

        //- Species molar weights
        List<scalar> W_;

        //- Reference to the temperature field
        volScalarField& T_;

        //- Reference to the VOF field
        volScalarField& alpha_;

        //- Reference to the density field
        volScalarField& rho_;

        //- Reactions
        // const PtrList<reactionMultiphase>& reactions_;
        const List<reactionMultiphase>& reactions_;

        //- Number of species
        label nSpecie_;

        //- Number of reactions
        label nReaction_;

        //- Temperature below which the reaction rates are assumed 0
        scalar Treact_;

        //- List of reaction rate per specie [kg/m3/s]
        PtrList<volScalarField::Internal> RR_;

        //- Temporary concentration field
        mutable scalarField c_;

        //- Temporary rate-of-change of concentration field
        mutable scalarField dcdt_;


    // Protected Member Functions

        //- Write access to chemical source terms
        //  (e.g. for multi-chemistry model)
        inline PtrList<volScalarField::Internal>& RR();


public:
    // Constructors

        //- Construct from reactionsList, species volScalarField, and VOF function
        //- Rationale: distinguish liquid and gas phase reactions
        multiphaseChemistryModel
        (
            List<reactionMultiphase>& reactions,
            // PtrList<volScalarField>& Y,
            PtrList<volScalarField>& C,
            List<scalar> W,
            volScalarField& T,
            volScalarField& alpha,
            volScalarField& rho,
            const fvMesh& mesh
        );
        
    //- Destructor
    virtual ~multiphaseChemistryModel();


    // Member Functions
        //- from basicChemistryModel
        //- Return const access to the mesh database
        inline const fvMesh& mesh() const;

        //- The reactions
        inline const List<reactionMultiphase>& reactions() const;

        //- The number of species
        virtual inline label nSpecie() const;

        //- The number of reactions
        virtual inline label nReaction() const;

        //- Temperature below which the reaction rates are assumed 0
        inline scalar Treact() const;

        //- Define smoothstep function that takes alpha
        inline scalar smoothstep(scalar edge0, scalar edge1, scalar x) const;
        
        //- Temperature below which the reaction rates are assumed 0
        inline scalar& Treact();
        
        //- Return non-const access to the local temperature field [kg/m^3]
        inline tmp<volScalarField> T() const;

        //- Temperature for patch
        inline tmp<scalarField> T(const label patchi) const;

        //- Return non-const access to the local VOF field [-]
        inline tmp<volScalarField> alpha() const;

        //- VOF for patch
        inline tmp<scalarField> alpha(const label patchi) const;

        //- Density [kg/m^3]
        inline tmp<volScalarField> rho() const;

        //- Density for patch [kg/m^3]
        inline tmp<scalarField> rho(const label patchi) const;

        //- Return non-const access to the local density field [kg/m^3]
        inline volScalarField& rho();

        //- Return const access to the species concentration fields
        inline PtrList<volScalarField>& C();

        //- dc/dt = omega, rate of change in concentration, for each species
        virtual void omega
        (
            const scalarField&  c,
            const scalar&       T,
            const scalar&       alpha,
            // const scalar     p,
            scalarField&        dcdt
        ) const;

        //- Return the reaction rate for reaction r and the reference
        //  species and characteristic times
        virtual scalar omega
        (
            const reactionMultiphase&   r,
            const scalarField&          c,
            const scalar&               T,
            const scalar&               alpha,
            // const scalar             p,
            scalar&                     pf,
            scalar&                     cf,
            label&                      lRef
        ) const;


        //- Return the reaction rate for iReaction and the reference
        //  species and characteristic times
        virtual scalar omegaI
        (
            label iReaction,
            const scalarField& c,
            const scalar& T,
            const scalar& alpha,
            // const scalar p,
            scalar& pf,
            scalar& cf,
            label& lRef
        ) const;

        //- Calculates the reaction rates
        virtual void calculate();


        // Chemistry model functions (overriding abstract functions in
        // basicChemistryModel.H)

            //- Return const access to the chemical source terms for specie, i
            inline const volScalarField::Internal& RR
            (
                const label i
            ) const;

            //- Return non const access to chemical source terms [kg/m3/s]
            virtual volScalarField::Internal& RR
            (
                const label i
            );

            //- Return reaction rate of the speciei in reactionI
            virtual tmp<volScalarField::Internal> calculateRR
            (
                const label reactionI,
                const label speciei
            ) const;

            //- Fuel consumption rate matrix, i.e. source term for fuel equation
            tmp<fvScalarMatrix> R(volScalarField& C, label specieI);


        // ODE functions (overriding abstract functions in ODE.H)

            //- Number of ODE's to solve
            inline virtual label nEqns() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multiphaseChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
