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

Class
    Foam::RASModels::Chien

Group
    grpRASTurbulence

Description
    Chien's k-epsilon turbulence model for incompressible and compressible
    flows.

    Reference:
    \verbatim
        Chien, K.-Y. (1982)
        “Predictions of Channel and Boundary-Layer Flows with a Low-
        Reynolds-Number Turbulence Model,” 
        AIAA Journal, vol. 20, pp. 33–38, 1982
    \endverbatim

    The default model coefficients are
    \verbatim
        ChienCoeffs
        {
            Cmu         0.09;
            C1          1.35;
            C2          1.8;
            C3          0.0;
            sigmak      1.0;
            sigmaEps    1.3;
        }
    \endverbatim

SourceFiles
    Chien.C

\*---------------------------------------------------------------------------*/

#ifndef Chien_H
#define Chien_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class Chien Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class Chien
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        Chien(const Chien&) = delete;

        //- No copy assignment
        void operator=(const Chien&) = delete;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar C3_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;

        // Fields

            volScalarField k_;
            volScalarField epsilon_;


    // Protected Member Functions
            tmp<volScalarField> fMu() const;
            tmp<volScalarField> f2() const;
            virtual void correctNut();
            virtual tmp<fvScalarMatrix> kSource() const;
            virtual tmp<fvScalarMatrix> epsilonSource() const;

        public:
            typedef typename BasicTurbulenceModel::alphaField alphaField;
            typedef typename BasicTurbulenceModel::rhoField rhoField;
            typedef typename BasicTurbulenceModel::transportModel transportModel;

            //- Runtime type information
            TypeName("Chien");

            // Constructors

            //- Construct from components
            Chien(
                const alphaField &alpha,
                const rhoField &rho,
                const volVectorField &U,
                const surfaceScalarField &alphaRhoPhi,
                const surfaceScalarField &phi,
                const transportModel &transport,
                const word &propertiesName = turbulenceModel::propertiesName,
                const word &type = typeName);

            //- Destructor
            virtual ~Chien()
            {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();
    
        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    (this->nut_/sigmak_ + this->nu())
                )
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DepsilonEff",
                    (this->nut_/sigmaEps_ + this->nu())
                )
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Chien.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
