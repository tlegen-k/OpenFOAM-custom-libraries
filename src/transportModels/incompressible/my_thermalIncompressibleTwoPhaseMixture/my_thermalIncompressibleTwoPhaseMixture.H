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
    Foam::my_thermalIncompressibleTwoPhaseMixture

Description
    A two-phase thermal incompressible transportModel

SourceFiles
    my_thermalIncompressibleTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef my_thermalIncompressibleTwoPhaseMixture_H
#define my_thermalIncompressibleTwoPhaseMixture_H

#include "incompressible/transportModel/transportModel.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "twoPhaseMixture.H"
#include "conductivityModel.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class my_thermalIncompressibleTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class my_thermalIncompressibleTwoPhaseMixture
:
    public IOdictionary,
    public transportModel,
    public twoPhaseMixture
{
protected:

    // Protected data

        autoPtr<viscosityModel> nuModel1_;
        autoPtr<viscosityModel> nuModel2_;

        dimensionedScalar rho1_;
        dimensionedScalar rho2_;

        //--modified_start--//
        dimensionedScalar cp1_;
        dimensionedScalar cp2_;
        dimensionedScalar Pr1_;
        dimensionedScalar Pr2_;
        //--modified_end--//

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        // Added by ASR
        autoPtr<conductivityModel> lambdaModel1_;
        autoPtr<conductivityModel> lambdaModel2_;

        // Added by ASR
        volScalarField lambda_;
        volScalarField cp_;
        volScalarField rho_;

        // Added by ASR
        // const fvMesh& mesh_;

        // //--modified_start--//
        // const volScalarField& T_;
        // //--modified_end--//

        volScalarField nu_;

        // //--modified_start--//
        	// volScalarField kappa_;
        // //--modified_end--//

    // Protected Member Functions

        //- Calculate and return the laminar viscosity
        void calcNu();

        //- ASR Calculate and return the laminar conductivity
        void calcLambda();

        //- ASR Calculate specific heat
        void calcCp();

        //- ASR Calculate local density
        void calcRho();

public:

    // TypeName("my_thermalIncompressibleTwoPhaseMixture");


    // Constructors

        //- Construct from components
        my_thermalIncompressibleTwoPhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
            // //--modified_start--//
            // const volScalarField& T_;
            // //--modified_end--//
        );


    //- Destructor
    virtual ~my_thermalIncompressibleTwoPhaseMixture()
    {}


    // Member Functions

        //- Return const-access to phase1 viscosityModel
        const viscosityModel& nuModel1() const
        {
            return *nuModel1_;
        }

        //- Return const-access to phase2 viscosityModel
        const viscosityModel& nuModel2() const
        {
            return *nuModel2_;
        }

        //- RAS Return const-access to phase1 conductivityModel
        const conductivityModel& lambdaModel1() const
        {
            return lambdaModel1_();
        }

        //- RAS Return const-access to phase2 conductivityModel
        const conductivityModel& lambdaModel2() const
        {
            return lambdaModel2_();
        }

        //- Return const-access to phase1 density
        const dimensionedScalar& rho1() const
        {
            return rho1_;
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& rho2() const
        {
            return rho2_;
        };

        //-----------modified_start------------//
        //- Return const-access to phase1 density
        const dimensionedScalar& cp1() const
        {
            return cp1_;
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& cp2() const
        {
            return cp2_;
        };

        //- Prandtl number and kappaf() are not used with phase change implementation

        //- Return const-access to phase1 density
        const dimensionedScalar& Pr1() const
        {
            return Pr1_;
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& Pr2() const
        {
            return Pr2_;
        };
        //-----------modified_end------------//

        //- Return const-access to the mixture velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const;

        //- RAS Return the conductivity
        tmp<volScalarField> lambda() const;

        //- Return the dynamic laminar viscosity on patch
        tmp<scalarField> mu(const label patchI) const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;

        //- RAS Return the face-interpolated dynamic laminar conductivity
        tmp<surfaceScalarField> lambdaf() const;

        //- Return the kinematic laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the kinematic laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Return the face-interpolated kinematic laminar viscosity
        tmp<surfaceScalarField> nuf() const;

        //-----------modified_start------------//
        	//- Return the face-interpolated conductivity
        tmp<surfaceScalarField> kappaf() const;
        //-----------modified_end------------//

        //- ASR Return the specific heats
        tmp<volScalarField> cp() const;

    	//- Return the face-interpolated thermal capacity
        tmp<surfaceScalarField> cpf() const;

        //- ASR for calculating thermal diffusivity
        tmp<volScalarField> alpha() const;

        //- ASR for calculating alpha at faces
        tmp<surfaceScalarField> alphaf() const;

        //- ASR return the density
        tmp<volScalarField> rho() const;

        //- Correct the laminar viscosity, conductivity, thermal capacity and density
        virtual void correct()
        {
            calcNu();
            calcLambda();
            calcCp();
            calcRho();
        }

        //- Read base transportProperties dictionary
        virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
