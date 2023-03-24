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

\*---------------------------------------------------------------------------*/

#include "my_thermalIncompressibleTwoPhaseMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// namespace Foam
// {
//     defineTypeNameAndDebug(my_thermalIncompressibleTwoPhaseMixture, 0);
// }


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::my_thermalIncompressibleTwoPhaseMixture::calcNu()
{
    nuModel1_->correct();
    nuModel2_->correct();

    const volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu()/(limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_);
}

// Calculate and return the thermal conductivity
void Foam::my_thermalIncompressibleTwoPhaseMixture::calcLambda()
{
    // Apply thermal conductivity corrections
    lambdaModel1_->correct();
    lambdaModel2_->correct();

    // We may need to calculate lambda here, somehow
    lambda_ = lambda();
}

//- Calculate specific heat
void Foam::my_thermalIncompressibleTwoPhaseMixture::calcCp()
{
    cp_ = cp();
}

//- Calculate density
void Foam::my_thermalIncompressibleTwoPhaseMixture::calcRho()
{
    rho_ = rho();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::my_thermalIncompressibleTwoPhaseMixture::my_thermalIncompressibleTwoPhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    twoPhaseMixture(U.mesh(), *this),

    nuModel1_
    (
        viscosityModel::New
        (
            "nu1",
            subDict(phase1Name_),
            U,
            phi
        )
    ),
    nuModel2_
    (
        viscosityModel::New
        (
            "nu2",
            subDict(phase2Name_),
            U,
            phi
        )
    ),
    lambdaModel1_
    (
        conductivityModel::New
        (
            "lambda1",
            subDict(phase1Name_),
            U,
            phi
        )
    ),
    lambdaModel2_
    (
        conductivityModel::New
        (
            "lambda2",
            subDict(phase2Name_),
            U,
            phi
        )
    ),
    rho1_("rho", dimDensity, nuModel1_->viscosityProperties()),
    rho2_("rho", dimDensity, nuModel2_->viscosityProperties()),
    cp1_("cp", dimensionSet(0, 2, -2, -1, 0, 0, 0), nuModel1_->viscosityProperties().lookup("cp")),
    cp2_("cp", dimensionSet(0, 2, -2, -1, 0, 0, 0), nuModel2_->viscosityProperties().lookup("cp")),
    Pr1_("Pr", dimensionSet(0, 0, 0, 0, 0, 0, 0), nuModel1_->viscosityProperties().lookup("Pr")),
    Pr2_("Pr", dimensionSet(0, 0, 0, 0, 0, 0, 0), nuModel2_->viscosityProperties().lookup("Pr")),
    U_(U),
    phi_(phi),
//  T_(T),
    nu_
    (
        IOobject
        (
            "nu",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar(dimViscosity, Zero),
        calculatedFvPatchScalarField::typeName
    ),
    //for the thermal conductivity
    lambda_
    (
        IOobject
        (
            "lambda",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("lambda", dimensionSet(1, 1, -3, -1, 0), 0),
        calculatedFvPatchScalarField::typeName
    ),
    cp_
    (
        IOobject
        (
            "cp",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("cp", dimensionSet(0, 2, -2, -1, 0), 0),
        calculatedFvPatchScalarField::typeName
    ),
    rho_
    (
        IOobject
        (
            "rho",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("rho", dimensionSet(1, -3, 0, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    )


    // mesh_(alpha1_.mesh() ),

    //--modified_start--//
    // kappa_
    // (
    //     IOobject
    //     (
    //         "kappa",
    //         U_.time().timeName(),
    //         U_.db()
    //     ),
    //     U_.mesh(),
    //     dimensionedScalar("kappa", dimensionSet(1, 1, -3, -1, 0, 0, 0), 0),
    //     calculatedFvPatchScalarField::typeName
    // )
    //--modified_end--//
    {
        //Read fluid properties
        read();

        calcNu();
        //Initial calculation
        calcLambda();
        calcCp();
        calcRho();
    }


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::my_thermalIncompressibleTwoPhaseMixture::mu() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>::New
    (
        "mu",
        limitedAlpha1*rho1_*nuModel1_->nu()
      + (scalar(1) - limitedAlpha1)*rho2_*nuModel2_->nu()
    );
}


Foam::tmp<Foam::volScalarField>
Foam::my_thermalIncompressibleTwoPhaseMixture::lambda() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>::New
    (
        "lambda",
        limitedAlpha1*lambdaModel1_->lambda()
      + (scalar(1) - limitedAlpha1)*lambdaModel2_->lambda()
    );
}

Foam::tmp<Foam::volScalarField>
Foam::my_thermalIncompressibleTwoPhaseMixture::cp() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>::New
    (
        "cp",
        ( cp1_*rho1_*limitedAlpha1 + cp2_*rho2_*(scalar(1) - limitedAlpha1) )
       /( rho1_*limitedAlpha1 + rho2_*(scalar(1) - limitedAlpha1) )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::my_thermalIncompressibleTwoPhaseMixture::rho() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>::New
    (
        "rho",
        limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::my_thermalIncompressibleTwoPhaseMixture::alpha() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    //Calculate average alpha
    return tmp<volScalarField>::New
    (
        "alpha",
        lambda()/(cp()
       *(rho1_*limitedAlpha1 + rho2_*(scalar(1) - limitedAlpha1) ) )
    );
}


Foam::tmp<Foam::scalarField>
Foam::my_thermalIncompressibleTwoPhaseMixture::mu(const label patchI) const
{

    return mu()().boundaryField()[patchI];
}


Foam::tmp<Foam::surfaceScalarField>
Foam::my_thermalIncompressibleTwoPhaseMixture::muf() const
{
    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );

    return tmp<surfaceScalarField>::New
    (
        "muf",
        alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
      + (scalar(1) - alpha1f)*rho2_*fvc::interpolate(nuModel2_->nu())
    );
}

Foam::tmp<Foam::surfaceScalarField>
Foam::my_thermalIncompressibleTwoPhaseMixture::lambdaf() const
{
    const surfaceScalarField alpha1f
    (
    	min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );

    return tmp<surfaceScalarField>::New
    (
      "lambdaf",
      (
          alpha1f*fvc::interpolate(lambdaModel1_->lambda())
        + (scalar(1) - alpha1f)
        *fvc::interpolate(lambdaModel2_->lambda())
      )
    );
}

//- Return the face-interpolated VOF function
Foam::tmp<Foam::surfaceScalarField>
Foam::my_thermalIncompressibleTwoPhaseMixture::alphaf() const
{
    return tmp<surfaceScalarField>::New
    (
      "alphaf",
      fvc::interpolate(alpha())
    );
}

//- Return the face-interpolated thermal capacity
Foam::tmp<Foam::surfaceScalarField>
Foam::my_thermalIncompressibleTwoPhaseMixture::cpf() const
{
    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );

    return tmp<surfaceScalarField>::New
    (
      "cpf",
      fvc::interpolate(cp())
    );
}

Foam::tmp<Foam::surfaceScalarField>
Foam::my_thermalIncompressibleTwoPhaseMixture::nuf() const
{
    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );

    return tmp<surfaceScalarField>::New
    (
        "nuf",
        (
            alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
          + (scalar(1) - alpha1f)*rho2_*fvc::interpolate(nuModel2_->nu())
        )/(alpha1f*rho1_ + (scalar(1) - alpha1f)*rho2_)
    );
}

//- Not used with phase change model
Foam::tmp<Foam::surfaceScalarField>
Foam::my_thermalIncompressibleTwoPhaseMixture::kappaf() const
{
    const surfaceScalarField alpha1f
    (
    	min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );

    return tmp<surfaceScalarField>::New
    (
      "kappaf",
      (
      	alpha1f*rho1_*cp1_*(1/Pr1_)
      	*fvc::interpolate(nuModel1_->nu())
      	+ (scalar(1) - alpha1f)*rho2_*cp2_
      	*(1/Pr2_)*fvc::interpolate(nuModel2_->nu())
      )
    );
}

bool Foam::my_thermalIncompressibleTwoPhaseMixture::read()
{
    if (regIOobject::read())
    {
        if
        (
            nuModel1_().read
            (
                subDict(phase1Name_ == "1" ? "phase1": phase1Name_)
            )
         && nuModel2_().read
            (
                subDict(phase2Name_ == "2" ? "phase2": phase2Name_)
            )
        && lambdaModel1_().read
           (
               subDict(phase1Name_ == "1" ? "phase1": phase1Name_)
           )
        && lambdaModel2_().read
          (
              subDict(phase2Name_ == "2" ? "phase2": phase2Name_)
          )
        )
        {
            nuModel1_->viscosityProperties().readEntry("rho", rho1_);
            nuModel2_->viscosityProperties().readEntry("rho", rho2_);
            nuModel1_->viscosityProperties().lookup("cp") >> cp1_;
            nuModel2_->viscosityProperties().lookup("cp") >> cp2_;
            //- Not used in phase change model
            nuModel1_->viscosityProperties().lookup("Pr") >> Pr1_;
            nuModel2_->viscosityProperties().lookup("Pr") >> Pr2_;
            //-----------modified_end----------------------//

            return true;
        }
    }

    return false;
}


// ************************************************************************* //
