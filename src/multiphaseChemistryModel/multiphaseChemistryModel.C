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

\*---------------------------------------------------------------------------*/

#include "multiphaseChemistryModel.H"
#include "UniformField.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Constructor using VOF alpha field
Foam::multiphaseChemistryModel::multiphaseChemistryModel
(
    List<reactionMultiphase>&   reactions,
    // PtrList<volScalarField>& Y,
    PtrList<volScalarField>&    C,
    List<scalar>                W,
    volScalarField&             T,
    volScalarField&             alpha,
    volScalarField&             rho,
    const fvMesh&               mesh
)
:
    // ODESystem(),
    mesh_(mesh),
    // Y_(Y),
    C_(C),
    W_(W),
    T_(T),
    alpha_(alpha),
    rho_(rho),
    reactions_(reactions),
    nSpecie_(C_.size()),
    nReaction_(reactions_.size()),
    Treact_
    ( 0.0
    ),
    RR_(nSpecie_),
    c_(nSpecie_),
    dcdt_(nSpecie_)
{
    // Create the fields for the chemistry sources
    forAll(RR_, fieldi)
    {
        RR_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RR." + C_[fieldi].name(),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimMoles/dimVolume/dimTime, Zero)
                // dimensionedScalar(dimensionSet(0, 0, -1, 0, 0, 0, 0), Zero)
            )
        );
    }

    Info<< "multiphaseChemistryModel: Number of species = " << nSpecie_
        << " and reactions = " << nReaction_ << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiphaseChemistryModel::
~multiphaseChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- dc/dt = omega, rate of change in concentration, for each species
//- [dc/dt] = [kg / m3 / sec]
void
Foam::multiphaseChemistryModel::omega
(
    const scalarField&  c,
    const scalar&       T,
    const scalar&       alpha,
    // const scalar     p,
    scalarField&        dcdt
) const
{
    scalar pf, cf;
    label lRef;
    //- rate of change = 0 initially
    dcdt = Zero;

    //- Add the reaction rates of each reaction containing the species k
    //- to get the reaction rate of the specie k
    //- wk = Sum(wi * (nuki' - nuki''))
    forAll(reactions_, i)
    {
        const reactionMultiphase& R = reactions_[i];

        //- Return the reaction rate for reaction R and the reference
        //  species and characteristic times
        //- [omega] = [mole / m3 / sec]
        scalar omegai = omega
        (
            R, c, T, alpha, pf, cf, lRef
        );

        // Info<< "omegai = " << omegai << endl;

        //- Loop through species in LHS and RHS of current reaction and
        //- calculate rate of change for current specie
        forAll(R.lhs(), s)
        {
            const label si = R.lhs()[s].index;
            const scalar sl = R.lhs()[s].stoichCoeff;
            dcdt[si] -= sl*omegai;
        }

        forAll(R.rhs(), s)
        {
            const label si = R.rhs()[s].index;
            const scalar sr = R.rhs()[s].stoichCoeff;
            dcdt[si] += sr*omegai;
        }
    }
}

//- Return the reaction rate for iReaction and the reference
//  species and characteristic times
//- [omegaI] = [mole / m3 / sec]
Foam::scalar
Foam::multiphaseChemistryModel::omegaI
(
    const label         index,
    const scalarField&  c,
    const scalar&       T,
    const scalar&       alpha,
    // const scalar     p,
    scalar&             pf,
    scalar&             cf,
    label&              lRef
) const
{
    const reactionMultiphase& R = reactions_[index];
    scalar w = omega(R, c, T, alpha, pf, cf, lRef);
    return(w);
}

//- Return the reaction rate for reaction R and the reference
//  species and characteristic times
//- [omega] = [mole / m3 / sec]
Foam::scalar
Foam::multiphaseChemistryModel::omega
(
    const reactionMultiphase&   R,
    const scalarField&          c,
    const scalar&               T,
    const scalar&               alpha,
    // const scalar             p,
    scalar&                     pf,
    scalar&                     cf,
    label&                      lRef
) const
{
    scalar kf = 0.0;

    //- Forward Arrhenius reaction rate with taking into the accound the phase of reaction specified

    if (R.phase() == "liquid"){
        kf = R.kf(T, c) * smoothstep(0.4, 0.6, alpha);
    }
    else if (R.phase() == "gas"){
        kf = R.kf(T, c) * (1 - smoothstep(0.4, 0.6, alpha));
    }

    pf = 1.0;
    //- Number of species in left side of reaction
    const label Nl = R.lhs().size();
    // const label Nr = R.rhs().size();

    label slRef = 0;
    lRef = R.lhs()[slRef].index;

    pf = kf;
    for (label s = 1; s < Nl; s++)
    {
        const label si = R.lhs()[s].index;

        if (c[si] < c[lRef])
        {
            const scalar exp = R.lhs()[slRef].exponent;
            pf *= pow(max(c[lRef], 0.0), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            const scalar exp = R.lhs()[s].exponent;
            pf *= pow(max(c[si], 0.0), exp);
        }
    }
    cf = max(c[lRef], 0.0);

    {
        const scalar exp = R.lhs()[slRef].exponent;
        if (exp < 1.0)
        {
            if (cf > SMALL)
            {
                pf *= pow(cf, exp - 1.0);
            }
            else
            {
                pf = 0.0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1.0);
        }
    }

    return pf*cf;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::multiphaseChemistryModel::calculateRR
(
    const label ri,
    const label si
) const
{
    scalar pf, cf;
    label lRef;

    tmp<volScalarField::Internal> tRR
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "RR",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimMoles/dimVolume/dimTime, Zero)
        )
    );

    volScalarField::Internal& RR = tRR.ref();

    tmp<volScalarField> trho(this->rho());
    const scalarField& rho = trho();

    const scalarField& T = this->T();
    const scalarField& alpha = this->alpha();

    forAll(rho, celli)
    {
        const scalar Ti = T[celli];
        const scalar alphai = alpha[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Ci = C_[i][celli];
            c_[i] = Ci;
        }

        const scalar w = omegaI
        (
            ri,
            c_,
            Ti,
            alphai,
            pf,
            cf,
            lRef
        );

        RR[celli] = w*W_[si];
    }

    return tRR;
}


void
Foam::multiphaseChemistryModel::calculate()
{
    tmp<volScalarField> trho(this->rho());
    const scalarField& rho = trho();
    const scalarField& T = this->T();
    const scalarField& alpha = this->alpha();


    forAll(rho, celli)
    {
        const scalar Ti = T[celli];
        const scalar alphai = alpha[celli];
        
        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Ci = C_[i][celli];
            c_[i] = Ci;
        }

        omega(c_, Ti, alphai, dcdt_);

        for (label i=0; i<nSpecie_; i++)
        {
            //- Get mass concentration of specie
            RR_[i][celli] = dcdt_[i]*W_[i];
        }
    }
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::multiphaseChemistryModel::R(volScalarField& C, label specieI)
{
    this->calculate();
    // create temporary fvScalarMatrix tSu from C volScalarField
    // with the result of calculation
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(C, dimMoles/dimTime));

    // create pointer to tSu
    fvScalarMatrix& Su = tSu.ref();

    // Add to Su Total reaction rate for specieI
    Su += RR(specieI);

    // return tSu object
    return tSu;
}


// ************************************************************************* //
