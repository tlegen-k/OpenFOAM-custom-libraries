/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "volFields.H"
#include "zeroGradientFvPatchFields.H"
#include "multiphaseChemistryModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- from basicChemistryModel
inline const Foam::fvMesh& Foam::multiphaseChemistryModel::mesh() const
{
    return mesh_;
}


inline Foam::label
Foam::multiphaseChemistryModel::nEqns() const
{
    return nSpecie_ + 2;
}


inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::multiphaseChemistryModel::RR()
{
    return RR_;
}

inline const Foam::List<Foam::reactionMultiphase>&
Foam::multiphaseChemistryModel::reactions() const
{
    return reactions_;
}


inline Foam::label
Foam::multiphaseChemistryModel::nSpecie() const
{
    return nSpecie_;
}

inline Foam::label
Foam::multiphaseChemistryModel::nReaction() const
{
    return nReaction_;
}

inline Foam::scalar
Foam::multiphaseChemistryModel::Treact() const
{
    return Treact_;
}

inline Foam::scalar
Foam::multiphaseChemistryModel::smoothstep(scalar edge0, scalar edge1, scalar x) const
{
    if (x <= edge0){
        return 0.0;
    }
    else if (x >= edge1){
        return 1.0;
    }

    // Scale/bias into [0..1] range
    x = (x - edge0) / (edge1 - edge0);

    return x * x * (3 - 2 * x);
}

inline Foam::scalar&
Foam::multiphaseChemistryModel::Treact()
{
    return Treact_;
}

inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::multiphaseChemistryModel::RR
(
    const label i
) const
{
    return RR_[i];
}

Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::multiphaseChemistryModel::RR
(
    const label i
)
{
    return RR_[i];
}

inline Foam::tmp<Foam::volScalarField>
Foam::multiphaseChemistryModel::T() const
{
    return T_;
}

inline Foam::tmp<Foam::scalarField>
Foam::multiphaseChemistryModel::T(const label patchi) const
{
    return T_.boundaryField()[patchi];
}

inline Foam::tmp<Foam::volScalarField>
Foam::multiphaseChemistryModel::alpha() const
{
    return alpha_;
}

inline Foam::tmp<Foam::scalarField>
Foam::multiphaseChemistryModel::alpha(const label patchi) const
{
    return alpha_.boundaryField()[patchi];
}

inline Foam::tmp<Foam::volScalarField>
Foam::multiphaseChemistryModel::rho() const
{
    return rho_;
}

inline Foam::tmp<Foam::scalarField>
Foam::multiphaseChemistryModel::rho(const label patchi) const
{
    return rho_.boundaryField()[patchi];
}


inline Foam::volScalarField&
Foam::multiphaseChemistryModel::rho()
{
    return rho_;
}

inline PtrList<volScalarField>&
Foam::multiphaseChemistryModel::C()
{
    return C_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
