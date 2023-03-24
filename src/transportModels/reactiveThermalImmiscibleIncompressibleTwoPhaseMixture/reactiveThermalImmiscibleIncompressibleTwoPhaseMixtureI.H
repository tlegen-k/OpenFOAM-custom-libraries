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

\*---------------------------------------------------------------------------*/
// Inline functions for class reactiveThermalImmiscibleIncompressibleTwoPhaseMixture //

inline const Foam::speciesTable&
Foam::reactiveThermalImmiscibleIncompressibleTwoPhaseMixture::species() const
{
    return species_;
}


inline bool Foam::reactiveThermalImmiscibleIncompressibleTwoPhaseMixture::contains
(
    const word& specieName
) const
{
    return species_.contains(specieName);
}


inline bool Foam::reactiveThermalImmiscibleIncompressibleTwoPhaseMixture::active(label speciei) const
{
    return active_[speciei];
}


inline const Foam::List<bool>& Foam::reactiveThermalImmiscibleIncompressibleTwoPhaseMixture::active() const
{
    return active_;
}


inline void Foam::reactiveThermalImmiscibleIncompressibleTwoPhaseMixture::setActive(label speciei)
{
    active_[speciei] = true;
}


inline void Foam::reactiveThermalImmiscibleIncompressibleTwoPhaseMixture::setInactive(label speciei)
{
    active_[speciei] = false;
}


inline Foam::PtrList<Foam::volScalarField>&
Foam::reactiveThermalImmiscibleIncompressibleTwoPhaseMixture::Y()
{
    return Y_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::reactiveThermalImmiscibleIncompressibleTwoPhaseMixture::Y() const
{
    return Y_;
}


inline Foam::volScalarField& Foam::reactiveThermalImmiscibleIncompressibleTwoPhaseMixture::Y(const label i)
{
    return Y_[i];
}


inline const Foam::volScalarField& Foam::reactiveThermalImmiscibleIncompressibleTwoPhaseMixture::Y
(
    const label i
) const
{
    return Y_[i];
}


inline Foam::volScalarField& Foam::reactiveThermalImmiscibleIncompressibleTwoPhaseMixture::Y
(
    const word& specieName
)
{
    return Y_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::reactiveThermalImmiscibleIncompressibleTwoPhaseMixture::Y
(
    const word& specieName
) const
{
    return Y_[species_[specieName]];
}

// from StandardChemistryModel
inline Foam::label
Foam::reactiveThermalImmiscibleIncompressibleTwoPhaseMixture::nSpecie() const
{
        return nSpecie_;
}

inline Foam::label
Foam::reactiveThermalImmiscibleIncompressibleTwoPhaseMixture::nReaction() const
{
        return nReaction_;
}

inline Foam::label
Foam::reactiveThermalImmiscibleIncompressibleTwoPhaseMixture::nEqns() const
{
    // nEqns = number of species
    return nSpecie_;
}
