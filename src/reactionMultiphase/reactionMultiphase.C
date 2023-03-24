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

#include "reactionMultiphase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// const dataType Foam::reactionMultiphase::staticData();
Foam::label Foam::reactionMultiphase::nUnNamedReactions = 0;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::label Foam::reactionMultiphase::getNewReactionID()
{
    return nUnNamedReactions++;
}
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::reactionMultiphase::reactionStrLeft
(
    OStringStream& reaction
) const
{
    for (label i = 0; i < lhs_.size(); ++i)
    {
        if (i > 0)
        {
            reaction << " + ";
        }
        if (mag(lhs_[i].stoichCoeff - 1) > SMALL)
        {
            reaction << lhs_[i].stoichCoeff;
        }
        reaction << species_[lhs_[i].index];
        if (mag(lhs_[i].exponent - lhs_[i].stoichCoeff) > SMALL)
        {
            reaction << "^" << lhs_[i].exponent;
        }
    }
}

void Foam::reactionMultiphase::reactionStrRight
(
    OStringStream& reaction
) const
{
    for (label i = 0; i < rhs_.size(); ++i)
    {
        if (i > 0)
        {
            reaction << " + ";
        }
        if (mag(rhs_[i].stoichCoeff - 1) > SMALL)
        {
            reaction << rhs_[i].stoichCoeff;
        }
        reaction << species_[rhs_[i].index];
        if (mag(rhs_[i].exponent - rhs_[i].stoichCoeff) > SMALL)
        {
            reaction << "^" << rhs_[i].exponent;
        }
    }
}

//- originally was in private member function
Foam::string Foam::reactionMultiphase::reactionStr
(
    OStringStream& reaction
) const
{
    reactionStrLeft(reaction);
    reaction << " = ";
    reactionStrRight(reaction);
    return reaction.str();
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionMultiphase::reactionMultiphase
(
    speciesTable& species,
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs
)
:
    name_("un-named-reaction-" + Foam::name(getNewReactionID())),
    phase_("unknown-phase"),
    species_(species),
    lhs_(lhs),
    rhs_(rhs)
    {

    }

Foam::reactionMultiphase::reactionMultiphase
(
    const reactionMultiphase& r,
    speciesTable& species
)
:
    name_(r.name() + "Copy"),
    phase_(r.phase() + "Copy"),
    species_(species),
    lhs_(r.lhs_),
    rhs_(r.rhs_)
    {

    }


Foam::reactionMultiphase::specieCoeffs::specieCoeffs
(
    speciesTable& species,
    Istream& is
)
{
    token t(is);
    if (t.isNumber())
    {
        stoichCoeff = t.number();
        is >> t;
    }
    else
    {
        stoichCoeff = 1.0;
    }

    exponent = stoichCoeff;

    if (t.isWord())
    {
        word specieName = t.wordToken();

        const size_t i = specieName.find('^');

        if (i != word::npos)
        {
            exponent = atof(specieName.substr(i + 1).c_str());
            specieName.resize(i);
        }

        if (species.contains(specieName))
        {
            index = species[specieName];
        }
        else
        {
            index = -1;
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Expected a word but found " << t.info()
            << exit(FatalIOError);
    }
}

void Foam::reactionMultiphase::setLRhs
(
    Istream& is,
    speciesTable& species,
    List<specieCoeffs>& lhs,
    List<specieCoeffs>& rhs
)
{
    DynamicList<specieCoeffs> dlrhs;

    while (is.good())
    {
        dlrhs.append(specieCoeffs(species, is));

        if (dlrhs.last().index != -1)
        {
            token t(is);
            if (t.isPunctuation())
            {
                if (t == token::ADD)
                {
                }
                else if (t == token::ASSIGN)
                {
                    lhs = dlrhs.shrink();
                    dlrhs.clear();
                }
                else
                {
                    rhs = dlrhs.shrink();
                    is.putBack(t);
                    return;
                }
            }
            else
            {
                rhs = dlrhs.shrink();
                is.putBack(t);
                return;
            }
        }
        else
        {
            dlrhs.remove();
            if (is.good())
            {
                token t(is);
                if (t.isPunctuation())
                {
                    if (t == token::ADD)
                    {
                    }
                    else if (t == token::ASSIGN)
                    {
                        lhs = dlrhs.shrink();
                        dlrhs.clear();
                    }
                    else
                    {
                        rhs = dlrhs.shrink();
                        is.putBack(t);
                        return;
                    }
                }
            }
            else
            {
                if (!dlrhs.empty())
                {
                    rhs = dlrhs.shrink();
                }
                return;
            }
        }
    }

    FatalIOErrorInFunction(is)
        << "Cannot continue reading reaction data from stream"
        << exit(FatalIOError);
}

Foam::reactionMultiphase::reactionMultiphase
(
    speciesTable& species,
    const dictionary& dict,
    const scalar & A,
    const scalar & beta,
    const scalar & Ta,
    const word & phase
)
:
    name_(dict.dictName()),
    species_(species),
    A_(A),
    beta_(beta),
    Ta_(Ta),
    phase_(phase)
{
    setLRhs
    (
        IStringStream(dict.getString("reaction"))(),
        species_,
        lhs_,
        rhs_
    );

}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::reactionMultiphase::write(Ostream& os) const
{
    OStringStream reaction;
    os.writeEntry("reaction", reactionStr(reaction));
}

Foam::scalar Foam::reactionMultiphase::kf
(
    // const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    scalar ak = A_;
    // Info<< "ak_init = " << ak << endl;
    if (mag(beta_) > VSMALL)
    {
        ak *= pow(T, beta_);
        // Info<< "mag(beta_) > VSMALL : " << ak << endl;
    }

    if (mag(Ta_) > VSMALL)
    {
        ak *= exp(-Ta_/T);
        // Info<< "mag(Ta_) > VSMALL : " << ak << endl;
    }

    return ak;
}

const Foam::speciesTable& Foam::reactionMultiphase::species() const
{
    return species_;
}

// ************************************************************************* //
