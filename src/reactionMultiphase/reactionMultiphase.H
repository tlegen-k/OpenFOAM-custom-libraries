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
    Foam::reactionMultiphase

Description
    Library to hold multiphase reaction information

SourceFiles
    reactionMultiphaseI.H
    reactionMultiphase.C
    reactionMultiphaseIO.C

\*---------------------------------------------------------------------------*/

#ifndef reactionMultiphase_H
#define reactionMultiphase_H

#include "speciesTable.H"
#include "HashPtrTable.H"
#include "scalarField.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class Istream;
class Ostream;

class reactionMultiphase;
Istream& operator>>(Istream&, reactionMultiphase&);
Ostream& operator<<(Ostream&, const reactionMultiphase&);


/*---------------------------------------------------------------------------*\
                         Class reactionMultiphase Declaration
\*---------------------------------------------------------------------------*/

class reactionMultiphase
{
private:

    // Private data

        //- Name of reaction
        word name_;
        
        //- Name of phase of reaction
        word phase_;

        //- List of specie names present in reaction system
        speciesTable species_;

        scalar A_;
        scalar beta_;
        scalar Ta_;
    // Private Member Functions

        //- Return new reaction ID for un-named reactions
        label getNewReactionID();

protected:

    //- Return string representation of the left of the reaction
    void reactionStrLeft(OStringStream& reaction) const;

    //- Return string representation of the right of the reaction
    void reactionStrRight(OStringStream& reaction) const;

    //- Return string representation of the reactions
    string reactionStr(OStringStream& reaction) const;

public:
    // Public data types
        //- Number of un-named reactions
        static label nUnNamedReactions;
        //- Class to hold the specie index and its coefficients in the
        //  reaction rate expression
        struct specieCoeffs
        {
            label index;
            scalar stoichCoeff;
            scalar exponent;

            specieCoeffs()
            :
                index(-1),
                stoichCoeff(0),
                exponent(1)
            {}

            specieCoeffs(speciesTable& species, Istream& is);

            bool operator==(const specieCoeffs& sc) const
            {
                return index == sc.index;
            }

            bool operator!=(const specieCoeffs& sc) const
            {
                return index != sc.index;
            }

            friend Ostream& operator<<(Ostream& os, const specieCoeffs& sc)
            {
                os  << sc.index << token::SPACE
                    << sc.stoichCoeff << token::SPACE
                    << sc.exponent;
                return os;
            }
        };

    // Constructors

        //- Construct null
         reactionMultiphase() :
         name_(),
         phase_(),
         species_(),
         lhs_(),
         rhs_()
         {}


        //- Construct from components
        reactionMultiphase
        (
            speciesTable& species,
            const List<specieCoeffs>& lhs,
            const List<specieCoeffs>& rhs
        );

        //- Construct as copy given new speciesTable
        reactionMultiphase
        (
            const reactionMultiphase&,
            speciesTable& species
        );

        //- Construct from dictionary
        // NOTE: initReactionThermo is used by solidReaction where there is no
        // need of setting a lhs - rhs thermo type for each reaction. This is
        // needed for mechanism with reversible reactions
        reactionMultiphase
        (
            speciesTable& species,
            const dictionary& dict,
            const scalar & A,
            const scalar & beta,
            const scalar & Ta,
            const word & phase
        );
        //- Construct and return a clone
        virtual autoPtr<reactionMultiphase> clone() const //
        {
            return autoPtr<reactionMultiphase>::New(*this);
        }

    // Selectors

        //- Select null constructed
        static autoPtr<reactionMultiphase> New();


    //- Destructor
    virtual ~reactionMultiphase() = default;


    // Member Functions

        // Access

        inline const word& name() const;
        inline const word& phase() const;
        // Access Arhenius reaction constants
        inline const scalar& A() const;
        inline const scalar& beta() const;
        inline const scalar& Ta() const;

        // - Access to basic components of the reaction
        inline const List<specieCoeffs>& lhs() const;
        inline const List<specieCoeffs>& rhs() const;

        //- Access to specie list
        const speciesTable& species() const;

        //- Construct the left- and right-hand-side reaction coefficients
        void setLRhs
        (
            Istream&,
            speciesTable&,
            List<specieCoeffs>& lhs,
            List<specieCoeffs>& rhs
        );


    // Reaction rate coefficients

        //- Forward rate constant
        virtual scalar kf
        (
            // const scalar p,
            const scalar T,
            const scalarField& c
        ) const;

        virtual void write(Ostream&) const;

    // Member Operators

        // void operator=(const reactionMultiphase&);
        reactionMultiphase& operator=(const reactionMultiphase& other)
        {
            name_ = other.name();
            phase_ = other.phase();
            species_ = other.species();
            lhs_ = other.lhs();
            rhs_ = other.rhs();
            A_ = other.A();
            beta_ = other.beta();
            Ta_ = other.Ta();
            return *this;
        }

    // IOstream Operators

        // friend Istream& operator>>(Istream&, reactionMultiphase&);
        friend Ostream& operator<<(Ostream&, const reactionMultiphase&);
        //- Placed here after SpecieCoeffs structure definition
private:
        //- Specie info for the left-hand-side of the reaction
        List<specieCoeffs> lhs_;

        //- Specie info for the right-hand-side of the reaction
        List<specieCoeffs> rhs_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reactionMultiphaseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
//     #include "reactionMultiphase.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
