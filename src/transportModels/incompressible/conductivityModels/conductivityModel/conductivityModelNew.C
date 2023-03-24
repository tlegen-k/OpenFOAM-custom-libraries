/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 Alex Rattner
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

#include "conductivityModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::conductivityModel> Foam::conductivityModel::New
(
    const word& name,
    const dictionary& conductivityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
{
    const word modelType
    (
        conductivityProperties.get<word>("thermalTransportModel")
    ); //changed from transportModel for laminar fluids

    Info<< "Selecting incompressible thermal transport model "
        << modelType << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown thermalTransportModel type "
            << modelType << nl << nl
            << "Valid thermalTransportModels are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<conductivityModel>
        (
            cstrIter()
            (
                name,
                conductivityProperties,
                U,
                phi
            )
        );
}


// ************************************************************************* //
