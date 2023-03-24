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

Application
    interPlasmaFoam

Group
    grpMultiphaseSolvers

Description
    Solver for two incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

    Equations being solved:

    Species transport equation is solved using 3 formulations:
    1. Haroun's using Harmonic mean.
    2. Marschall's using Arithmetic mean.
    3. Deising's unified approach using Arithmetic mean.

    + 	Homo- and hetero- reactions.

    Energy conservation equation using Enthalpy formulation.

    Gauss'law for electrical potential and charges.

    Electric charge density advection equation.

    Momentum equation is updated to account for electrical forces roles on interface.

Authors (in chronological order):

    Tlegen Kamidollayev (tlegen_kamidollayev@uml.edu, main developer, University of Massachusetts Lowell)
    Juan Pablo Trelles (juan_trelles@uml.edu, refactoring - PI, University of Massachusetts Lowell)


Contact persons:

    Tlegen Kamidollayev, tlegen_kamidollayev@uml.edu, main developer, University of Massachusetts Lowell REng Lab

    Juan Pablo Trelles, Juan_Trelles@uml.edu, principal investigator, University of Massachusetts Lowell REng Lab

Affiliations:

    Affiliation X) REng Lab
                   Juan Pablo Trelles:
                   Mechanical Engineering Department
                   University of Massachusetts Lowell, United States

Acknowledgement:
    Financed by the U.S. Department of Energy through award DE-SC0018230 and the U.S. National Science Foundation through Award CBET-1552037.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H" //including basic libraries required for CFD
#include "dynamicFvMesh.H" //abstract base class for geometry and topology changing
#include "CMULES.H" //Multidimensional universal limiter for explicit corrected implicit solution.
#include "EulerDdtScheme.H"  //Euler time derivative shceme
#include "localEulerDdtScheme.H" // local Euler time derivative scheme
#include "CrankNicolsonDdtScheme.H" //Crank-Nicolson time derivative scheme
#include "subCycle.H" //perform subCycleTime on a field
#include "my_thermalImmiscibleIncompressibleTwoPhaseMixture.H" //library for immiscible incompressible two phase mixture
#include "turbulentTransportModel.H" // turbulent transport model library
#include "pimpleControl.H" // PIMPLE algorithm library
#include "fvOptions.H" // finitie volume options
#include "CorrectPhi.H" //correct fluxes
#include "fvcSmooth.H" //Provides functions smooth spread and sweep which use the FaceCellWave algorithm to smooth and redistribute the first field argument
#include "hashedWordList.H"
#include "reactionMultiphase.H" // multiphase reactions
#include "specie.H" // specie
#include "multiphaseChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, non-isothermal immiscible fluids"
        " using VOF phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing. Energy conservation equation, species transport" ", chemical reaction, electrostatic equations are modeled."
    );

    #include "postProcess.H" //Execute application functionObjects to post-process existing results.

    #include "addCheckCaseOptions.H" //check case set-up
    #include "setRootCaseLists.H" //check root case
    #include "createTime.H" //creates time
    #include "createDynamicFvMesh.H" //create dynamic mesh
    #include "initContinuityErrs.H" //initialize and declare the cumulative continuity error
    #include "createDyMControls.H" //create controls for dynamic mesh
    #include "createFields.H" //create fields
    #include "createAlphaFluxes.H" //create phase fraction fluxes
    #include "initCorrectPhi.H" //initialize correct fluxes
    #include "createUfIfPresent.H" // Creates and initialises the velocity field Uf if required.
    // check turbulence
    turbulence->validate();
    // if trasient equations being solved
    if (!LTS)
    {
        #include "CourantNo.H" // creates, calculates and outputs the mean and max Courant Numbers
        #include "setInitialDeltaT.H" //set initial timestep corresponding to the timestep adjustment algorithm in setDeltaT but only if it would reduce the timestep
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
    // main time advancing loop
    while (runTime.run())
    {
        #include "readDyMControls.H" //read dynamic mesh controls

        if (LTS)
        {
            #include "setRDeltaT.H" //initialize reciprocal time scale field
        }
        else
        {
            #include "CourantNo.H" //create and calculate Courant number
            #include "alphaCourantNo.H" //create VoF Courant number
            #include "setDeltaT.H" //reset timestep to maintain Co Number. Reduction is immediate, but increase is damped to avoid unstable oscillations
        }
        // advance time step
        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // update mesh
                mesh.update();
                // if mesh topology changes use
                if (mesh.changing())
                {
                    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H" //correct fluxes

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H" //calculates and outputs the mean and max Courant Numbers
                    }
                }
            }

            #include "alphaControls.H" //corrections of VoF field
            #include "alphaEqnSubCycle.H" //VOF equation subcycles

            //Update turbulence and two phase properties
            mixture.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

            // #include "elEqn.H" // solve ehd equations
            #include "UEqn.H" // solve momentum equation
            // #include "TEqn.H" // solve temperature equation

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H" //solve pressure equation
            }

            if (pimple.turbCorr()) // if have turbulence -> solve
            {
                turbulence->correct();
            }
        }

        // For now, the energy equation is only 1-way coupled with the
        // momentum/pressure equations, so it can be solved explicitly, and
        // separately here

        // Energy equation with phase change
        // #include "EEqn.H"

        // Solve temperature based energy equation without phase change
        // #include "TEqn.H"

        //- Moved out of PIMPLE loop
        if (transport_model_number == 0)
            #include "CEqn-H.H" // solve Haroun's formulation for species transport equation
        else if (transport_model_number == 1)
            #include "CEqn-M.H" // solve Marschall's formulation
        else if (transport_model_number == 2)
            #include "CEqn-D.H" // solve Deising's formulation
        else
            ;

        runTime.write(); // write data to output file

        runTime.printExecutionTime(Info); // write current time
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
