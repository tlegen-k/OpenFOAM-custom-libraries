
#include "createRDeltaT.H"
#include "DynamicList.H" // should I include it here?

Info<< "Reading field p_rgh\n" << endl; //create custom pressure field
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl; //create velocity field
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field T\n" << endl; //create temperature field
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        // IOobject::MUST_READ_IF_MODIFIED,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

#include "createPhi.H" //create fluxes of fields


Info<< "Reading transportProperties\n" << endl; //read transport properties
// thermalImmiscibleIncompressibleTwoPhaseMixture mixture(U, phi, T);
my_thermalImmiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1()); //create VoF of 1st phase
volScalarField& alpha2(mixture.alpha2()); //create VoF of 2nd phase

const dimensionedScalar& rho1 = mixture.rho1(); //set density of 1st phase
const dimensionedScalar& rho2 = mixture.rho2(); //set density of 2nd phase
// set thermal capacities of phases
const dimensionedScalar& cp1 = mixture.cp1();
const dimensionedScalar& cp2 = mixture.cp2();

// Need to store rho for ddt(rho, U)
volScalarField rho //set density of mixture
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi //set flux of density of mixture
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H" //read gravitational acceleration
#include "readhRef.H" //read reference height
#include "gh.H" //calculate g*h


volScalarField p //create pressure field
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

//- MULES compressed flux is registered in case scalarTransport FO needs it.
surfaceScalarField alphaPhiUn
(
    IOobject
    (
        "alphaPhiUn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(phi.dimensions(), Zero)
);

//- Create Moving Reference Frame options
#include "createMRF.H"

//- Create finite volume options
#include "createFvOptions.H"

//Create kEff here so it can be accessed globally
surfaceScalarField kEff
(
    IOobject
    (
        "kEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate( mixture.lambda() )
  + fvc::interpolate( rho*mixture.cp()*turbulence->nut() )
);


Info<< "Creating chemistry model.\n" << endl;

IOdictionary chemDict
(
    IOobject
    (
        "chemistryProperties",               //dictionary name
        runTime.constant(),     //dict is found in constant
        mesh,                   //registry for the dict
        IOobject::MUST_READ,    //must exist, otherwise failure
        IOobject::NO_WRITE      //dict is only read by the solver
    )
);


//- Look for species in dictionary values
const wordList specie_names
(
    chemDict.lookup("species")
);

Info<< "Choosing species transport model.\n" << endl;

//- Look for species transport formulation model
const word transport_model
(
    transportProperties.lookup("transportModel")
);

//- Choosing species transport formulation
//- 0 - Haroun's
//  1 - Marschall's
//  2 - Deising's
int transport_model_number;
if (transport_model == "Haroun")
{
    transport_model_number = 0;
    Info<< "Haroun's species transport formulation is chosen.\n" <<endl;
}
else if (transport_model == "Marschall")
{
    transport_model_number = 1;
    Info<< "Marschall's species transport formulation is chosen.\n" <<endl;
}
else if (transport_model == "Deising")
{
    transport_model_number = 2;
    Info<< "Deising's species transport formulation is chosen.\n" <<endl;
}
else {
    Info<< "Transport model not specified in constant/transportProperties.\n"   << endl;
}

//- Read and print number of species
scalar n = specie_names.size();
Info<< "Number of species = " << n << "\n" << endl;
Info<< "Species are: ";
forAll(specie_names, i)
{
    const word& current_specie = specie_names[i];
    Info<< current_specie << ", ";
}

//- Look for elements in dictionary
const wordList elements
(
    chemDict.lookup("elements")
);
//- Create list of species
PtrList<volScalarField> C(specie_names.size());
//- Create list of species volume fratcions
PtrList<volScalarField> Y(specie_names.size());
//- Create list with Henry constants for species
PtrList<dimensionedScalar> He(specie_names.size());
//- Create list with diffusivity in liquid
PtrList<dimensionedScalar> DL(specie_names.size());
//- Create list with diffusivity in gas
PtrList<dimensionedScalar> DG(specie_names.size());
//- Create list with Molar weights of species
PtrList<dimensionedScalar> W(specie_names.size());
//- Create list with dimensionless Molar weights of species
List<scalar> W_scalar(specie_names.size());

forAll(C, i)
{
    const word& currentSpecie = specie_names[i];

    Info << "Reading " << currentSpecie << endl;
    C.set //-HM ‘set or hook’ (OF-version dependent)
    (
        i,
        new volScalarField
        (
            IOobject
            (
                // "C" + currentSpecie,
                currentSpecie,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    const dictionary& subDict = chemDict.subDict(currentSpecie);

    He.set
    (
         i,
         new dimensionedScalar(subDict.lookup("He"))
    );
    Info << He[i] << endl;
    DL.set
    (
         i,
         new dimensionedScalar(subDict.lookup("DL"))
    );
    Info << DL[i] << endl;
    DG.set
    (
         i,
         new dimensionedScalar(subDict.lookup("DG"))
    );
    Info << DG[i] << endl;
    W.set
    (
        i,
        new dimensionedScalar(subDict.lookup("W"))
    );
    Info << W[i] << "\n" << endl;
    W_scalar[i] = W[i].value();
}

forAll(W_scalar, i)
{
    Info<< "Molecular weight of " << specie_names[i]
        << " is " << W_scalar[i] << endl;
}

//- Read and print number of elements
scalar n_elem = elements.size();
Info<< "\nNumber of elements = " << n_elem << endl;
Info<< "Elements are: ";
forAll(elements, i)
{
    const word& current_element = elements[i];
    Info<< current_element;
}

//- Read and print reactions from dictionary
const wordList reaction_names
(
    chemDict.lookup("reaction_names")
);
Info<< "\nNumber of reactions = " << reaction_names.size() << endl;
Info<< "Reactions are: " << endl;
forAll(reaction_names, i)
{
    const word& current_reaction = reaction_names[i];
    Info<< current_reaction << endl;
}

//- Create list with A, beta, Ta constants for reactions
PtrList<dimensionedScalar> A(reaction_names.size());
PtrList<dimensionedScalar> beta(reaction_names.size());
PtrList<dimensionedScalar> Ta(reaction_names.size());
stringList reactions_string(reaction_names.size());
stringList reactions_phases_string(reaction_names.size());

// - List of specie names present in reaction system
hashedWordList species_(specie_names);

speciesTable species_table(specie_names);

List<reactionMultiphase> reactions(reaction_names.size());

forAll(reaction_names, i)
{
    const word& currentReaction = reaction_names[i];

    // Info << currentReaction << endl;
    Info<< "Reading reaction info for: " << currentReaction << endl;

    const dictionary& subDict = chemDict.subDict(currentReaction);

    // Reading Arrhenius reaction rate constants for each reaction
    A.set
    (
         i,
         new dimensionedScalar(subDict.lookup("A"))
    );

    beta.set
    (
         i,
         new dimensionedScalar(subDict.lookup("beta"))
    );

    Ta.set
    (
         i,
         new dimensionedScalar(subDict.lookup("Ta"))
    );

    reactions_string[i] = string(subDict.lookup("reaction"));
    Info<< reactions_string[i] << endl;

    reactions_phases_string[i] = string(subDict.lookup("phase"));

    // Creating model for multiphase reactions list
    reactionMultiphase reaction_i(species_table, subDict,  A[i].value(),        beta[i].value(), Ta[i].value(), reactions_phases_string[i]);

    reactions[i] = reaction_i;

}

Info << "Testing List<reactionMultiphase>:" << endl;
Info<< "List size after filling " << reactions.size() << endl;
for(int i=0; i < reactions.size(); i++)
// forAll(reactions, i)
{
    Info<< "For " << reactions[i] << endl;

    Info<< "Lhs:" << endl;
    Info<< reactions[i].lhs() << endl;
    Info<< "Rhs:" <<  endl;
    Info<< reactions[i].rhs() << endl;

    Info<< "A = " << reactions[i].A() << endl;
    Info<< "beta = " << reactions[i].beta() << endl;
    Info<< "Ta = " << reactions[i].Ta() << endl;
    Info<< "Phase = " << reactions[i].phase() << endl;
}

// TODO: add alpha1 to the constructor
Info<< "Create class multiphaseChemistryModel" << endl;
multiphaseChemistryModel multiphaseChemistry
(
    reactions, C, W_scalar, T, alpha1, rho, mesh
);

//- Read electrical properties and create fields
Info<< "Reading electrical properties\n" << nl << endl;
word phase1Name_(transportProperties.get<wordList>("phases")[0]);
word phase2Name_(transportProperties.get<wordList>("phases")[1]);
// Electrical Properties
dimensionedScalar eps1
(
    transportProperties.subDict(phase1Name_).lookup("eps")
);

dimensionedScalar eps2
(
    transportProperties.subDict(phase2Name_).lookup("eps")
);

dimensionedScalar sigmaEl1
(
    transportProperties.subDict(phase1Name_).lookup("sigmaEl")
);

dimensionedScalar sigmaEl2
(
    transportProperties.subDict(phase2Name_).lookup("sigmaEl")
);

Info<<"Found:\neps1=" << eps1 << "\neps2=" << eps2 << "\nsigmaEl1=" << sigmaEl1 << "\nsigmaEl2=" << sigmaEl2 << endl;

//- Averaged electrical conductivity
volScalarField sigmaEl
(
    IOobject
    (
        "sigmaEl",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    //- Harmonic mean
    scalar(1) / (alpha1/sigmaEl1 + (scalar(1) - alpha1)/sigmaEl2),
    //- Arithmetic mean
    // alpha1*sigmaEl1 + (scalar(1) - alpha1)*sigmaEl2,
    alpha1.boundaryField().types()
);
sigmaEl.oldTime();

//- Averaged electrical permitivity
volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    //- Harmonic mean
    // scalar(1) / (alpha1/eps1 + (scalar(1) - alpha1)/eps2)
    //- Arithmetic mean
    alpha1*eps1 + (scalar(1) - alpha1)*eps2,
    alpha1.boundaryField().types()
);
eps.oldTime();

Info<< "Reading field Electric Potential Ue\n" << endl;
//- Electrical potential
volScalarField Ue
(
    IOobject
    (
        "Ue",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field Electric Density rhoE\n" << endl;
//- Electric charge density
volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
//- Electric field
volVectorField E
(
    IOobject
    (
        "E",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(Ue)
);

//- Maxwell stress tensor
volTensorField MST
(
    IOobject
    (
        "MST",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    // eps*(fvc::grad(Ue)*fvc::grad(Ue) - scalar(0.5)*magSqr(fvc::grad(Ue))*I)
    eps*((E*E) - scalar(0.5)*magSqr(E)*I)
);
