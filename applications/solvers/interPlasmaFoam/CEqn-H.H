//- Solve transport equation for species in Haroun's formulation
for (label i=0; i<C.size(); i++)
{
    //- Create reference to specie concentration label i
    volScalarField& Ci=C[i];

    //- Create reference to specie Henry constant
    dimensionedScalar& Hei=He[i];

    //- Create reference to specie diffusivity in gas
    dimensionedScalar& DGi=DG[i];

    //- Create reference to specie diffusivity in liquid
    dimensionedScalar& DLi=DL[i];

    //- Create Harmonic mean diffusivity
    //- Used in Haroun's model
    volScalarField DHi
    (
        DLi*DGi/(alpha1*DGi + (1 - alpha1)*DLi)
    );

    //- Flux of specie concentration
    surfaceScalarField phiCi(
     (
        (
             fvc::interpolate(DHi)*(1 - Hei)
            /(fvc::interpolate(alpha1) + (1 - fvc::interpolate(alpha1))*Hei)
        )
 			*fvc::snGrad(alpha1)
     ) *mesh.magSf()
     );

    //- Solve transport equation for specie i
    fvScalarMatrix CiEqn
    (
        fvm::ddt(Ci) // [1/s]
       +fvm::div(phi, Ci, "div(phi,C)") // [1/m * ()]
       -fvm::laplacian(fvc::interpolate(DHi), Ci, "laplacian(C)")
       +fvm::div(phiCi, Ci, "div(phi,C)")
      ==multiphaseChemistry.R(Ci,i)
    );

    CiEqn.solve(mesh.solver("Ci"));

    Ci.max(0.0); //- Set lower value of concentration to 0

    //- Output min and mac concentrations
    Info<< "Min(" << specie_names[i] << ")=" << min(Ci).value()
        << " Max(" << specie_names[i] << ")=" << max(Ci).value() << endl;
}
