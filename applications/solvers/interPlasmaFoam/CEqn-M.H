//- Solve transport equation for species in Marschall's formulation
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

    //- Create Arithmetic mean diffusivity
    //- Used in Marschall's model
    volScalarField DAi
    (
        alpha1*DLi + (1 - alpha1)*DGi
    );

    //- Create flux of specie concentration
    surfaceScalarField phiCi(
     (
    //  		fvc::interpolate((DLi-DGi/Hei)/(alpha1+(1-alpha1)/Hei))
    // 		*fvc::snGrad(alpha1)
        (DLi - DGi/Hei)
       /(fvc::interpolate(alpha1) + (1 - fvc::interpolate(alpha1))/Hei)
       *fvc::snGrad(alpha1)
     ) *mesh.magSf()
     );

    // Additional term from Deising for equation of species transport
    surfaceScalarField phiCurv(
        (
          fvc::interpolate(1/(alpha1 + (1 - alpha1)/Hei)*
          (1/Hei*(DLi - DGi)/(alpha1 + (1 - alpha1)/Hei) - (DLi - DGi/Hei)))
         *fvc::snGrad(alpha1)
        )*mesh.magSf()
     );

    //- Create diffusivity flux
    surfaceScalarField phiD(
    (
        fvc::snGrad(DAi)
    ) * mesh.magSf()
    );


    //- Solve transport equation for specie i
    fvScalarMatrix CiEqn
    (
    	fvm::ddt(Ci) // [1/s]
       +fvm::div(phi, Ci, "div(phi,C)") // [1/m * ()]
       -fvm::laplacian(fvc::interpolate(DAi), Ci, "laplacian(C)")
       -fvm::div(phiD, Ci, "div(phi,C)")
       +fvm::div(phiCi, Ci, "div(phi,C)")
        // -fvm::div(phiCurv, Ci, "div(phi,C)")
      ==multiphaseChemistry.R(Ci,i)
    );

    CiEqn.solve(mesh.solver("Ci"));

    Ci.max(0.0); //- Set lower value of concentration to 0

    //- Output min and mac concentrations
    Info<< "Min(" << specie_names[i] << ")=" << min(Ci).value()
        << " Max(" << specie_names[i] << ")=" << max(Ci).value() << endl;
}
