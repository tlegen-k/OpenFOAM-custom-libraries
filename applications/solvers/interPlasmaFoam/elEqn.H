//- Solve Gauss's law for electric potential and charge
surfaceScalarField rhoEFlux(
   -fvc::interpolate(sigmaEl)
   *fvc::snGrad(Ue)
   *mesh.magSf()
);

fvScalarMatrix rhoEEqn
(
    fvm::ddt(rhoE)
  + fvm::div(phi, rhoE)
  + fvc::div(rhoEFlux)
);

// rhoEEqn.solve(mesh.solver("rhoEEqn"));
rhoEEqn.solve();

//- Solve charge fluid advection equation (Lopez-Herrera et al. 2011)

fvScalarMatrix UeEqn
(
    fvm::laplacian(eps,Ue) == -rhoE
);

// UeEqn.solve(mesh.solver("UeEqn"));
UeEqn.solve();

Info<< "Min(Ue)=" << min(Ue).value()
    << " Max(Ue)=" << max(Ue).value() << endl;


Info<< "Min(rhoE)=" << min(rhoE).value()
    << " Max(rhoE)=" << max(rhoE).value() << endl;
