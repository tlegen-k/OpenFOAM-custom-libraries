    MRF.correctBoundaryVelocity(U); //in case of Moving Reference Frame update velocity
    //- Maxwell stress tensor
    MST = eps*(fvc::grad(Ue)*fvc::grad(Ue)
        - scalar(0.5)*magSqr(fvc::grad(Ue))*I);

    volVectorField elForce
    (
        fvc::div(MST)
    );

    fvVectorMatrix UEqn // create implicit velocity equation matrix
    (
        fvm::ddt(rho, U) // transient term
      + fvm::div(rhoPhi, U) // convective term
      + MRF.DDt(rho, U) // MRF term
      + turbulence->divDevRhoReff(rho, U) // turbulence term
      - elForce // Electic force term
     ==
        fvOptions(rho, U) //source terms
    );

    UEqn.relax(); //relax matrix

    fvOptions.constrain(UEqn); //apply constrains on matrix

    if (pimple.momentumPredictor()) // solve matrix
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce() //surface tension force
                  - ghf*fvc::snGrad(rho) // gravity terms
                  - fvc::snGrad(p_rgh) //
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U); //correct solution
    }
    // returns max Ux,Uy in liquid phase
    // possibly can be expanded for air phase and for min values as well
    Info << "maxUx: " << max((pos((alpha1-0.5))*U.component(0))).value()
         << " minUx: " << min((pos((alpha1-0.5))*U.component(0))).value() << endl;
    Info << "maxUy: " << max((pos((alpha1-0.5))*U.component(1))).value()
         << " minUy: " << min((pos((alpha1-0.5))*U.component(1))).value() << endl;
