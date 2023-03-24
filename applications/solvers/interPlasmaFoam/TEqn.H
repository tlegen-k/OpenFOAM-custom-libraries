{
	// const surfaceScalarField& kappaf = mixture.kappaf()(); //create thermal conductivity field
    surfaceScalarField kappaf = mixture.kappaf()();

    surfaceScalarField cpf = mixture.cpf()();
    // volScalarField Cp = mixture.alpha1()*mixture.cp1() + (scalar(1) - mixture.alpha1()) * mixture.cp2();

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho,T)
        + fvm::div(rhoPhi,T)
        - fvm::laplacian((kappaf/cpf),T)
        // == fvm::Sp(c_source,T);
        == fvOptions(rho, T)
    );
    TEqn.solve(); //solve matrix


    Info << "min/max(kappaf) = " << min(kappaf).value() << ", " << max(kappaf).value() << endl;

    Info << "min/max(cpf) = " << min(cpf).value() << ", " << max(cpf).value() << endl;

    Info << "min/max(T) = " << min(T).value() << ", " << max(T).value() << endl;
}
