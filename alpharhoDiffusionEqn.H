{
    fvScalarMatrix alpha1Eqn
    (
        fvm::ddt(alpha1)
      + fvm::div(phi, alpha1)
      - fvm::laplacian
        (
            volScalarField("Dab", 
            (rho2/rho)*(Dab + alphatab*turbulence->nut())),
            alpha1
        )
    );

    alpha1Eqn.relax();
    alpha1Eqn.solve();

    alpha1 = max(min(alpha1, alphamax), alphamin);

    alpha2 = 1.0 - alpha1;
}

rho = alpha1*rho1 + alpha2*rho2;

{
    fvScalarMatrix alpha1Eqn
    (
        fvm::ddt(alpha1)
      + fvm::div(phi, alpha1)
      - fvm::laplacian
        (
            volScalarField("Dab", 
            (rho2/rho)*(Dab + alphatab*turbulence->nut())),
            alpha1
        )
    );

    //alpha1Eqn.relax();
    alpha1Eqn.solve();

    alpha1 = max(min(alpha1, alphamax), alphamin);

    alpha2 = 1.0 - alpha1;
    rhoPhi = alpha1Eqn.flux()*(rho1 - rho2)+phi*rho2;
}

rho = alpha1*rho1 + alpha2*rho2;

Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
