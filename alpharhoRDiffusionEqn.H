// Solve alpha equation

#include "updaterhofields.H"
#include "fluxYTildeEqn.H"

if (turbdiff.value() == 1.0)
{
    fvScalarMatrix alpha1Eqn
    (
        fvm::ddt(alpha1)
      + fvm::div(phi, alpha1)
      - fvm::laplacian
        (
            volScalarField("Dab",
            (Dab + (2.0/3.0)*((Cy*rho2)/(0.09*rho))*turbulence->nut())),
            alpha1
        )
      - fvc::div(
          ((Cy*rho2/rho)*(turbulence->k()/turbulence->epsilon())*
          dev(turbulence->R() & I)) & fvc::grad(alpha1)
                )
    );

    alpha1Eqn.relax();
    alpha1Eqn.solve();

    alpha1 = max(min(alpha1, alphamax), alphamin);

    alpha2 = 1.0 - alpha1;
    rhoPhi = alpha1Eqn.flux()*(rho1 - rho2)+phi*rho2;
}

else if (turbdiff.value() == 2.0)
{
    fvScalarMatrix alpha1Eqn
    (
        fvm::ddt(alpha1)
      + fvm::div(phi, alpha1)
      - fvm::laplacian
        (
            volScalarField("Dab",
            (Dab + (2.0/3.0)*((rho2/Cy2)/(0.09*rho))*turbulence->nut())),
            alpha1
        )
      - fvc::div(
          ((rho2/(rho*Cy2))*(turbulence->k()/turbulence->epsilon())*
          dev(turbulence->R())) & fvc::grad(alpha1)
                )
      - fvc::div(
          (((rho*(1.0-Cy3))/(Cy2*rho1))*(turbulence->k()/turbulence->epsilon()))*
          (fluxYtilde & fvc::grad(U))
                )
      - fvc::div(
          (((Ybar*(1.0-Cy4))/(Cy2*rho1))*(turbulence->k()/turbulence->epsilon()))*
          (fvc::grad(p_rgh))
                )
    );

    alpha1Eqn.relax();
    alpha1Eqn.solve();

    //alpha1.relax();

    alpha1 = max(min(alpha1, alphamax), alphamin);
    alpha2 = 1.0 - alpha1;
    rhoPhi = alpha1Eqn.flux()*(rho1 - rho2)+phi*rho2;
}

else if (turbdiff.value() == 4.0)
{
    fvScalarMatrix alpha1Eqn
    (
        fvm::ddt(alpha1)
      + fvm::div(phi, alpha1)
      - fvm::laplacian
        (
            volScalarField("Dab",
            (Dab + (2.0/3.0)*(((Cym4+(Cpm4/Cym4)*rho*(1.0/rho2-1.0/rho1))*rho2)/(0.09*rho))*turbulence->nut())),
            alpha1
        )
      - fvc::div(
          (((Cy4+(Cpm4/Cym4)*rho*(1.0/rho2-1.0/rho1))*rho2/rho)*(turbulence->k()/turbulence->epsilon())*
          dev(turbulence->R() & I)) & fvc::grad(alpha1)
                )

    );

    alpha1Eqn.relax();
    alpha1Eqn.solve();


    //alpha1.relax();

    alpha1 = max(min(alpha1, alphamax), alphamin);
    alpha2 = 1.0 - alpha1;
    rhoPhi = alpha1Eqn.flux()*(rho1 - rho2)+phi*rho2;
}

else if (turbdiff.value() == 5.0)
{
    fvScalarMatrix alpha1Eqn
    (
       fvm::ddt(alpha1)
     + fvm::div(phi, alpha1)
     - fvm::laplacian
       (
           volScalarField("Dab",
           CFblend*(rho2/rho)*(Dab + alphatab*turbulence->nut())),
           alpha1
       )
     + fvc::div((1.0/rho1)*(1.0-CFblend)*rho*F)
    );

    alpha1Eqn.relax();
    alpha1Eqn.solve();

    alpha1 = max(min(alpha1, alphamax), alphamin);
    alpha2 = 1.0 - alpha1;
    rhoPhi = alpha1Eqn.flux()*(rho1 - rho2)+phi*rho2;
}

else
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
    rhoPhi = alpha1Eqn.flux()*(rho1 - rho2)+phi*rho2;
}

rho = alpha1*rho1 + alpha2*rho2;

Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
