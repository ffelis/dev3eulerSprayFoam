// Update rho values

Ytilde = alpha1*rho1/rho;
Ytilde = max(min(Ytilde, 0.9999999), 1.0e-7);

Ybar = rho*Ytilde*(1.0-Ytilde)*(1.0/rho1-1.0/rho2);

fluxYtilde_base =
      - (alphatab*turbulence->nut())*(rho1*rho2/(rho*rho))*(fvc::grad(alpha1));

if (turbdiff.value() == 1.0)
{
    fluxYtilde.internalField() =
          - (((Cy*rho2*rho1/(rho*rho))*(turbulence->k()/turbulence->epsilon())*
            (turbulence->R() & I)) & (fvc::grad(alpha1)));
    fluxYtilde.boundaryField() = fluxYtilde_base.boundaryField();
}

else if (turbdiff.value() == 2.0)
{
    fluxYtilde.internalField() =
          - (((rho2*rho1/(rho*rho*Cy2))*(turbulence->k()/turbulence->epsilon())*
            (turbulence->R())) & (fvc::grad(alpha1)))
          - (((1.0-Cy3)/(Cy2))*(turbulence->k()/turbulence->epsilon()))*
            (fluxYtilde & fvc::grad(U))
          - (((Ybar*(1.0-Cy4))/(Cy2*rho))*(turbulence->k()/turbulence->epsilon()))*(fvc::grad(p_rgh));
    fluxYtilde.boundaryField() = fluxYtilde_base.boundaryField();
}

else if (turbdiff.value() == 4.0)
{
    fluxYtilde.internalField() =
          - ((((Cym4+(Cpm4/Cym4)*rho*(1.0/rho2-1.0/rho1))*rho2*rho1/(rho*rho))*(turbulence->k()/turbulence->epsilon())*
            (turbulence->R() & I)) & (fvc::grad(alpha1)));
    fluxYtilde.boundaryField() = fluxYtilde_base.boundaryField();
}

else if (turbdiff.value() == 5.0)
{
    fluxYtilde.internalField() = CFblend*fluxYtilde_base + (1.0-CFblend)*F;
    fluxYtilde.boundaryField() = fluxYtilde_base.boundaryField();
}

else
{
    fluxYtilde = fluxYtilde_base;
    //F = fluxYtilde_base;
}

// Bound Mass Flux
Rke.internalField() = turbulence->R();

forAll(diagR,cellI)
{
    diagR[cellI].x() = Foam::mag(Rke[cellI].xx());
    diagR[cellI].y() = Foam::mag(Rke[cellI].yy());
    diagR[cellI].z() = Foam::mag(Rke[cellI].zz());
}

forAll(sqrtDiagR,cellI)
{
    sqrtDiagR[cellI].x() = Foam::sqrt(mag(diagR[cellI].x()));
    sqrtDiagR[cellI].y() = Foam::sqrt(mag(diagR[cellI].y()));
    sqrtDiagR[cellI].z() = Foam::sqrt(mag(diagR[cellI].z()));
}

forAll(sqrtDiagR,cellI)
{
    sqrtDiagR[cellI].x() = max(min(sqrtDiagR[cellI].x(), 20.0), 1e-14);
    sqrtDiagR[cellI].y() = max(min(sqrtDiagR[cellI].y(), 20.0), 1e-14);
    sqrtDiagR[cellI].z() = max(min(sqrtDiagR[cellI].z(), 20.0), 1e-14);
}
  
//sqrtDiagR.component(vector::X) = min(sqrt(mag(Rke.component(tensor::XX))), intensitymax);
//sqrtDiagR.component(vector::Y) = min(sqrt(mag(Rke.component(tensor::YY))), intensitymax);
//sqrtDiagR.component(vector::Z) = min(sqrt(mag(Rke.component(tensor::ZZ))), intensitymax);

fluxYtilde.internalField() = min(fluxYtilde.internalField(), sqrtDiagR.internalField()*Ytilde.internalField()*(1.0-Ytilde.internalField()));
fluxYtilde.internalField() = max(fluxYtilde.internalField(), -sqrtDiagR.internalField()*Ytilde.internalField()*(1.0-Ytilde.internalField()));

//Info<< "Flux values:"
//    << "  Min(Flux) = " << min(sqrtDiagR).value()
//    << "  Max(Flux) = " << max(sqrtDiagR).value()
//    << endl;
