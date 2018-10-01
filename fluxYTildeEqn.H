//Quasi-Multiphase Eulerian Equation for mass-fluxes

//Droplet's local Reynolds number

//Info<< "Debug 1" << endl;

ELSAD32 = 6.0*alpha1*(1.0-alpha1)/(ELSASigmaMin+ELSASigmaPrime+ELSATiny);
volVectorField gF(F-fluxYtilde_base);
gF = min(gF, sqrtDiagR*Ytilde*(1.0-Ytilde));
gF = max(gF, -sqrtDiagR*Ytilde*(1.0-Ytilde));
volScalarField Red("Red", (mag(gF)/(1.0/(Ytilde*(1.0-Ytilde))))*(ELSAD32/nu2));
Red = max(min(Red, 100000.0), 1.0);
//volScalarField Tcoef("Tcoef", 18.0*rho2*nu2*alpha1*(1.0+0.15*pow(Red,0.687))/(ELSAD32*ELSAD32+ELSAD32Tiny*ELSAD32Tiny));
Tcoef = 18.0*rho2*nu2*alpha1*(1.0+0.15*pow(Red,0.687))/(ELSAD32*ELSAD32+ELSAD32Tiny*ELSAD32Tiny);

{
    fvVectorMatrix FEqn
    (
        fvm::ddt(rho, F)
      + fvm::div(rhoPhi, F)
      - fvm::laplacian
        (
            volScalarField("Dab",
            rho*(Dab + Ftab*turbulence->nut())),
            F
        )
      ==
      - CF1*rho*(F & fvc::grad(U))
      - CF2*(rho1*rho2/rho)*(turbulence->R() & fvc::grad(alpha1))
      - CF3*Ybar*fvc::grad(p_rgh)
      + CF4*Tcoef*(1.0/(Ytilde*(1.0-Ytilde)+YtildeTiny))*fluxYtilde_base
      - fvm::Sp(CF4*Tcoef*(1.0/(Ytilde*(1.0-Ytilde)+YtildeTiny)), F)
    );

    FEqn.relax();
    FEqn.solve();
    
    F.internalField() = min(F, sqrtDiagR*Ytilde*(1.0-Ytilde));
    F.internalField() = max(F, -sqrtDiagR*Ytilde*(1.0-Ytilde));

}
