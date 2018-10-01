// Update IO fields

// 1) Rho related

Ytilde = alpha1*rho1/rho;
Ybar = rho*Ytilde*(1.0-Ytilde)*(1.0/rho1-1.0/rho2);
gradalpha1 = fvc::grad(alpha1);
Uibar = -((1.0/rho2-1.0/rho1)*rho)*fluxYtilde_turb;

// 2) Pressure related

gradprgh = fvc::grad(p_rgh);
gradp = fvc::grad(p);

// 3) Velocity related 

divU = fvc::div(phi);
UL = U + (fluxYtilde/Ytilde);
UG = U + (fluxYtilde/(Ytilde-1.0));

// 4) Turbulence related

Rke.internalField() = turbulence->R();
io_devR = dev(turbulence->R());

io_Pk = turbulence->R() && fvc::grad(U);
io_sigma = Uibar & gradprgh;

io_Pij = twoSymm(turbulence->R() & fvc::grad(U));
io_sigmat = twoSymm(Uibar*gradprgh);

io_devPij = dev(io_Pij);
io_devSigma = dev(io_sigmat);

volSymmTensorField S = dev(symm(fvc::grad(U)));
volScalarField S2(2*magSqr(dev(symm(fvc::grad(U)))));
volScalarField magS(sqrt(S2));
volScalarField W
    (
        (2.0*1.41421356)*((S&S)&&S)
       /(
            magS*S2
          + dimensionedScalar("small", dimensionSet(0, 0, -3, 0, 0), SMALL)
        )
    );

volScalarField phis
    (
        (1.0/3.0)*acos(min(max(2.44948974*W, -scalar(1)), scalar(1)))
    );

volScalarField As(2.44948974*cos(phis));
volScalarField Us(sqrt(S2/2.0 + magSqr(skew(fvc::grad(U)))));

volScalarField eta(magS*(turbulence->k())/(turbulence->epsilon()));

io_C1 = max(eta/(scalar(5) + eta), scalar(0.43));
io_rCmu = 1.0/(scalar(4.0) + As*Us*(turbulence->k())/(turbulence->epsilon()));
io_magS = magS;

// 5) Mass Flux F related

io_F1 = rho*(F & fvc::grad(U));
io_F2 = (rho1*rho2/rho)*(turbulence->R() & gradalpha1);
io_F3 = Ybar*gradprgh;
io_F4 = Tcoef*(1.0/(Ytilde*(1.0-Ytilde)+YtildeTiny))*(fluxYtilde_base-F);

