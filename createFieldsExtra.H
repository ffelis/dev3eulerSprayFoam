// Create extra fields

    dimensionedScalar YtildeTiny("YtildeTiny",dimensionSet(0,0,0,0,0,0,0), 1.0e-07);
    dimensionedScalar ELSAD32Tiny("ELSAD32Tiny",dimensionSet(0,1,0,0,0,0,0),1.0e-05);

// Mass flux terms

    volScalarField Ytilde
    (
        IOobject
        (
            "Ytilde",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1*rho1/rho
    );
    
    volScalarField Ybar
    (
        IOobject
        (
            "Ybar",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho*Ytilde*(1.0-Ytilde)*(1.0/rho1-1.0/rho2)
    );
    
    Ytilde = max(min(Ytilde, 0.999999), 1.0e-7);

    volVectorField fluxYtilde
    (
        IOobject
        (
            "fluxYtilde",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
      0.01*U
    );
    
    volVectorField fluxYtilde_s("fluxYtilde_s", fluxYtilde);

    volVectorField fluxYtilde_turb("fluxYtilde_turb", fluxYtilde);

    volVectorField fluxYtilde_base
    (
        IOobject
        (
            "fluxYtilde_base",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
      0.01*U
    );

    Info<< "Reading field F\n" << endl;
    volVectorField F
    (
        IOobject
        (
            "F",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
      mesh
    );
    
    volVectorField Uibar
    (
        IOobject
        (
            "Uibar",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
      0.01*U
    );

// Mixture pressure gradient

    volVectorField gradp
    (
        IOobject
        (
            "gradp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(p_rgh)
    );

    volVectorField gradprgh
    (
        IOobject
        (
            "gradprgh",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(p_rgh)
    );

    volVectorField gradalpha1
    (
        IOobject
        (
            "gradalpha1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(alpha1)
    );
    

// Pressure solver

    volScalarField io_sigma
    (
        IOobject
        (
            "io_sigma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.01*(Uibar & gradprgh)
    );

    volSymmTensorField io_sigmat
    (
        IOobject
        (
            "io_sigmat",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.01*twoSymm(fluxYtilde*gradprgh)
    );

    volScalarField divU
    (
        IOobject
        (
            "divU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::div(phi) // This is the real expression
    );


    volScalarField io_divutilde_source
    (
        IOobject
        (
            "io_divutilde_source",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::div(phi) // Dummy value
    );

// Liquid and Gas velocity fields

    volVectorField UL
    (
        IOobject
        (
            "UL",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U
    );

    volVectorField UG
    (
        IOobject
        (
            "UG",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U
    );



// IO k equation


    volScalarField io_Pk
    (
        IOobject
        (
            "io_Pk",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        turbulence->R() && fvc::grad(U)
    );


    volSymmTensorField Rke
    (
        IOobject
        (
            "Rke",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        //turbulence->R()
        (2.0/3.0)*(turbulence->k())*symm(tensor::one)
    );

    volScalarField io_rCmu
    (
        IOobject
        (
            "io_rCmu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho/rho1 // Dummy unit consistent value
    );

    volScalarField io_C1
    (
        IOobject
        (
            "io_C1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho/rho1 // Dummy unit consistent value
    );

    volScalarField io_magS
    (
        IOobject
        (
            "io_magS",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt(2*magSqr(dev(symm(fvc::grad(U)))))
    );

// Extra create fields for Rij Equation Budget

    

    volSymmTensorField io_Pij
    (
        IOobject
        (
            "io_Pij",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        twoSymm(turbulence->R() & fvc::grad(U))
    );

    volSymmTensorField io_devR
    (
        IOobject
        (
            "io_devR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dev(turbulence->R())
    );

    volSymmTensorField io_devPij
    (
        IOobject
        (
            "io_devPij",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dev(io_Pij)
    );

    volSymmTensorField io_devSigma
    (
        IOobject
        (
            "io_devSigma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dev(io_sigmat)
    );
    
    volVectorField io_F1
    (
        IOobject
        (
            "io_F1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.0001*rho*(F & fvc::grad(U)) // dummy
    );
    
    volVectorField io_F2
    (
        IOobject
        (
            "io_F2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.0001*rho*(F & fvc::grad(U)) // dummy
    );
    
    volVectorField io_F3
    (
        IOobject
        (
            "io_F3",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.0001*rho*(F & fvc::grad(U)) // dummy
    );
    
    volVectorField io_F4
    (
        IOobject
        (
            "io_F4",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.0001*rho*(F & fvc::grad(U)) // dummy
    );
    
    volScalarField Tcoef
    (
        IOobject
        (
            "Tcoef",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        0.00001*rho*((turbulence->epsilon())/(turbulence->k())) // dummy
    );
    
    // Pressure related fields
    
    
    volVectorField sqrtDiagR
    (
        IOobject
        (
            "sqrtDiagR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        0.0001*U
    );
    
    volVectorField diagR
    (
        IOobject
        (
            "diagR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0,2,-2,0,0,0,0),
        "zeroGradient"
    );
    
    dimensionedScalar intensitymax("intensitymax",dimensionSet(0,1,-1,0,0,0,0),20.0);
    dimensionedScalar intensitymin("intensitymax",dimensionSet(0,1,-1,0,0,0,0),1e-14);
    
    

//mesh,
//       dimensionSet(0,2,-2,0,0,0,0),
//       "zeroGradient"


