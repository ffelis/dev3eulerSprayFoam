    volScalarField ELSARho
    (
        IOobject
        (
            "ELSARho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho
    );

    Info<< "Reading field ELSASigmaPrime\n" << endl;
    volScalarField ELSASigmaPrime
    (
        IOobject
        (
            "ELSASigmaPrime",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField ELSASigmaMin
    (
        IOobject
        (
            "ELSASigmaMin",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        ELSASigmaPrime
/*        mesh,
       dimensionSet(0,-1,0,0,0,0,0),
       "zeroGradient"*/
    );

    volScalarField ELSASigma
    (
        IOobject
        (
            "ELSASigma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        ELSASigmaPrime
/*        mesh,
       dimensionSet(0,-1,0,0,0,0,0),
       "zeroGradient"*/
    );

    volScalarField ELSASigmaEquil
    (
        IOobject
        (
            "ELSASigmaEquil",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        ELSASigmaPrime
/*        mesh,
       dimensionSet(0,-1,0,0,0,0,0),
       "zeroGradient"*/
    );

    volScalarField ELSAD32
    (
        IOobject
        (
            "ELSASD32",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
       dimensionSet(0,1,0,0,0,0,0),
       "zeroGradient"
    );



    volScalarField ELSASigmaRatio
    (
        IOobject
        (
            "ELSASigmaRatio",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
       dimensionSet(0,0,0,0,0,0,0),
       "zeroGradient"
    );

    dimensionedScalar CstELSASigmaMin(twoPhaseProperties.lookup("CstELSASigmaMin"));

    dimensionedScalar ELSASurfaceTensionCoef(twoPhaseProperties.lookup("ELSASurfaceTensionCoef"));
    // model for equilibriume surface density choice between 
    //  Lebas IJMF 2009 ==> webEquilMod=0
    //  and Duret IJMF 2013  ==> webEquilMod=1
    dimensionedScalar webEquilMod(twoPhaseProperties.lookup("webEquilMod"));

    dimensionedScalar ELSAWeberEquil(twoPhaseProperties.lookup("ELSAWeberEquil"));

    dimensionedScalar CTurbELSASigmaPrime(twoPhaseProperties.lookup("CTurbELSASigmaPrime"));

    dimensionedScalar ELSATiny("ELSATiny",dimensionSet(0,-1,0,0,0,0,0),1.e-1);
    dimensionedScalar kTiny("kTiny",dimensionSet(0,2,-2,0,0,0,0),1.e-1);
    dimensionedScalar sigmaTiny("sigmaTiny",dimensionSet(0,-1,0,0,0,0,0),1.e-1);
    
    //Info<< "Debug 1" << endl;
