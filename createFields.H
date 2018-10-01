    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

    Info<< "Reading transportProperties\n" << endl;
    incompressibleTwoPhaseMixture twoPhaseProperties(U, phi);

    volScalarField& alpha1(twoPhaseProperties.alpha1());
    volScalarField& alpha2(twoPhaseProperties.alpha2());

    const dimensionedScalar& rho1 = twoPhaseProperties.rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties.rho2();

    dimensionedScalar nu1(twoPhaseProperties.subDict("water").lookup("nu"));
    dimensionedScalar nu2(twoPhaseProperties.subDict("air").lookup("nu"));

    dimensionedScalar Dab(twoPhaseProperties.lookup("Dab"));

    // Read the parameters for Sigma equation
    dimensionedScalar alphatr(twoPhaseProperties.lookup("alphatr"));
    dimensionedScalar Dr(twoPhaseProperties.lookup("Dr"));
    dimensionedScalar a0(twoPhaseProperties.lookup("a0"));
    dimensionedScalar a1(twoPhaseProperties.lookup("a1"));
    dimensionedScalar a2(twoPhaseProperties.lookup("a2"));
    dimensionedScalar Cr(twoPhaseProperties.lookup("Cr"));
    dimensionedScalar sigma(twoPhaseProperties.lookup("sigma"));
    dimensionedScalar RDmin(twoPhaseProperties.lookup("RDmin"));
    dimensionedScalar RDmax(twoPhaseProperties.lookup("RDmax"));
    dimensionedScalar alphamin(twoPhaseProperties.lookup("alphamin"));
    dimensionedScalar alphamax(twoPhaseProperties.lookup("alphamax"));

    dimensionedScalar activdiv(twoPhaseProperties.lookup("activdiv"));

    // Diffusion model to use in :
    dimensionedScalar tmod(twoPhaseProperties.lookup("tmod")); //pressure
    dimensionedScalar tturb(twoPhaseProperties.lookup("tturb")); //Turbulence

    // Read the reciprocal of the turbulent Schmidt number
    dimensionedScalar alphatab(twoPhaseProperties.lookup("alphatab"));
    dimensionedScalar turbdiff(twoPhaseProperties.lookup("turbdiff")); //YTilde-Eq
    dimensionedScalar Cy(twoPhaseProperties.lookup("Cy"));
    dimensionedScalar Cy2(twoPhaseProperties.lookup("Cy2"));
    dimensionedScalar Cy3(twoPhaseProperties.lookup("Cy3"));
    dimensionedScalar Cy4(twoPhaseProperties.lookup("Cy4"));

    dimensionedScalar Cym4(twoPhaseProperties.lookup("Cym4"));
    dimensionedScalar Cpm4(twoPhaseProperties.lookup("Cpm4"));

    dimensionedScalar Ftab(twoPhaseProperties.lookup("Ftab"));
    dimensionedScalar CF1(twoPhaseProperties.lookup("CF1"));
    dimensionedScalar CF2(twoPhaseProperties.lookup("CF2"));
    dimensionedScalar CF3(twoPhaseProperties.lookup("CF3"));
    dimensionedScalar CF4(twoPhaseProperties.lookup("CF4"));

    dimensionedScalar CFblend(twoPhaseProperties.lookup("CFblend"));

    // Need to store rho for ddt(rho, U)
    volScalarField rho("rho", alpha1*rho1 + alpha2*rho2);
    rho.oldTime();

    // Mass fraction for new alpha equation
    //volScalarField alphamass("alphamass", alpha1*rho1/rho);

    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha1 solution before it is used in the U equation.
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho1*phi
    );

    // Construct incompressible turbulence model ++ rho variable
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, twoPhaseProperties)
    );

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    #include "createFieldsExtra.H"
    #include "createELSA.H"
