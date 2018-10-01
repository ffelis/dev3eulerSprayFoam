/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    dev3eulerSprayFoam

Description
    Solver for mixing 2 incompressible fluids.
    With an specific surface transport equation model.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "incompressibleTwoPhaseMixture.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        twoPhaseProperties.correct();

        // --- Activate divergence into pressure equation ---
        activdiv = twoPhaseProperties.lookup("activdiv");
        // --------------------------------------------------

        // --- Activate on-the-fly diffusion model change ---
        turbdiff = twoPhaseProperties.lookup("turbdiff");
        Cy = twoPhaseProperties.lookup("Cy");
        Cy2 = twoPhaseProperties.lookup("Cy2");
        Cy3 = twoPhaseProperties.lookup("Cy3");
        Cy4 = twoPhaseProperties.lookup("Cy4");
        // --------------------------------------------------

        //#include "alpharhoRDiffusionEqn.H"

        // Update divergence *new*
        //#include "updatedivutilde.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "alpharhoRDiffusionEqn.H"

            // Update divergence *new*
            #include "updatedivutilde.H"

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                // --- update uprim2.gradp for turbulence source *new*
                #include "updateup2.H"
                turbulence->correct();
            }
        }

        #include "updateiofields.H"

        // --- Specific liquid-gaz surface transport equation
        #include "sigmaPrimeEqn.H"

        Info<< "U values:"
            << "  Min(U) = " << min(mag(U)).value()
            << "  Max(U) = " << max(mag(U)).value()
            << "  Min(F) = " << min(mag(F)).value()
            << "  Max(F) = " << max(mag(F)).value()
            << "  Act Y-Tilde = " << activdiv.value()
            << "  Y-Flux mode = " << turbdiff.value()
            << endl;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
