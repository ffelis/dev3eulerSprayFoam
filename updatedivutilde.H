// Update alpha1-rho dependant values

#include "updaterhofields.H"



if (tmod.value() == turbdiff.value())
    {
        fluxYtilde_s = fluxYtilde;
    }
    
else if (tmod.value() == 0.0)

    {
        fluxYtilde.internalField() = CFblend*fluxYtilde_base + -(1.0-CFblend)*F;
        fluxYtilde.boundaryField() = fluxYtilde_base.boundaryField();
    }

else
    {
        fluxYtilde_s = fluxYtilde_base;
    }

// Continuity-equation based divergence source term
io_divutilde_source = ((1.0/rho1)-(1.0/rho2))*fvc::div(-rho*fluxYtilde_s);
