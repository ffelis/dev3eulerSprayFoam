// K equation source term - turbulent mass flux / pressure gradient correlation

#include "updaterhofields.H"

if (tturb.value() == turbdiff.value())
    {
        fluxYtilde_turb = fluxYtilde;
    }
    
else if (tturb.value() == 0.0)

    {
        fluxYtilde_turb.internalField() = CFblend*fluxYtilde_base + (1.0-CFblend)*F;
        fluxYtilde_turb.boundaryField() = fluxYtilde_base.boundaryField();
    }
    
else if (tturb.value() == 10)

    {
        fluxYtilde_turb.internalField() = 1.0*F;
        fluxYtilde_turb.boundaryField() = fluxYtilde_base.boundaryField();
    }

else
    {
        fluxYtilde_turb = fluxYtilde_base;
    }
    
fluxYtilde_turb.internalField() = min(fluxYtilde_turb, sqrtDiagR*Ytilde*(1.0-Ytilde));
fluxYtilde_turb.internalField() = max(fluxYtilde_turb, -sqrtDiagR*Ytilde*(1.0-Ytilde));
    
    
Uibar = -((1.0/rho2-1.0/rho1)*rho)*fluxYtilde_turb;
gradprgh = fvc::grad(p_rgh);

io_sigma = Uibar & gradprgh;
io_sigmat = twoSymm(Uibar*gradprgh);
