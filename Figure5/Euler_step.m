function [Cci_output, Cmi_output, Cma_output, Gc_output, Gm_output, GmCma_output, Ii_output, Ia_output] = Euler_step(dt, Cci, Cmi, Cma, Gc, Gm, GmCma, Ii, Ia, D_Cci, D_Cmi, D_Cma, D_Gc, D_Gm, D_GmCma, D_Ii, D_Ia, Flk,...
    k8a,k8b,k10,k2,k1,k9a,k9b,k5a,k5b,k4,k3,k6,k7)
% Reaction rate calculation

Cci_output = Cci + dt*(k5b*Cmi - k5a*Cci);
Cmi_output = Cmi + dt*((k2 + k3.*Ia).*Cma + k5a*Cci - (k10.*Gm + k1.*GmCma + k5b).*Cmi);
Cma_output = Cma + dt*((k10.*Gm + k1.*GmCma).*Cmi + k9b*GmCma - (k2 + k9a.*Gm + k7.*Gc + k3.*Ia).*Cma);
Gc_output = Gc + dt*(k8b*Gm - (k8a + k7.*Cma).*Gc);
Gm_output = Gm + dt*(k8a*Gc + k9b*GmCma - (k8b + k9a.*Cma).*Gm);
GmCma_output = GmCma + dt*((k9a.*Gm + k7.*Gc).*Cma - k9b.*GmCma);
Ii_output = Ii + dt*(k6*Ia - k4.*Cma.*Ii);
Ia_output = Ia + dt*(k4.*Cma.*Ii - k6*Ia);

% Diffusion using Fourier approach
Cci_output = diffusion_step(dt, Cci_output, D_Cci, Flk);
Cmi_output = diffusion_step(dt, Cmi_output, D_Cmi, Flk);
Cma_output = diffusion_step(dt, Cma_output, D_Cma, Flk);
Gc_output = diffusion_step(dt, Gc_output, D_Gc, Flk);
Gm_output = diffusion_step(dt, Gm_output, D_Gm, Flk);
GmCma_output = diffusion_step(dt, GmCma_output, D_GmCma, Flk);
Ii_output = diffusion_step(dt, Ii_output, D_Ii, Flk);
Ia_output = diffusion_step(dt, Ia_output, D_Ia, Flk);
end