function [Cci_output, Cmi_output, Cma_output, Gc_output, Gm_output, GmCma_output, Ii_output, Ia_output] = Euler_step(dt, Cci, Cmi, Cma, Gc, Gm, GmCma, Ii, Ia, D_Cci, D_Cmi, D_Cma, D_Gc, D_Gm, D_GmCma, D_Ii, D_Ia, Flk,...
    k1a,k1b,k2a,k2b,k3,k4a,k4b,k5a,k5b,k6a,k6b,k6c,k7)
% Reaction rate calculation

Cci_output = Cci + dt*(k5b*Cmi - k5a*Cci);
Cmi_output = Cmi + dt*((k2b + k6b.*Ia).*Cma + k5a*Cci - (k2a.*Gm + k3.*GmCma + k5b).*Cmi);
Cma_output = Cma + dt*((k2a.*Gm + k3.*GmCma).*Cmi + k4b*GmCma - (k2b + k4a.*Gm + k7.*Gc + k6b.*Ia).*Cma);
Gc_output = Gc + dt*(k1b*Gm - (k1a + k7.*Cma).*Gc);
Gm_output = Gm + dt*(k1a*Gc + k4b*GmCma - (k1b + k4a.*Cma).*Gm);
GmCma_output = GmCma + dt*((k4a.*Gm + k7.*Gc).*Cma - k4b.*GmCma);
Ii_output = Ii + dt*(k6c*Ia - k6a.*Cma.*Ii);
Ia_output = Ia + dt*(k6a.*Cma.*Ii - k6c*Ia);

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