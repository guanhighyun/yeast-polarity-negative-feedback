function [Ua, Ui, X] = Euler_step(dt,Ua,Ui,X,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,n,h,noise,Flk)
% Rate of all species
rate = (k0 + (k1+s+noise).*Ua.^2).*Ui - (k2+k3*X).*Ua;
rate_X = k4*Ua - k5*X;

% Forward Euler method
X = X + dt*rate_X;
Ua = Ua + dt*(rate);
Ui = Ui + dt*(-rate);

% Diffusion using Fourier method
Ua_ff   = fft(Ua);
Ua_nn = Ua_ff.* exp(-dt*Flk'*D_Ua); 
Ua = real(ifft(Ua_nn));

Ui_ff   = fft(Ui);
Ui_nn = Ui_ff.* exp(-dt*Flk'*D_Ui); 
Ui = real(ifft(Ui_nn));

X_ff   = fft(X);
X_nn = X_ff.* exp(-dt*Flk'*D_Ui); 
X = real(ifft(X_nn));
end