% This code is adapted from Pablo M, Ramirez SA, Elston TC (2018) Particle-based simulations of polarity establishment reveal stochastic promotion of Turing pattern formation. PLOS Computational Biology 14(3): e1006016. https://doi.org/10.1371/journal.pcbi.1006016
clear;

n = 250; % Number of grids
tfinal = 10000; % Final time point
dt = 0.01; % Time step

nt = tfinal/dt; % Number of time steps
L = 5*pi; % Domain length
h = L/n; % spatial discretization
dt2 = dt/10; % Initial time step

% Scaling factor to convert zeptomole to molecules
SF = 10^(-21)*6.02*10^(23);

% Model parameters
% Reaction rates, not scaled by SF. 
k0 = 0; 
k1 = 0.24; % um^2/(s*zeptomole^2)
k2 = 0.14; % /s
k3 = 5; % um/(s*zeptomole)
k4 = 0.01; % /s
k5 = 0.01; % /s

% Diffusion coefficient
D_Ua = 0.0025; % um^2/s
D_Ui = 0.25; %um^2/s

% Fourier space setup
wx = 2 * pi / L * [0:n/2 -n/2+1:-1];
Flk = wx.^2;

% No negative regulator for the phase space scan.
X = 0;

% signal
s_vals = 0.56:-0.02:0; 

% Noise setup. No noise.
noise = 0;

i = 0;
uniform = 1;
for s = s_vals
    i = i+1;
    sum_Ua_uni(i) = run_Euler_step(uniform,dt,dt2,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,noise,Flk,tfinal,nt,n);
end

i = 0;
uniform = 0; % Uniform initial condition
sum_Ua_prepo = nan(numel(s_vals),1);
for s = s_vals
    i = i+1;
    sum_Ua_prepo(i) = run_Euler_step(uniform,dt,dt2,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,noise,Flk,tfinal,nt,n);
end
% Plot the figures
sum_Ua_uni(sum_Ua_uni>=30) = nan; % Eliminate repeated regions
sum_Ua_prepo(sum_Ua_prepo<30) = nan; % Eliminate repeated regions

% Scale the signal by 1/(SF*SF), and scale the output species by SF*h to
% convert from zeptomole to molecules
figure; plot(s_vals/SF/SF,sum_Ua_uni*SF*h,'r','LineWidth',4); hold on
plot(s_vals/SF/SF,sum_Ua_prepo*SF*h,'r','LineWidth',4)
ylabel('Total C_{GTP} (molecules)')
set(gca,'fontsize',20);
axis square
xlim([0,Inf])

function Ua_output = run_Euler_step(uniform,dt,dt2,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,noise,Flk,tfinal,nt,n)
U_concentration = 1; X = 0*ones(n,1);
if uniform == 1
    % Uniform initial condition. Random noise.
    noise_Ua = 0.1*abs(randn(1,n)); % 0.3x or 1x produce same results
    noise_Ua(noise_Ua>0.7)=0.7;
    Ua = noise_Ua';
    Ui = U_concentration*ones(n,1)-Ua;
else
    % Pre-polarized initial condition.
    sigma = n*0.01;
    a =230;
    % Gaussian peak of active Cdc42
    Ua = (a*1/(sigma*sqrt(2*pi))*exp(-1/2*((1:n)-n/2).^2/sigma^2))';
    Ui = (U_concentration*n-sum(Ua))/n*ones(n,1);
end

% Set up vectors to store outputs
Ua_output = nan(n,tfinal);
Ui_output = nan(n,tfinal);
X_output = nan(n,tfinal);

% Main simulation loop. Adapted from code written by Dr. Michael Pablo.
curr_t = 0;
for tidx = 1:nt
    curr_t = curr_t+dt;  % update simulation clock
    %% run the actual simulation
    lp = 0;      % Track status of simulation
    react_t = 0; % Internal reaction clock
    et = 1e-04; % Error tolerance

    while  lp ~= 2

        [Uah, Uih, Xh] = Euler_step(dt2,Ua,Ui,X,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,noise,Flk);
        [Uah2, Uih2, Xh2] = Euler_step(dt2/2,Ua,Ui,X,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,noise,Flk);
     
        % error checking
        Uae    = max(max(abs(Uah   - Uah2   )));
        Uie  = max(max(abs(Uih - Uih2 )));
        Xe  = max(max(abs(Xh - Xh2 )));
        
        r = max([Uae,Uie,Xe]);
        
        % doing the step again with smaller stepsize
        if r > et && lp == 0
            dt2 = dt2/1.5;     
        else
            Ua  = Uah;
            Ui  = Uih;
            X  = Xh2;
            % adding the step on 
           react_t = react_t + dt2;
            % setting new stepsize for next step
            dt2 = 1.5*dt2;

            if lp == 1
                lp = 2;
            end
        end   

        % Making sure we stay within the dt interval
        if react_t + dt2 > dt && lp ~= 2
            dt2 = dt - react_t;
            lp = 1; 
        end
    end

    if mod(round(tidx*dt,5),1) == 0
    Ua_output(:,int64(tidx*dt)) = Ua;
    Ui_output(:,int64(tidx*dt)) = Ui;
    X_output(:,int64(tidx*dt)) = X;
    end
end
% Ouput the total number of active Cdc42
Ua_output = sum(Ua);
end

% Perform forward Euler method
function [Ua, Ui, X] = Euler_step(dt,Ua,Ui,X,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,noise,Flk)
rate = (k0 + (k1+s+noise).*Ua.^2).*Ui - (k2+k3*X).*Ua;
rate_X = k4*Ua - k5*X;

X = X + dt*rate_X;
Ua = Ua + dt*(rate);
Ui = Ui + dt*(-rate);

% Diffusion using Fourier approach
Ua = diffusion_step(dt, Ua, D_Ua, Flk);
Ui = diffusion_step(dt, Ui, D_Ui, Flk);
X = diffusion_step(dt, X, D_Ui, Flk);
end

function result = diffusion_step(dt, data, diffusion_rate, Flk)
    data_fft = fft(data);
    data_fft = data_fft .* exp(-dt * Flk' * diffusion_rate);
    result = real(ifft(data_fft));
end