% Fig 4D
clear;
rng(1);

% Number of grids
n = 250;

% Final time point
tfinal = 2000;

% Time step
dt = 0.01;

% Number of time points
nt = tfinal/dt;

% Domain length
L = 5*pi;

% Grid size
h = L/n;

% Initial time step (for checking consistency)
dt2 = dt/10;

% Reaction rates
k0 = 0; 
k1 = 0.24;
k2 = 0.14;
k3 = 5;
k4 = 0.01; 
k5 = 0.01;

% Diffusion rates
D_Ua = 0.0025;
D_Ui = 0.25; % Same as D_X

% Fourier space setup
wx = 2 * pi / L * [0:n/2 -n/2+1:-1];
Flk = wx.^2;

% Noise setup
noise = 0; % No noise to reaction rates were set.

% Pheromone signal setup
t1 = 1000/dt; % Time at which we reduced pheromone signal
s = zeros(1,nt);
s(1:t1) = 0.92;
s(t1+1:nt) = 0.01;

% PDE with prepolarized initial condition and periodic boundary.
[Ua_prepo_IC,X_prepo_IC] = run_Euler_step(dt,dt2,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,n,h,noise,Flk,tfinal,nt);

% Figure plotting
% Active X
figure('Position',[0 0 1000 300]); SF = 10^(-21)*6.02*10^(23);
plot(1:tfinal,sum(X_prepo_IC)*h*SF,'r','linewidth',3); 
xlabel('Time (min)'); title('Total X'); xticks((0:600:tfinal)); xticklabels((0:600:tfinal)/60);
set(gca,'fontsize',25); ylabel('X (molecules)'); set(gca,'linewidth',2);

% Active Cdc42, not used in the final figure version
figure('Position',[0 0 1000 300]); SF = 10^(-21)*6.02*10^(23);
plot(1:tfinal,sum(Ua_prepo_IC)*h*SF,'r','linewidth',3); 
xlabel('Time (min)'); ylabel('Total active Cdc42'); xticks((0:600:tfinal)); xticklabels((0:600:tfinal)/60);
set(gca,'fontsize',25); title('Total active Cdc42')
set(gca,'linewidth',2);

% The pheromone signal s
figure('Position',[0 0 1000 300]);
grid off;
box on;
plot(dt*(1:nt), s, 'Color', [1 0 0], 'LineWidth', 3);
ylabel('$s \left( \frac{\mu m^2}{s \cdot mol^2} \right)$', 'Interpreter','latex');
xticks((0:600:tfinal));
xticklabels((0:600:tfinal)/60);
xlabel('Time (min)');
title('Pheromone signal', 'Interpreter','none');
set(gca,'fontsize',25);
ylim([0, 0.1+max(s)]);
yticklabels(round(yticks/SF/SF,7))
set(gca,'linewidth',2);