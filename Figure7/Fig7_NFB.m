clear;

% Random seeds control
rng(1);

% Grid number
n = 250;

% Final time point (secs)
tfinal = 2000;

% Time step (secs)
dt = 0.01;

% Number of time points
nt = tfinal/dt;

% Domain length
L = 5*pi;

% Grid size
h = L/n;

% Initial time step (secs)
dt2 = dt/10;

% Rate constants
k0 = 0; 
k1 = 0.24;
k2 = 0.14;
k3 = 5;
k4 = 0.01; 
k5 = 0.01;

% Diffusion rates (um^2/secs)
D_Ua = 0.0025;
D_Ui = 0.25;


% Fourier space setup
wx = 2 * pi / L * [0:n/2 -n/2+1:-1];
Flk = wx.^2;

% Noise setup
noise = 0; % no noise to s

% Pheromone signal setup
t1 = 800/dt;
t2 = t1+360/dt;

s = zeros(1,nt);
s(1:t1) = 0.02;
s(t1+1:t2) = 0.92;
s(t2+1:nt) = 0.01;

% Uniform initial condition
uniform = 1;
[Ua_uni_IC,X_uni_IC] = run_Euler_step(uniform,dt,dt2,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,n,h,noise,Flk,tfinal,nt,t1);

% Prepolarized initial condition
uniform = 0;
[Ua_prepo_IC,X_prepo_IC] = run_Euler_step(uniform,dt,dt2,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,n,h,noise,Flk,tfinal,nt,t1);

% Figure plotting
% Compare X
figure('Position',[0 0 1000 300]); SF = 10^(-21)*6.02*10^(23);
plot(1:tfinal,sum(X_prepo_IC)*h*SF,'linewidth',3); hold on; plot(1:tfinal,sum(X_uni_IC)*h*SF,'--','LineWidth',3)
xlabel('Time (min)'); ylabel('Total active X'); xticks((0:600:tfinal)); xticklabels((0:600:tfinal)/60);
set(gca,'fontsize',25);

% Compare active Cdc42
figure('Position',[0 0 1000 300]); SF = 10^(-21)*6.02*10^(23);
plot(1:tfinal,sum(Ua_prepo_IC)*h*SF,'linewidth',3); hold on; plot(1:tfinal,sum(Ua_uni_IC)*h*SF,'--','LineWidth',3)
xlabel('Time (min)'); ylabel('Total active Cdc42'); xticks((0:600:tfinal)); xticklabels((0:600:tfinal)/60);
set(gca,'fontsize',25);

% Plot s
figure('Position',[0 0 1000 300]);
grid off;
box on;
plot(dt*(1:nt), s, 'Color', [1 0 0], 'LineWidth', 2);
ylabel('$s \left( \frac{\mu m^2}{s \cdot mol^2} \right)$', 'Interpreter','latex');
xticks((0:600:tfinal));
xticklabels((0:600:tfinal)/60);
xlabel('Time (min)');
title('Pheromone signal', 'Interpreter','none');
set(gca,'fontsize',25);
ylim([0, 0.1+max(s)]);
yticklabels(round(yticks/SF/SF,7))
set(gca,'linewidth',2);


function [Ua_output,X_output] = run_Euler_step(uniform,dt,dt2,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,n,h,noise,Flk,tfinal,nt,t1)
U_concentration = 1; X = 0*ones(n,1);
Ua_output = nan(n,tfinal);
Ui_output = nan(n,tfinal);
X_output = nan(n,tfinal);
sigma = n*0.01;
a = 230;

sigma2 = n*0.01;
a2 = 1;

if uniform == 1
    % Uniform IC
    noise_Ua = 0*abs(randn(1,n));
    noise_Ua(noise_Ua>1)=1;
    Ua = noise_Ua';
    Ui = U_concentration*ones(n,1)-Ua;
else
    % Pre-polarized IC
    % Ua = zeros(n,1);
    Ua = (a*1/(sigma*sqrt(2*pi))*exp(-1/2*((1:n)-n/2).^2/sigma^2))';
    Ui = (U_concentration*n-sum(Ua))/n*ones(n,1);
end

% Main simulation loop
curr_t = 0;
for tidx = 1:nt
    curr_t = curr_t+dt;  % update simulation clock
    %% run the actual simulation
    lp = 0;      % Track status of simulation
    react_t = 0; % Internal reaction clock
    et = 1e-04; % Error tolerance
    curr_s = s(tidx);
    
    if tidx == t1 && uniform == 1
        noise_Ua2 = a2*(1/(sigma2*sqrt(2*pi))*exp(-1/2*((1:n)-n/2).^2/sigma2^2))';
        noise_Ua = 0.1*abs(randn(1,n));
        noise_Ua(noise_Ua>0.7)= 0.7;
        Ua = Ua + noise_Ua' + noise_Ua2;
        Ui = Ui - noise_Ua' - noise_Ua2;
    end

    while  lp ~= 2
        %dbstop if naninf

        [Uah, Uih, Xh] = Euler_step(dt2,Ua,Ui,X,k0,k1,k2,k3,k4,k5,curr_s,D_Ua,D_Ui,n,h,noise,Flk);
        [Uah2, Uih2, Xh2] = Euler_step(dt2/2,Ua,Ui,X,k0,k1,k2,k3,k4,k5,curr_s,D_Ua,D_Ui,n,h,noise,Flk);
     
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
figure('Position',[0 0 1000 300]);
imagesc(Ua_output)
SF = 10^(-21)*6.02*10^(23);
hold on;
colormap parula; set(gca,'ydir','normal')
xticks((0:600:tfinal));
xticklabels((0:600:tfinal)/60);
yticks([0,n/2,n])
yticklabels({'0','2.5\pi','5\pi'})
imagesc(Ua_output*SF*h);
colormap(hot);
set(gca,'fontsize',25);


colorbar
end

function [Ua, Ui, X] = Euler_step(dt,Ua,Ui,X,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,n,h,noise,Flk)

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