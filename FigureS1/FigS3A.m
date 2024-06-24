clear;

n = 250;
tfinal = 10000;
dt = 0.01;

nt = tfinal/dt;
L = 5*pi;
h = L/n;
dt2 = dt/10;

k0 = 0; % might need to reduce
k1 = 0.24;
k2 = 0.14;
k3 = 5;
k4 = 0.01; %0.2
k5 = 0.01; %1.2e-03;

D_Ua = 0.0025;
D_Ui = 0.25;

% Fourier space setup
wx = 2 * pi / L * [0:n/2 -n/2+1:-1];
Flk = wx.^2;

X = 0;

% signal
s_vals = 1.92:-0.02:0; 

% Noise setup
noise = 0;

i = 0;
uniform = 1;
for s = s_vals
    i = i+1;
    sum_Ua_uni(i) = run_Euler_step(uniform,dt,dt2,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,n,h,noise,Flk,tfinal,nt);
end

i = 0;
uniform = 0;
for s = s_vals
    i = i+1;
    sum_Ua_prepo(i) = run_Euler_step(uniform,dt,dt2,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,n,h,noise,Flk,tfinal,nt);
end
save('save_IFF_bifurcation.mat')

sum_Ua_uni(1:find(sum_Ua_uni>=30,1,'last')) = nan;
sum_Ua_prepo(sum_Ua_prepo<30) = nan;

% Scaling factor to molecules
SF = 10^(-21)*6.02*10^(23);

figure; plot(s_vals/SF/SF,sum_Ua_uni*SF*h,'r','LineWidth',4); hold on
plot(s_vals/SF/SF,sum_Ua_prepo*SF*h,'r','LineWidth',4)
ylabel('Total C_{GTP} (molecules)'); xlabel('$s(\frac{\mum^2}{s \cdot mol^2})$','Interpreter','latex')
set(gca,'fontsize',20); xlim([0,Inf]); ylim([0,Inf])
axis square

function Ua_output = run_Euler_step(uniform,dt,dt2,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,n,h,noise,Flk,tfinal,nt)
U_concentration = 1; X = 0*ones(n,1);
if uniform == 1
    % Uniform IC
    noise_Ua = 0.1*abs(randn(1,n)); % 0.3x or 1x produce same results
    noise_Ua(noise_Ua>0.7)=0.7;
    Ua = noise_Ua';
    Ui = U_concentration*ones(n,1)-Ua;
else
    % Pre-polarized IC
    % Ua = zeros(n,1);
    sigma = n*0.01;
    a =230;
    Ua = (a*1/(sigma*sqrt(2*pi))*exp(-1/2*((1:n)-n/2).^2/sigma^2))';
    Ui = (U_concentration*n-sum(Ua))/n*ones(n,1);
end


Ua_output = nan(n,tfinal);
Ui_output = nan(n,tfinal);
X_output = nan(n,tfinal);

% Main simulation loop
curr_t = 0;
for tidx = 1:nt
    curr_t = curr_t+dt;  % update simulation clock
    %% run the actual simulation
    lp = 0;      % Track status of simulation
    react_t = 0; % Internal reaction clock
    et = 1e-04; % Error tolerance

    while  lp ~= 2
        %dbstop if naninf

        [Uah, Uih, Xh] = Euler_step(dt2,Ua,Ui,X,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,n,h,noise,Flk);
        [Uah2, Uih2, Xh2] = Euler_step(dt2/2,Ua,Ui,X,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,n,h,noise,Flk);
     
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
Ua_output = sum(Ua);
end

function [Ua, Ui, X] = Euler_step(dt,Ua,Ui,X,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,n,h,noise,Flk)

rate = (k0 + (k1+s+noise).*Ua.^2).*Ui - (k2+k3*X).*Ua;
rate_X = k4*s - k5*X;

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