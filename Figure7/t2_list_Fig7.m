clear;
rng(1);

n = 250;
tfinal = 2000;
dt = 0.01;

nt = tfinal/dt;
L = 5*pi;
h = L/n;
dt2 = dt/10;

k0 = 0; % might need to reduce
k1 = 0.16;
k2 = 0.14;
k3 = 5;
k4 = 0.01; %0.2
k5 = 0.01; %1.2e-03;

D_Ua = 0.0025;
D_Ui = 0.25;


% Fourier space setup
wx = 2 * pi / L * [0:n/2 -n/2+1:-1];
Flk = wx.^2;

% Noise setup
noise = 0;
t1 = 800/dt;
t2 = t1+400/dt;nt-1;t1+(20*60)/dt;

t2_list = t1+(10:10:550)/dt;

i = 1:numel(t2_list);
for i = 1:numel(t2_list)
    t2 = t2_list(i);

    s = zeros(1,nt);
    s(1:t1) = 0.1;
    s(t1+1:t2) = 1;%;
    s(t2+1:nt) = 0.09;
    
    uniform = 1;
    [Ua_uni_IC(:,i),X_uni_IC(:,i)] = run_Euler_step(uniform,dt,dt2,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,n,h,noise,Flk,tfinal,nt,t1);
    
    uniform = 0;
    [Ua_prepo_IC(:,i),X_prepo_IC(:,i)] = run_Euler_step(uniform,dt,dt2,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,n,h,noise,Flk,tfinal,nt,t1);

    diff(i) = sum(sum(Ua_uni_IC(:,i))) - sum(sum(Ua_prepo_IC(:,i)));
end
diff(diff>=30)
save("t2_list_Fig_7.mat")


function [Ua_end,X_end] = run_Euler_step(uniform,dt,dt2,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,n,h,noise,Flk,tfinal,nt,t1)
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
Ua_end = Ua_output(:,end);
Ui_end = Ui_output(:,end);
X_end = X_output(:,end);
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