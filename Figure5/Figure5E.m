clear;

uniform = 0;
L = 5*pi;
n = 250;
tfinal = 5000;
dt = 0.01;
dt2 = dt/10;
nt = tfinal/dt;
h = L/n;

k8a = 10; 
k8b = 40; 
k10 = 0.0382;
k2 = 1;
k9a = 0.059;
k9b = 34.5262;
k5a = 36;
k5b = 13;
k7 = 1.8474;
k4 = 0.001;
k3 = 0.8;
k6 = 0.01;

C_tc = 63.662;
G_tc = 6.3662;
I_tc = 6.3662*2;

Dm = 0.0025;      % membrane diffusion coefficient
Dc = 15;          % cytoplasmic diffusion coefficient

D_Cci = Dc;
D_Cmi = Dm;
D_Cma = Dm;
D_Gc = Dc;
D_Gm = Dm;
D_GmCma = Dm;
D_Ii = Dc;
D_Ia = Dc;

% Fourier space setup
wx = 2 * pi / L * [0:n/2 -n/2+1:-1];
Flk = wx.^2;

% Preallocate output arrays
C_output = zeros(n,tfinal);
G_output = zeros(n,tfinal);
% X_output = zeros(n,tfinal);

% k3 setup
t1 = 200/dt; % Time at which we started to reduce pheromone signal

k1_1 = 10;
k1_2 = 0.6;

k1_eff = zeros(1, nt);
for i = 1:nt
    if (i < t1)
        k1_eff(i) = k1_1;
    end
    if (i >= t1)
        k1_eff(i) = k1_2;
    end
end

% Initial conditions
if uniform == 1
    noise_Ua = 0.3*abs(randn(1,n));
    noise_Ua(noise_Ua>0.7)=0.7;
    Cma = noise_Ua';
    Cci = (C_tc*n-sum(Cma))/n*ones(n,1);
else
    sigma = n*0.05;
    a = 15000;
    Cma = (a*1/(sigma*sqrt(2*pi))*exp(-1/2*((1:n)-n/2).^2/sigma^2))';
    Cci = (C_tc*n-sum(Cma))/n*ones(n,1);
end

Cmi = zeros(n,1);
Gc = G_tc*ones(n,1);
Gm = zeros(n,1);
GmCma = zeros(n,1);
Ii = I_tc*ones(n,1);
Ia = zeros(n,1);

[Cma_output, Ia_output] = run_Euler_step_detailed_model(dt2, Cci, Cmi, Cma, Gc, Gm, GmCma, Ii, Ia, D_Cci, D_Cmi, D_Cma, D_Gc, D_Gm, D_GmCma, D_Ii, D_Ia, Flk,...
            k8a,k8b,k10,k2,k1_eff,k9a,k9b,k5a,k5b,k4,k3,k6,k7,n,tfinal,nt,dt);

SF = 10^(-21)*6.02*10^(23);
figure('Position',[50 30 1000 900]); 
subplot(2,1,1)
hold on;
grid off;
box on;
plot((1:nt)*dt/60, k1_eff/SF, 'Color', [1 0 0], 'LineWidth', 3);
xlim([0,Inf])

xlabel('Time (min)', 'Interpreter', 'Latex');
ylabel('$k_1 \left( \frac{\mu m^2}{s} \right)$', 'Interpreter','latex');
set(gca,'fontsize', 28); 
set(gca,'linewidth', 3);

ax = subplot(2,1,2);
imagesc(Cma_output*h*SF); set(gca,'ydir','normal')
colorbar;
colormap(hot)
xticks(0:600:tfinal)
xticklabels(xticks/60)
xlabel('Time (min)');
ylabel({'Distance (\mum)'});
yticks([0,n/2,n]);
yticklabels({'0','2.5\pi','5\pi'}); 
set(gca,'fontsize', 28); 