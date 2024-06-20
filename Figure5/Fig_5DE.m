clear; clc;

rng('shuffle');

save_data = false;

r_yc = 2.5;

m = 0;
x_len = r_yc*pi;
x_num = 250;
t_max = 50000;
t_num = 5000;
x = linspace(0,x_len,x_num);

k1a = 10; 
k1b = 40; 
k2a = 0.0382; 
k2b = 1;
k3_1 = 10; 
k3_2 = 0.38;
t1 = 5000;
t = [linspace(0,t_max,t_num)];

k4a = 0.059;
k4b = 34.5262;
k5a = 36;
k5b = 13;
k7 = 1.8474;

k6a = 0.001;
k6b = 0.8;
k6c = 0.01;

C_tc = 63.662;
G_tc = 6.3662;
I_tc = 5;

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
 

C_all = C_tc*x_len;
ex_C_start = x_len/2;
ex_C_end = x_len;

%initial conditions 1
Cmi_init_min = 43.2719;
%Cmi_init_max = 0*C_all/(ex_C_end-ex_C_start);
Cmi_init_max = 43.2719;
Cma_init_min = 4.8459;
%Cma_init_max = 0*C_all/(ex_C_end-ex_C_start);
Cma_init_max = 4.8459;
Cci_init = (C_all - Cmi_init_max*(ex_C_end-ex_C_start) - Cmi_init_min*(x_len-ex_C_end+ex_C_start) ...
                  - Cma_init_max*(ex_C_end-ex_C_start) - Cma_init_min*(x_len-ex_C_end+ex_C_start))/x_len;

%initial conditions 2
Cmi_init_min_2 = 43.2719;
%Cmi_init_max = 0*C_all/(ex_C_end-ex_C_start);
Cmi_init_max_2 = 43.2719;
Cma_init_min_2 = 0;
%Cma_init_max = 0*C_all/(ex_C_end-ex_C_start);
Cma_init_max_2 = 30;
Cci_init_2 = (C_all - Cmi_init_max_2*(ex_C_end-ex_C_start) - Cmi_init_min_2*(x_len-ex_C_end+ex_C_start) ...
                  - Cma_init_max_2*(ex_C_end-ex_C_start) - Cma_init_min_2*(x_len-ex_C_end+ex_C_start))/x_len;

Gc_init = 3.6650;
Gm_init = 1.7365;
GmCma_init = G_tc - Gc_init - Gm_init;

Ii_init = 4.7689;
Ia_init = I_tc - Ii_init;


opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

polarize = 0;
P1 = [x_len, t_max, k1a, k1b, k2a, k2b, k3_1 k4a, k4b, k5a, k5b, k7, k6a, k6b, k6c, ...
         D_Cci, D_Cmi, D_Cma, D_Gc, D_Gm, D_GmCma, D_Ii, D_Ia, ...
         ex_C_start, ex_C_end, Cmi_init_min_2, Cmi_init_max_2, Cma_init_min_2, ...
         Cma_init_max_2, Cci_init_2, Gc_init, Gm_init, GmCma_init, Ii_init, Ia_init,polarize,k3_2,t1];
sol1 = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t,opts,P1);

polarize = 1;
P2 = [x_len, t_max, k1a, k1b, k2a, k2b, k3_1 k4a, k4b, k5a, k5b, k7, k6a, k6b, k6c, ...
         D_Cci, D_Cmi, D_Cma, D_Gc, D_Gm, D_GmCma, D_Ii, D_Ia, ...
         ex_C_start, ex_C_end, Cmi_init_min_2, Cmi_init_max_2, Cma_init_min_2, ...
         Cma_init_max_2, Cci_init_2, Gc_init, Gm_init, GmCma_init, Ii_init, Ia_init,polarize,k3_2,t1];
sol2 = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t,opts,P2);

SF = 10^(-21)*6.02*10^(23);
k3_eff = zeros(1, t_num);
for i = 1:length(t)
    ti = t(i);
    if (ti < t1)
        k3_eff(i) = k3_1;
    end
    if (ti >= t1)
        k3_eff(i) = k3_2;
    end
end
figure('Position',[50 30 1000 900]); 
subplot(2,1,1)
hold on;
grid off;
box on;
plot(t, k3_eff/SF/SF, 'Color', [1 0 0], 'LineWidth', 3);
delta_t = 6000; 

xticks(0:delta_t:t_max);
xticklabels((0:delta_t:t_max)/60);
xlabel('Time (min)', 'Interpreter', 'Latex');
ylabel('$k_1 \left( \frac{\mu m^2}{s \cdot mol^2} \right)$', 'Interpreter','latex');
set(gca,'fontsize', 28); xlim([0,Inf])
set(gca,'linewidth', 3);

ax = subplot(2,1,2);
Ua_vals = sol1(:,:,3)'+sol1(:,:,6)';
%Ua_vals = [Ua_vals; flip(Ua_vals,1)];
imagesc((sol1(:,:,3)'+sol1(:,:,6)').*x_len/x_num*SF); set(gca,'ydir','normal')
colormap(hot); 

xlabel('Time (min)');
dt = t_max/t_num; 
xticks((0:delta_t:t_max)/dt);
xticklabels(0:100:t_max/60);

ylabel({'Distance (\mum)'});
yticks([0,x_num/2,x_num]);
yticklabels({'0','2.5\pi','5\pi'}); 

title('Active Cdc42','Color', [0,0,0]);
set(ax,'fontsize',28);

k3_eff = zeros(1, t_num);
for i = 1:length(t)
    ti = t(i);
    if (ti < t1)
        k3_eff(i) = k3_1;
    elseif (ti >= t1) && (ti <= t_max-10000)
        k3_eff(i) = k3_1 + (k3_2 - k3_1)/(t_max-t1-10000)*(ti-t1);
    else
        k3_eff(i) = k3_2;
    end
end
figure('Position',[50 30 1000 900]); 
subplot(2,1,1)
hold on;
grid off;
box on;
plot(t, k3_eff/SF/SF, 'Color', [1 0 0], 'LineWidth', 3);
delta_t = 6000; 

xticks(0:delta_t:t_max);
xticklabels((0:delta_t:t_max)/60);
xlabel('Time (min)', 'Interpreter', 'Latex');
ylabel('$k_1 \left( \frac{\mu m^2}{s \cdot mol^2} \right)$', 'Interpreter','latex');
set(gca,'fontsize', 28); xlim([0,Inf])
set(gca,'linewidth', 3);

ax = subplot(2,1,2);
Ua_vals = sol2(:,:,3)'+sol2(:,:,6)';
%Ua_vals = [Ua_vals; flip(Ua_vals,1)];
imagesc(Ua_vals.*x_len/x_num*SF); set(gca,'ydir','normal')
caxis([0,23000])
colormap(hot);

xlabel('Time (min)');
dt = t_max/t_num; 
xticks((0:delta_t:t_max)/dt);
xticklabels(0:100:t_max/60);

ylabel({'Distance (\mum)'});
yticks([0,x_num/2,x_num]);
yticklabels({'0','2.5\pi','5\pi'}); 

title('Active Cdc42','Color', [0,0,0]);
set(ax,'fontsize',28);
cb = colorbar('Position', [0.89 0.1100 0.01 0.1967], ...
    'Ticks', c_ticks, ...
    'TickLabels', c_ticklabels);
set(sp3, 'InnerPosition', [0.1450 0.1100 0.74 0.1964]);



% --------------------------------------------------------------
function [c,f,s] = pdex4pde(x,t,u,DuDx,P)

% P1 = [x_len, t_max, k1a, k1b, k2a, k2b, k3(j), k4a, k4b, k5a, k5b, k7, k6a, k6b, k6c, ...
%              D_Cci, D_Cmi, D_Cma, D_Gc, D_Gm, D_GmCma, D_Ii, D_Ia, ...
%              ex_C_start, ex_C_end, Cmi_init_min, Cmi_init_max, Cma_init_min, ...
%              Cma_init_max, Cci_init, Gc_init, Gm_init, GmCma_init, Ii_init, Ia_init];

x_len = P(1);
t_max = P(2);
k1a = P(3);
k1b = P(4);
k2a = P(5);
k2b = P(6);
k3_1 = P(7);
k4a = P(8);
k4b = P(9);
k5a = P(10);
k5b = P(11);
k7 = P(12);
k6a = P(13);
k6b = P(14);
k6c = P(15);
D_Cci = P(16);
D_Cmi = P(17);
D_Cma = P(18);
D_Gc = P(19);
D_Gm = P(20);
D_GmCma = P(21);
D_Ii = P(22);
D_Ia = P(23);
k3_2 = P(37);
t1 = P(38);
polarize = P(36);

c = [1; 1; 1; 1; 1; 1; 1; 1];
f = [D_Cci; D_Cmi; D_Cma; D_Gc; D_Gm; D_GmCma; D_Ii; D_Ia] .* DuDx; 

Cci = u(1);
Cmi = u(2);
Cma = u(3);
Gc = u(4);
Gm = u(5);
GmCma = u(6);
Ii = u(7);
Ia = u(8);

% r = 0.000001*(x_len-x);
% k5a = k5a + r;
% k2a = k2a + r;

if polarize == 1
    if (t < t1)
        k3 = k3_1;
    elseif (t >= t1) && (t <= t_max-10000)
        k3 = k3_1 + (k3_2 - k3_1)/(t_max-t1-10000)*(t-t1);
    else
        k3 = k3_2;
    end
else
    if (t < t1)
        k3 = k3_1;
    else
        k3 = k3_2;
    end
end


F1 = k5b*Cmi - k5a*Cci;
F2 = (k2b + k6b*Ia)*Cma + k5a*Cci - (k2a*Gm + k3*GmCma + k5b)*Cmi;
F3 = (k2a*Gm + k3*GmCma)*Cmi + k4b*GmCma - (k2b + k4a*Gm + k7*Gc + k6b*Ia)*Cma;
F4 = k1b*Gm - (k1a + k7*Cma)*Gc;
F5 = k1a*Gc + k4b*GmCma - (k1b + k4a*Cma)*Gm;
F6 = (k4a*Gm + k7*Gc)*Cma - k4b*GmCma;
F7 = k6c*Ia - k6a*Cma*Ii;
F8 = k6a*Cma*Ii - k6c*Ia;

s = [F1; F2; F3; F4; F5; F6; F7; F8];

end

% --------------------------------------------------------------
function u0 = pdex4ic(x,P)

ex_C_start = P(24);
ex_C_end = P(25);
Cmi_init_min = P(26);
Cmi_init_max = P(27);
Cma_init_min = P(28);
Cma_init_max = P(29);
Cci_init = P(30);
Gc_init = P(31);
Gm_init = P(32);
GmCma_init = P(33);
Ii_init = P(34);
Ia_init = P(35);


if (x >= ex_C_start) && ( x <= ex_C_end) && (ex_C_start ~= ex_C_end)
    Cmi_init = Cmi_init_max;
    Cma_init = Cma_init_max;
else
    Cmi_init = Cmi_init_min;
    Cma_init = Cma_init_min;
end


r_G = 0;
r_C = 0;

u0 = [Cci_init-r_C; Cmi_init; Cma_init+r_C; Gc_init-r_G; Gm_init; GmCma_init+r_G; Ii_init; Ia_init];

end

% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,t,P)
pl = [0;0;0;0;0;0;0;0];
ql = [1;1;1;1;1;1;1;1];
pr = [0;0;0;0;0;0;0;0];
qr = [1;1;1;1;1;1;1;1];
end
