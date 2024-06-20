function [Ua_output,X_output] = run_Euler_step(dt,dt2,k0,k1,k2,k3,k4,k5,s,D_Ua,D_Ui,n,h,noise,Flk,tfinal,nt)
U_concentration = 1; 
Ua_output = nan(n,tfinal);
Ui_output = nan(n,tfinal);
X_output = nan(n,tfinal);
sigma = n*0.01;
a = 230;

% Initial conditions
Ua = (a*1/(sigma*sqrt(2*pi))*exp(-1/2*((1:n)-n/2).^2/sigma^2))';
Ui = (U_concentration*n-sum(Ua))/n*ones(n,1);
X = zeros(n,1);

% Main simulation loop. Part of the code is adapted from the code of Mike Pablo, https://github.com/mikepab/yeastpolariz_particlesim
% Pablo M, Ramirez SA, Elston TC. PLoS Comput Biol. 2018.

curr_t = 0;
for tidx = 1:nt
    % Record simulation time
    curr_t = curr_t+dt; 

    % Internal reaction clock
    react_t = 0; 

    % Check if the internal reaction clock is within dt
    lp = 0;  

    % Error tolerance
    et = 1e-04; 

    % Current pheromone signal magnitude
    curr_s = s(tidx);

    while  lp ~= 2
        % Run forward-Euler simulation
        [Uah, Uih, Xh] = Euler_step(dt2,Ua,Ui,X,k0,k1,k2,k3,k4,k5,curr_s,D_Ua,D_Ui,n,h,noise,Flk);
        [Uah2, Uih2, Xh2] = Euler_step(dt2/2,Ua,Ui,X,k0,k1,k2,k3,k4,k5,curr_s,D_Ua,D_Ui,n,h,noise,Flk);
     
        % Error
        Uae    = max(max(abs(Uah   - Uah2   )));
        Uie  = max(max(abs(Uih - Uih2 )));
        Xe  = max(max(abs(Xh - Xh2 )));
        
        % Maximum error of all
        r = max([Uae,Uie,Xe]);
        
        % Conduct the simulation with smaller time step
        if r > et && lp == 0
            dt2 = dt2/1.5;     
        else
            Ua  = Uah;
            Ui  = Uih;
            X  = Xh2;
           react_t = react_t + dt2;

           % New time step
            dt2 = 1.5*dt2;

            if lp == 1
                lp = 2;
            end
        end   

        % Check if the simulation is within an individual time step
        if react_t + dt2 > dt && lp ~= 2
            dt2 = dt - react_t;
            lp = 1; 
        end            
    end

    % Record the concentration of the species
    if mod(round(tidx*dt,5),1) == 0
        Ua_output(:,int64(tidx*dt)) = Ua;
        Ui_output(:,int64(tidx*dt)) = Ui;
        X_output(:,int64(tidx*dt)) = X;
    end
end

% Plot spatiotemporal species amount
figure('Position',[0 0 1000 300]);
SF = 10^(-21)*6.02*10^(23);
imagesc(Ua_output*SF*h)
colormap parula; set(gca,'ydir','normal')
xticks((0:600:tfinal));
xticklabels((0:600:tfinal)/60);
yticks([0,n/2,n])
yticklabels({'0','2.5\pi','5\pi'})
colormap(hot);
set(gca,'fontsize',25);
xlabel('Time (min)');
ylabel('Distance (\mum)')
title('Active Cdc42')
caxis([0,1950])
colorbar
end