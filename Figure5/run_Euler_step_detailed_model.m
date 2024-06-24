function [Cma_output, Ia_output] = run_Euler_step_detailed_model(dt2, Cci, Cmi, Cma, Gc, Gm, GmCma, Ii, Ia, D_Cci, D_Cmi, D_Cma, D_Gc, D_Gm, D_GmCma, D_Ii, D_Ia, Flk,...
            k8a,k8b,k10,k2,k1_eff,k9a,k9b,k5a,k5b,k4,k3,k6,k7,n,tfinal,nt,dt)

Cci_output = nan(n,tfinal);
Cmi_output = nan(n,tfinal);
Cma_output = nan(n,tfinal);
Gc_output = nan(n,tfinal);
Gm_output = nan(n,tfinal);
GmCma_output = nan(n,tfinal);
Ii_output = nan(n,tfinal);
Ia_output = nan(n,tfinal);

% Main simulation loop
curr_t = 0;
for tidx = 1:nt
    curr_t = curr_t+dt;  % update simulation clock
    %% run the actual simulation
    lp = 0;      % Track status of simulation
    react_t = 0; % Internal reaction clock
    et = 1e-04; % Error tolerance
    k1 = k1_eff(tidx);
    while  lp ~= 2
        %dbstop if naninf

        [Ccih, Cmih, Cmah, Gch, Gmh, GmCmah, Iih, Iah] = Euler_step(dt2, Cci, Cmi, Cma, Gc, Gm, GmCma, Ii, Ia, D_Cci, D_Cmi, D_Cma, D_Gc, D_Gm, D_GmCma, D_Ii, D_Ia, Flk,...
            k8a,k8b,k10,k2,k1,k9a,k9b,k5a,k5b,k4,k3,k6,k7);
        [Ccih2, Cmih2, Cmah2, Gch2, Gmh2, GmCmah2, Iih2, Iah2] = Euler_step(dt2/2, Cci, Cmi, Cma, Gc, Gm, GmCma, Ii, Ia, D_Cci, D_Cmi, D_Cma, D_Gc, D_Gm, D_GmCma, D_Ii, D_Ia, Flk,...
            k8a,k8b,k10,k2,k1,k9a,k9b,k5a,k5b,k4,k3,k6,k7);
     
        % error checking
        Ccie = max(max(abs(Ccih   - Ccih2   )));
        Cmie    = max(max(abs(Cmih   - Cmih2   )));
        Cmae  = max(max(abs(Cmah - Cmah2 )));
        Gce  = max(max(abs(Gch - Gch2 )));
        Gme  = max(max(abs(Gmh - Gmh2 )));
        GmCmae  = max(max(abs(GmCmah - GmCmah2 )));
        Iie  = max(max(abs(Iih - Iih2 )));
        Iae  = max(max(abs(Iah - Iah2 )));
        
        r = max([Ccie, Cmie, Cmae, Gme, GmCmae, Iie, Iae]);
        
        % doing the step again with smaller stepsize
        if r > et && lp == 0
            dt2 = dt2/1.5;     
        else
            Cci  = Ccih;
            Cmi  = Cmih;
            Cma = Cmah;
            Gc = Gch;
            Gm = Gmh;
            GmCma = GmCmah;
            Ii = Iih;
            Ia = Iah;

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
        Cci_output(:,int64(tidx*dt)) = Cci;
        Cmi_output(:,int64(tidx*dt)) = Cmi;
        Cma_output(:,int64(tidx*dt)) = Cma;
        Gc_output(:,int64(tidx*dt)) = Gc;
        Gm_output(:,int64(tidx*dt)) = Gm;
        GmCma_output(:,int64(tidx*dt)) = GmCma;
        Ii_output(:,int64(tidx*dt)) = Ii;
        Ia_output(:,int64(tidx*dt)) = Ia;
    end
end

figure; imagesc(Cma_output)
figure; plot(sum(Cma_output))
figure; plot(sum(Ia_output))
end
