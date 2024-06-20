%% Calculate time series of Ripley's K values
% We will read the sample data file which includes coordinates of all molecules.
% Lambda3 was suddenly reduced from 500 /s to 95/s.
filename = 'FigureData/Sample_file_suddenly_reduce_positive_feedback.xyz';
% Domain length
L = 8.8623;
% Maximum time frames
nframes = 801;
% The range of search radius for calculation of Ripley's K-function values
r_vec = 0.1:0.1:2.5; 

% Read molecular coordinates from each line of the .xyz file.
% nframes: total number of time points in the simulation
[t,coordinates]=read_molPos3(filename,nframes);
finalidx = find(~isnan(t),1,'last');
% Extract coordinates of active Cdc42 (in the form of Cdc42-GTP and
% the complex Bem1-GEF-Cdc42-GTP)
[active_xyz] = get_active_Cdc42_distr(coordinates,finalidx);
% Remove the first 400 seconds of calibration time
K = nan(1,finalidx-40);
lambda3 = nan(1,finalidx-40);
for i = 41:finalidx
    xyz = active_xyz{i};
    % Calculate the centroid of Cdc42-GTP cluster
    centroid = mean(xyz);
    x = xyz(:,1);
    y = xyz(:,2);
    % Move the center of the domain to the centroid of Cdc42-GTP cluster
    x(x-centroid(1)>L/2) = x(x-centroid(1)>L/2) - L;
    x(x-centroid(1)<-L/2) = x(x-centroid(1)<-L/2) + L;
    y(y-centroid(2)>L/2) = y(y-centroid(2)>L/2) - L;
    y(y-centroid(2)<-L/2) = y(y-centroid(2)<-L/2) + L;
    % Calculate Ripley's K-function values
    [H]=compute_Kr(x,y,r_vec,L);
    % Use the maximum value as the clustering score
    [maxH,~] = max(H);
    K(i-40) = maxH;
    
    if t(i) >= 600
        lambda3(i-40) = 95;
    else
        lambda3(i-40) = 500;
    end
end

%% Plot time series of Ripley's K values
% The nth time point is (n-1)*10 seconds because we stored molecular coordinates
% every other 10 seconds and the first time point is 0 seconds.
figure('position',[300 300 400 400]); hold on
time = 0:10:(finalidx-41)*10; % Remove the first 400 seconds of calibration time
time = time/60; % convert from seconds to minutes
yyaxis left
plot(time, K, 'linewidth', 4); set(gca,'fontsize',20);
xlabel('Time (min)'); ylabel('Clustering score'); xlim([0,Inf]); ylim([0,4.5]); set(gca,'linewidth',3);box off;
% This plot is the same as the black line in Figure 3C, which was generated
% with the same .xyz data file.

%% Plot k3 across time
yyaxis right
plot(time,lambda3,'linewidth', 4); ylim([0,550]); ylabel('\lambda_1(s^{-1})');