% Plot the bifurcation diagram with varied k2b for the model without GAP
load('FigureData/Figure8A.mat');
figure('position',[300 300 600 500])
shadedErrorBar(k2b(k2b>0.1),nanmean(K_pre_polarized_initial_condition(:,k2b>0.1)),nanstd(K_pre_polarized_initial_condition(:,k2b>0.1)),'LineProps',{'linewidth',3,'color','R'},'transparent',true,'patchSaturation',0.3); hold on;
shadedErrorBar(k2b(k2b>0.1),nanmean(K_uniform_initial_condition(:,k2b>0.1)),nanstd(K_uniform_initial_condition(:,k2b>0.1)),'LineProps',{'linewidth',3,'color','B'},'transparent',true,'patchSaturation',0.3); 
xlabel('k_2 (s^{-1})')
ylabel('Clustering score')
xlim([6.2,8.2]);
ylim([0,5]);
set(gca,'linewidth',3);
set(gca,'fontsize',22);
