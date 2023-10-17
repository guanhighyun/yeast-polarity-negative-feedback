% Plot the bifurcation diagram with varied k2b for the model without GAP
load('FigureData/Figure4A.mat');
figure('position',[300 300 400 400])
shadedErrorBar(k2b(k2b<=0.1),nanmean(K_pre_polarized_initial_condition(:,k2b<=0.1)),nanstd(K_pre_polarized_initial_condition(:,k2b<=0.1)),'LineProps',{'linewidth',3,'color','R'},'transparent',true,'patchSaturation',0.3); hold on;
shadedErrorBar(k2b(k2b<=0.1),nanmean(K_uniform_initial_condition(:,k2b<=0.1)),nanstd(K_uniform_initial_condition(:,k2b<=0.1)),'LineProps',{'linewidth',3,'color','B'},'transparent',true,'patchSaturation',0.3); 
xlabel('\lambda_3 (s^{-1})')
ylabel('Clustering score')
xlim([0.015,0.065]);
ylim([0,5]);
set(gca,'linewidth',3);
set(gca,'fontsize',22);

figure('position',[300 300 400 400])
shadedErrorBar(k2b(k2b>0.1),nanmean(K_pre_polarized_initial_condition(:,k2b>0.1)),nanstd(K_pre_polarized_initial_condition(:,k2b>0.1)),'LineProps',{'linewidth',3,'color','R'},'transparent',true,'patchSaturation',0.3); hold on;
shadedErrorBar(k2b(k2b>0.1),nanmean(K_uniform_initial_condition(:,k2b>0.1)),nanstd(K_uniform_initial_condition(:,k2b>0.1)),'LineProps',{'linewidth',3,'color','B'},'transparent',true,'patchSaturation',0.3); 
xlabel('\lambda_3 (s^{-1})')
ylabel('Clustering score')
xlim([6.2,8.2]);
ylim([0,5]);
set(gca,'linewidth',3);
set(gca,'fontsize',22);
