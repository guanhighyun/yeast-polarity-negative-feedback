% Plot the bifurcation diagram with varied k3 for the model with GAP
load('FigureData/Figure4B.mat')
figure('position',[300 300 400 400])
shadedErrorBar(Lambda3_pre_polarized_initial_condition,mean(K_pre_polarized_initial_condition),std(K_pre_polarized_initial_condition),'LineProps',{'linewidth',3,'color','r'},'transparent',true,'patchSaturation',0.3); hold on;
shadedErrorBar(Lambda3_uniform_initial_condition,mean(K_uniform_initial_condition),std(K_uniform_initial_condition),'LineProps',{'linewidth',3,'color','b'},'transparent',true,'patchSaturation',0.3)
xlabel('\lambda_3 (s^{-1})')
ylabel('Clustering score')
xlim([50,110]);
ylim([0,Inf]);
set(gca,'linewidth',3);
set(gca,'fontsize',22);