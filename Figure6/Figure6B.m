% Plot the bifurcation diagram with varied k3 for the model with GAP
load('../FigureData/Figure6B.mat')
figure('position',[300 300 600 500])
shadedErrorBar(Lambda3_pre_polarized_initial_condition,mean(K_pre_polarized_initial_condition),std(K_pre_polarized_initial_condition),'LineProps',{'linewidth',3,'color','r'},'transparent',true,'patchSaturation',0.3); hold on;
shadedErrorBar(Lambda3_uniform_initial_condition,mean(K_uniform_initial_condition),std(K_uniform_initial_condition),'LineProps',{'linewidth',3,'color','b'},'transparent',true,'patchSaturation',0.3)
xlabel('\lambda_1 (s^{-1})')
ylabel('Clustering score')
xlim([50,110]);
ylim([0,6]);
set(gca,'linewidth',3);
set(gca,'fontsize',22);
title('Bifurcation diagram with varying \lambda_1','FontSize',30)
legend('Gradually deceasing \lambda_1','Gradually increasing \lambda_1')