load('../FigureData/Figure5B.mat');

% scaling factor to convert 1/5 molecules to 1 molecule
SF = 5;
% Calculate the total number of Cma
% Uniform initial condition
sum_uniform_IC = sum(Cma_vals_1,2)*x_len/x_num*SF;
% Prepolarized initial condition
sum_prepolarized_IC = sum(Cma_vals_2,2)*x_len/x_num*SF;
% Eliminate error point
idx=find(k1==0.74|k1==0.68);
sum_prepolarized_IC(idx) = [];
sum_uniform_IC(idx) = [];
k1(idx)=[];

sum_prepolarized_IC(k1<=0.49) = nan;
sum_uniform_IC(k1>=0.67) = nan;
figure('position',[300 300 500 500])
plot(k1/SF,sum_uniform_IC,'r','LineWidth',3); hold on; plot(k1/SF,sum_prepolarized_IC,'r','linewidth',3);
ylabel('Total active Cdc42 (molecules)')
xlabel('$k_1 \left( \frac{\mu m}{s} \right)$', 'Interpreter','latex');
set(gca,'fontsize',22)
box on; axis square;
