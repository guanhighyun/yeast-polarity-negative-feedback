load('../FigureData/Figure5B.mat');
% scaling factor to convert zepto mole (10^(-21) mol) to molecules
SF = 10^(-21)*6.02*10^(23);
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
figure;
plot(k1/SF,sum_uniform_IC,'r','LineWidth',3); hold on; plot(k1/SF,sum_prepolarized_IC,'r','linewidth',3);
ylabel('Total active Cdc42 (molecules)')
xlabel('$k1 \left( \frac{\mu m^2}{s} \right)$', 'Interpreter','latex');
set(gca,'fontsize',18)
box on; axis square;
