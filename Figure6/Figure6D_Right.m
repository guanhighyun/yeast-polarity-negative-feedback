% Plot Cdc42-GTP distributions
load('FigureData/Figure4D_coordinates.mat')
figure('position',[300 300 1200 400]); L = 8.8623;
tiledlayout(1,3,'tilespacing','compact');
nexttile;
plot(x_1, y_1, 'k.');
set(gca,'linewidth',3);
xticks([]); yticks([]); xlim([0,L]); ylim([0,L])

nexttile;
plot(x_2, y_2, 'k.');
set(gca,'linewidth',3)
xticks([]); yticks([]); xlim([0,L]); ylim([0,L])

nexttile;
plot(x_3, y_3, 'k.');
set(gca,'linewidth',3)
xticks([]); yticks([]); xlim([0,L]); ylim([0,L])

