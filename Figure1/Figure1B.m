% Add the folder "errorbar" to the path.
addpath('../errorbar')
% Load data from Excel files
DfU_step_down = readtable('../FigureData/DfU_step_down.xlsx');
DfU_step_up = readtable('../FigureData/DfU_step_up.xlsx');

% Create figure and tiled layout
figure;
t = tiledlayout(1, 4, 'TileSpacing', 'compact');

% Create background axis
bgAx = axes(t, 'XTick', [], 'YTick', [], 'Box', 'off');
bgAx.Layout.TileSpan = [1 2];

% First subplot (ax1)
ax1 = nexttile(1,[1,3]);
shadedErrorBar(DfU_step_down.x__fr_nM_, DfU_step_down.AveDfU, DfU_step_down.StdError, 'lineprops', {'-r.','linewidth',3,'markersize',18}); hold on;
shadedErrorBar(DfU_step_up.x__fr_nM_, DfU_step_up.AveDfU, DfU_step_up.StdError, 'lineprops', {'.-k','linewidth',3,'markersize',18});
ax1.Box = 'off';
xline(ax1, 20, ':', 'LineWidth', 4);
xlim(ax1, [0 20]); set(gca,'fontsize',18)
ylabel('Deviation from uniformity')
xlabel('Pheromone (nM)')

% Second subplot (ax2)
ax2 = nexttile;
shadedErrorBar(DfU_step_down.x__fr_nM_, DfU_step_down.AveDfU, DfU_step_down.StdError, 'lineprops', {'-r.','linewidth',3,'markersize',18}); hold on;
shadedErrorBar(DfU_step_up.x__fr_nM_, DfU_step_up.AveDfU, DfU_step_up.StdError, 'lineprops', {'.-k','linewidth',3,'markersize',18});
xline(ax2, 45, ':', 'LineWidth', 4);
ax2.YAxis.Visible = 'off';
ax2.Box = 'off';
xlim(ax2, [45 50]);
xticks(ax2, [45 50]);
ylim([0.3, Inf]);
ax2.Position(3) = ax2.Position(3) * 0.1;
ylim(ax2,[0.25,0.55])

% Link y-axes of the two subplots
linkaxes([ax1, ax2], 'y'); set(gca,'fontsize',18)



