clear; clc;
% Load figure data
inp_file = '../FigureData/Figure5B.mat';
load(inp_file);

% Use the sum of active Cdc42 to draw regime boundary.
% Threshold for unpolarized (<=100) vs bistable and polarized (>100)
th1 = 100;
% Threshold for bistable (<=300) vs polarized (>300)
th2 = 300;

% compare total number of active Cdc42 and the threshold
im1 = zeros(size(data_ampl_sum));
im1(data_ampl_sum > th1 ) = 1;
im1(6:11, 16) = 0; im1(5:7, 17) = 0; im1(5:6, 18) = 0; im1(5, 19:20) = 0; im1(4, 21:24) = 0;
im2 = zeros(size(data_ampl_sum));
im2(data_ampl_sum > th2) = 1;
im2 = im2 + cell_outline(im2);
im = zeros(size(im1));
im(im1 == 1) = 3; % im1 = 1, Bistable regime & Polarized regime
im(im2 == 1) = 1.3; % im2 = 1, Polarized regime

SF = 10^(-21)*6.02*10^(23);
s = size(data_ampl_sum);

fig = figure('Position', [50 50 8*s(2) 8*s(1)]);
hold on;
axis xy;
colormap hot;
imagesc(im);

xlabel('$k_1 \left( \frac{\mu m^2}{s \cdot mol^2} \right)$', 'Interpreter','latex');
x_ticks = 1:20:length(k3_vals);
x_ticklabels = string(round(k3_vals(x_ticks)/SF/SF,7));
xticks(x_ticks);
xticklabels(x_ticklabels);
xlim([0, length(k3_vals)]);

ylabel('$k_{2} \left( \frac{1}{s} \right)$', 'Interpreter','latex');
y_ticks = 1:20:length(k2b_vals);
y_ticklabels = string(round(k2b_vals(x_ticks),1));
yticks(y_ticks);
yticklabels(y_ticklabels);
ylim([0, length(k2b_vals)]);

set(gca, 'FontSize', 33);

function outline = cell_outline(Im)
%CELL_OUTLINE - returns outline of the sell
%   Im - cell matrix (1s and 0s)

outline = zeros(size(Im)+2);
Im = [zeros(1,size(Im,1)+2);zeros(size(Im,1),1), Im, zeros(size(Im,1),1);zeros(1,size(Im,1)+2)];
outline(2:end-1,2:end-1) = Im(1:end-2,2:end-1) + Im(3:end,2:end-1) + ...
    Im(2:end-1,1:end-2) + Im(2:end-1,3:end);
outline = outline - 4*Im;
o = outline>0;
o = o(2:end-1, 2:end-1);
outline = o;

end
