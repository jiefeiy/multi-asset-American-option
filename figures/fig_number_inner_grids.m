clear, clc, close all
addpath ../src_figures
addpath ../src
%%% Plotting defaults %%%
[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();

MaxDepth = 5;
d = 2:10;
N_full = (pow2(MaxDepth)-1).^d;
N_CC = [ 145 441 1105 2433 4865 9017 15713 26017 41265];
N = [ 81 151 241 351 481 631 801 991 1201];

h1 = figure(1);
plot(d', N_full', 'Color',colors{1}, 'Marker' ,markers{1},'markersize',ms,'MarkerFaceColor',colors{1},'MarkerEdgeColor',colors{1},'LineWidth',lw);
hold on; plot(d', N_CC', 'Color',colors{2}, 'Marker',markers{2},'markersize',ms,'MarkerFaceColor',colors{2},'MarkerEdgeColor',colors{2},'LineWidth',lw);
hold on; plot(d', N', 'Color',colors{3}, 'Marker',markers{3},'markersize',ms,'MarkerFaceColor',colors{3},'MarkerEdgeColor',colors{3},'LineWidth',lw);

set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
xla = xlabel('$d$');
set(xla,'Interpreter','latex');
set(xla,'FontSize',20);
% ylabel('The number of points','Fontsize', 12);
box on;
ax = gca;  ax.FontSize = 20; 
leg1 = legend('$\tilde{N}_{full}$','$\tilde{N}_{CGL}$','$N$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',20);
set(leg1, 'Location', 'northwest')
set(h1, 'Position', [100 100 700 650])