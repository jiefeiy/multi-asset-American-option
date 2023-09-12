clear, clc, close all
addpath ../src_figures
%%% Plotting defaults %%%
[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();

h = figure(1);
%% GH
qlevel_GH = 5:16;
nb_quad_GH = [89         137         201         281         381 ...
    501         645         813        1009        1233 ...
    1489       1777];
maxerror_GH = [0.035622007802824   0.056285885739402   0.026295436347631   0.043635591760874...
    0.020839318339188   0.035605322689474   0.017258082946806   0.030055509400462...
    0.014726996079518   0.025991273542133   0.012843210167045   0.022887026130985];
hold on;
plot(nb_quad_GH, maxerror_GH,'Color',colors{1}, 'Marker' ,markers{1},...
    'markersize',ms,'MarkerFaceColor',colors{1},'MarkerEdgeColor',colors{1},'LineWidth',lw,'DisplayName','Gauss Hermite');

%% normal leja
qlevel_leja = 3:7;
nb_quad_leja = [29    65   145   321   705];
maxerr_leja = [0.074017162953470   0.025321628102900   0.028336642677640...
    0.006774116804772   0.005590932483823];
hold on;
plot(nb_quad_leja,maxerr_leja,'Color',colors{3}, 'Marker' ,markers{3},...
    'markersize',ms,'MarkerFaceColor',colors{3},'MarkerEdgeColor',colors{3},'LineWidth',lw,'DisplayName','normal Leja');

%% GK
qlevel_GK = 2:4;
nb_quad_GK = [21    65   173];
maxerror_GK = [0.013122442603773   0.023705397969289   0.005148659781106];
hold on;
plot(nb_quad_GK, maxerror_GK,'Color',colors{5}, 'Marker' ,markers{5},...
    'markersize',ms,'MarkerFaceColor',colors{5},'MarkerEdgeColor',colors{5},'LineWidth',lw,'DisplayName','Genz-Keister');

xla = xlabel('$M$');
set(xla,'Interpreter','latex');
set(xla,'FontSize',20);
yla = ylabel('$\|F_{K-1} - \hat{F}_{K-1}\|_\infty$');
set(yla,'Interpreter','latex');
set(yla,'FontSize',20);
legend('FontSize',20);
box on;
ax = gca;  ax.FontSize = 20; 

set(h, 'Position', [100 100 700 650]);
