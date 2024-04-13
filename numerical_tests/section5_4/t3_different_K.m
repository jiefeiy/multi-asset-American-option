clear, clc, close all
%%% Plotting defaults %%%
[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();

% convergence with respect to the number of time steps K
ref = 3.18469;

% L=2
price_i8 = [3.176067726777214   3.180118249497212   3.180998436220298   3.181382339175402...
    3.182113441358622   3.182456838147611   3.182757135342523   3.182988400955672...
    3.183109660188392   3.183205188353993   3.183380533973295   3.183471009399982...
    3.183588968742941   3.183721070600160   3.183875453296925   3.184001729476704...
    3.184091104241385   3.184161227647431   3.184206057902266   3.184227459130841];

err_i8 = abs(price_i8 - ref)./ref;
h3 = figure(3);
plot(10:5:105, err_i8,'DisplayName','L=2',...
     'Color',colors{1}, 'Marker' ,markers{1},'markersize',ms,'MarkerFaceColor',colors{1},'MarkerEdgeColor',colors{1},'LineWidth',lw);

xref = 10.^(1:.1:2);
yref = xref.^(-1)*exp(-4);
hold on; plot(xref, yref, '--k', 'DisplayName','slope = -1','LineWidth',lw);

set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');

xlabel('$K$','Interpreter','latex');
ylabel('relative error','Interpreter','latex');

legend('$L_I = 8, \beta = 1, L=2$','slope = -1',...
    'Interpreter','latex','FontSize',fs);

box on;
% axis tight
axis square;

ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.FontSize = fs;
ax.LineWidth = lw;
ax.FontWeight = 'bold';
ax.FontName = 'Times New Roman';
exportgraphics(ax,'fig_stable_K.eps','Resolution',300)