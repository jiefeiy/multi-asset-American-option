clear, clc, close all
%%% Plotting defaults %%%
[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();

ref = 3.183099555546359;
% test with different transformation parameter L
L = [1   1.5   2   2.5   3   3.5   4   4.5   5   5.5   6   6.5   7];
% ilevel = 7;
% dimension = 2;
% beta = 2;
price = [3.184247431472246   3.183962337722752   3.183324528947399   3.183519286035514...
    3.183718955693198   3.183101461474545   3.183088831206060   3.183142453428990...
    3.183081843029207   3.183196001668409   3.183512039018474   3.183655485441005...
    3.184078239557709];

err = abs(price - ref)./ref;

plot(L,err,...
    'Color',colors{1}, 'Marker' ,markers{1},'markersize',ms,'MarkerFaceColor',colors{1},'MarkerEdgeColor',colors{1},'LineWidth',lw);

xlabel('$L$','Interpreter','latex');
ylabel('relative error','Interpreter','latex');

box on;
axis tight
axis square;

ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.FontSize = fs;
ax.LineWidth = lw;
ax.FontWeight = 'bold';
ax.FontName = 'Times New Roman';
ax.XTick = [1 2 3 4 5 6 7];
exportgraphics(ax,'fig_stable_L.eps','Resolution',300)

