clear, clc, close all
%%% Plotting defaults %%%
[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();

ref = 3.18310;
h = figure(1);
%% quadrature convergence GH
% MaxDepth = 6; 
qlevel = 5:16;
V2_Q = [3.187801958432400   3.178410136319316   3.185784940669255...
    3.180773902730623   3.183970370718688   3.181859987646454...
    3.183655111175553   3.182162281747087   3.183964394583793...
    3.181973706278835   3.184348420207433   3.181507416242187];
N_Q = [89   137   201   281   381   501   645   813   1009 ...   
    1233   1489   1777];
err2_Q = abs(V2_Q - ref(1))./ref(1);

plot(N_Q, err2_Q,'DisplayName','Gauss Hermite, L_I = 6',...
    'Color',colors{1}, 'Marker' ,markers{1},'markersize',ms,'MarkerFaceColor',colors{1},'MarkerEdgeColor',colors{1},'LineWidth',lw);

% MaxDepth = 7;
V2_Q_i7 = [3.183578058863261   3.183065670172444   3.183869276533816...
    3.182755478943469   3.183279269999189   3.182709530544524...
    3.183503435402705   3.182409614434933   3.183836006561798...
    3.181955303012386   3.184163121892517   3.181379836915442];
err2_Q_i7 = abs(V2_Q_i7 - ref(1))./ref(1);
N_Q_i7 = [89   137   201   281   381   501   645   813   1009 ...
    1233   1489   1777];
hold on; plot(N_Q_i7, err2_Q_i7,'DisplayName','Gauss Hermite, L_I = 7',...
    'Color',colors{2}, 'Marker' ,markers{2},'markersize',ms,'MarkerFaceColor',colors{2},'MarkerEdgeColor',colors{2},'LineWidth',lw);

%% quadrature convergence Leja
% MaxDepth = 6;
qlevel_leja = 3:7;
V2_leja = [3.176171640203659   3.181310682108717   3.181128592517177...
    3.182701659955648   3.182814879223781];
N_Q_leja = [29    65   145   321   705];

err2_leja = abs(V2_leja - ref(1))./ref(1);
hold on; plot(N_Q_leja, err2_leja,'DisplayName','normal Leja, L_I = 6',...
    'Color',colors{3}, 'Marker' ,markers{3},'markersize',ms,'MarkerFaceColor',colors{3},'MarkerEdgeColor',colors{3},'LineWidth',lw);

% MaxDepth = 7;
V2_leja_i7 = [3.181890150281967   3.181462592589240   3.181294312739215 ...
    3.182681054992591   3.182789344551375];
N_Q_leja_i7 = [29    65   145   321   705];

err2_leja_i7 = abs(V2_leja_i7 - ref(1))./ref(1);
hold on; plot(N_Q_leja_i7, err2_leja_i7,'DisplayName','normal Leja, L_I = 7',...
    'Color',colors{4}, 'Marker' ,markers{4},'markersize',ms,'MarkerFaceColor',colors{4},'MarkerEdgeColor',colors{4},'LineWidth',lw);

%% quadrature convergence GK
% qlevel_GK = 2:4;
V2_GK = [3.183642009824087   3.182175081482277   3.183899095187631];
N_Q_GK = [21    65   173];
err2_GK = abs(V2_GK - ref(1))./ref(1);
hold on; plot(N_Q_GK, err2_GK,'DisplayName','Genz-Keister, L_I = 6',...
    'Color',colors{5}, 'Marker' ,markers{5},'markersize',ms,'MarkerFaceColor',colors{5},'MarkerEdgeColor',colors{5},'LineWidth',lw);

% MaxDepth = 7;
V2_GK_i7 = [3.181722588270920   3.181729395153899   3.183145582488626];
N_Q_GK = [21    65   173];
err2_GK_i7 = abs(V2_GK_i7 - ref(1))./ref(1);
hold on; plot(N_Q_GK, err2_GK_i7, 'DisplayName','Genz-Keister, L_I = 7',...
    'Color',colors{6}, 'Marker' ,markers{6},'markersize',ms,'MarkerFaceColor',colors{6},'MarkerEdgeColor',colors{6},'LineWidth',lw);

xlabel('$M$','Interpreter','latex');
ylabel('relative error','Interpreter','latex');

legend('Gauss Hermite, $L_I = 6$','Gauss Hermite, $L_I = 7$',...
    'normal Leja, $L_I = 6$', 'normal Leja, $L_I = 7$', ...
    'Genz-Keister, $L_I = 6$', 'Genz-Keister, $L_I = 7$',...
    'Interpreter','latex','FontSize',fs);

box on;
axis square;

ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.FontSize = fs;
ax.LineWidth = lw;
ax.FontWeight = 'bold';
ax.FontName = 'Times New Roman';
exportgraphics(ax,'fig_err_quad_b.eps','Resolution',300)