clear, clc, close all
% compare relative errors using RQMC and RQMC with preintegration
% example: dim = 5
%%% Plotting defaults %%%
[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();

ref = 2.84994;
% MaxDepth = 6 

% rmse of preintegration and plain rqmc
N1 = pow2(5);
NN = pow2(3:7);
M = N1*NN;
V_preint = [2.8504, 2.8504, 2.8503, 2.8503, 2.8503;
2.8504, 2.8504, 2.8504, 2.8504, 2.8503;
2.8506, 2.8506, 2.8503, 2.8504, 2.8503;
2.8506, 2.8505, 2.8503, 2.8504, 2.8503;
2.8509, 2.8504, 2.8504, 2.8504, 2.8503;
2.8508, 2.8505, 2.8505, 2.8503, 2.8503;
2.8503, 2.8505, 2.8504, 2.8504, 2.8503;
2.8512, 2.8504, 2.8503, 2.8503, 2.8503;
2.8506, 2.8504, 2.8504, 2.8503, 2.8503;
2.8506, 2.8508, 2.8504, 2.8504, 2.8503;
2.8505, 2.8505, 2.8504, 2.8503, 2.8503;
2.8508, 2.8507, 2.8503, 2.8504, 2.8503;
2.8504, 2.8504, 2.8503, 2.8503, 2.8503;
2.8504, 2.8504, 2.8504, 2.8504, 2.8503;
2.8506, 2.8506, 2.8503, 2.8504, 2.8503;
2.8506, 2.8505, 2.8503, 2.8504, 2.8503;
2.8509, 2.8504, 2.8504, 2.8504, 2.8503;
2.8508, 2.8505, 2.8505, 2.8503, 2.8503;
2.8503, 2.8505, 2.8504, 2.8504, 2.8503;
2.8512, 2.8504, 2.8503, 2.8503, 2.8503];

V_QMC = [2.8078, 2.8804, 2.8447, 2.8471, 2.8495;
2.8459, 2.8339, 2.8482, 2.8467, 2.8485;
2.8195, 2.8457,  2.861,  2.854, 2.8476;
2.8601, 2.8485, 2.8667, 2.8474, 2.8485;
2.7918,  2.837, 2.8516, 2.8511, 2.8476;
2.8235, 2.8317, 2.8575, 2.8573, 2.8485;
2.9297, 2.8329, 2.8553, 2.8475,  2.851;
2.8199, 2.8468, 2.8382, 2.8469, 2.8504;
2.8525, 2.8709, 2.8464, 2.8505, 2.8496;
2.8463, 2.8556, 2.8428, 2.8472, 2.8506;
2.8124,  2.854, 2.8399,  2.847, 2.8495;
2.8167,   2.84, 2.8456, 2.8441, 2.8487;
2.8078, 2.8804, 2.8447, 2.8471, 2.8495;
2.8459, 2.8339, 2.8482, 2.8467, 2.8485;
2.8195, 2.8457,  2.861,  2.854, 2.8476;
2.8601, 2.8485, 2.8667, 2.8474, 2.8485;
2.7918,  2.837, 2.8516, 2.8511, 2.8476;
2.8235, 2.8317, 2.8575, 2.8573, 2.8485;
2.9297, 2.8329, 2.8553, 2.8475,  2.851;
2.8199, 2.8468, 2.8382, 2.8469, 2.8504];

rmse_rel_QMC = sqrt(sum((V_QMC-ref).^2)./size(V_QMC,1));
rmse_rel_preint = sqrt(sum((V_preint - ref).^2)./size(V_preint,1));

h1 = figure(1);
plot( log2(M(1:end)), rmse_rel_QMC(1:end), 'DisplayName', 'RQMC',...
     'Color',colors{1}, 'Marker' ,markers{1},'markersize',ms,'MarkerFaceColor',colors{1},'MarkerEdgeColor',colors{1},'LineWidth',lw);

hold on;
plot(log2(M(1:end)), rmse_rel_preint(1:end), 'DisplayName','RQMC with preintegration',...
     'Color',colors{2}, 'Marker' ,markers{2},'markersize',ms,'MarkerFaceColor',colors{2},'MarkerEdgeColor',colors{2},'LineWidth',lw);
set(gca, 'YScale', 'log');

xlabel('$\log_2(M)$','Interpreter','latex');
ylabel('RMSE','Interpreter','latex');

legend('RQMC','RQMC with preintegration',...
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
exportgraphics(ax,'RMSE_preint.eps','Resolution',300)