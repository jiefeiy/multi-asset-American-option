clear, clc
%%% Plotting defaults %%%
[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();

MaxDepth = 4;

knots = @(n) knots_CC(n,-1,1);                                     % knots type: CC
S = mysmolyak_grid(2, MaxDepth, knots, @lev2knots_doubling);   % generate grids 
Sr = myreduce_sparse_grid(S);                                      % reduced grids
num_iknots = Sr.size; 
num_iknots_in = Sr.size_in;
zknots_in = Sr.knots_in; 
zknots = Sr.knots;

figure(1);
scatter(zknots(1,:), zknots(2,:),[], 'k',"filled");
title('CGL grid',['Points: ' num2str(num_iknots) ' Level:' num2str(MaxDepth)]);

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
ax.XTick = [-1,0,1];
ax.YTick = [-1,0,1];
exportgraphics(ax,'fig_grid_CGL.eps','Resolution',300)


%%
figure(2);
scatter(zknots_in(1,:), zknots_in(2,:),[], 'k',"filled");
title('Inner CGL grid', ['Points: ' num2str(num_iknots_in) ' Level:' num2str(MaxDepth)]);

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
ax.XTick = [-1,0,1];
ax.YTick = [-1,0,1];
exportgraphics(ax,'fig_grid_inner.eps','Resolution',300)