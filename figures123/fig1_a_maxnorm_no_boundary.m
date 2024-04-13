clear, clc, close all 
%%% requirement: Algorithm 847: Spinterp: piecewise multilinear hierarchical sparse grid interpolation in MATLAB 
% which can be found in https://dl.acm.org/doi/abs/10.1145/1114268.1114275.

[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();
MaxDepth = 4;
% sp grids
range = [-1, 1; -1, 1];
xknots = cell(1,MaxDepth+1);
numPts = cell(1,MaxDepth+1);
for l = 0:MaxDepth
    [xknots{1,l+1},numPts{1,l+1}] = mygengrid(l,2,range);  % sp grids in the range
    scatter(xknots{1,l+1}(:,1), xknots{1,l+1}(:,2), [], 'k', 'filled'); hold on;
end
numTotal = sum(cell2mat(numPts));

title('Maximum-norm-based no boundary grid',['Points: ' num2str(numTotal) ' Level:' num2str(MaxDepth)]);

box on;
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
exportgraphics(ax,'fig_grid_noboundary.eps','Resolution',300)


function [X,N] = mygengrid(k,d,range)
options = [];
options = spset(options, 'GridType', 'noboundary');
x = spgrid(k,d,options);
X = zeros(size(x));
for l = 1:d
    X(:,l) = x(:,l).*(range(l,2)-range(l,1)) + range(l,1);
end
N = size(X,1);
end