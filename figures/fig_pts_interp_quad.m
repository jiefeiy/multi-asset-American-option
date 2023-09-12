clear, clc, close all
addpath ../src_figures
addpath ../src
%%% Plotting defaults %%%
[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();

% set parameters
p.S0 = 100; p.strike = 100; p.rate = 0.03; p.dividend = 0;
p.volatility = 0.2; p.expiration = 0.25;
p.dim = 2;                                          % asset number
p.correlation = 0.5*eye(p.dim) + 0.5*ones(p.dim);
p.numTimeStep = 50;
type = 'geoBaskPut';                                % 'geoBaskPut' or 'arithmBaskPut'
ilevel = 4;                                         % interpolation level
beta = 1;                                           % bubble parameter
L = 2;                                              % scale parameter

dt = p.expiration./p.numTimeStep;
disc = exp(-p.rate*dt);
MaxDepth = 4;                                         % interpolation level
qMaxDepth = 3;                                        % GK level <= 4
% bubble function 
bubble = @(x) prod(1-x.^2).^beta;
% compute transformation                                                           % scale parameter
[Q,Lambda] = eig(p.volatility^2*p.correlation);                    % eigen pairs of \Sigma*P*\Sigma^\top
% The order of eigenvalues do not impact the sparse grid quadrature, but will impact QMC.
qmean = Q'*(p.rate - p.dividend - p.volatility^2/2)*ones(p.dim,1)*dt;
qstd = sqrt(dt)*sqrt(Lambda);

% interpolation knots %%%
knots = @(n) knots_CC(n,-1,1);                                     % knots type: CC
S = mysmolyak_grid(p.dim, MaxDepth, knots, @lev2knots_doubling);   % generate grids 
Sr = myreduce_sparse_grid(S);                                      % reduced grids
num_iknots = Sr.size; 
num_iknots_in = Sr.size_in;
zknots_in = Sr.knots_in; 
xknots_in = atanh(zknots_in)./L;
zbubble = bubble(zknots_in);
b0 = bubble(zeros(p.dim,1));

% quadrature knots %%%
yknots = cell(1,p.dim);
for k = 1:p.dim
    yknots{k} = @(n) knots_GK(n,qmean(k),qstd(k,k));
end
yS = smolyak_grid(p.dim, qMaxDepth,yknots,@lev2knots_GK);
ySr = reduce_sparse_grid(yS); 
qknots = ySr.knots;
num_qknots = ySr.size;  
% pts next step
xqknots = zeros(p.dim, num_qknots, num_iknots_in);
for n = 1:num_iknots_in
    xqknots(:,:,n) = xknots_in(:,n) + qknots;
end
zqknots = tanh(L*xqknots);
xqknots_mat = reshape(xqknots, [p.dim, num_qknots*num_iknots_in]);
zqknots_mat = reshape(zqknots, [p.dim, num_qknots*num_iknots_in]);
zqbubble = reshape(bubble(zqknots_mat), [num_qknots, num_iknots_in]);

Sknots_in = p.S0*exp(Q*xknots_in);
Sqknots_mat = p.S0*exp(Q*xqknots_mat);

%% plot
h = figure;
subplot(1,2,1);
scatter(zqknots_mat(1,:), zqknots_mat(2,:), [],colors{6},'.','DisplayName','quadrature');
hold on;
scatter(zknots_in(1,:), zknots_in(2,:),[],'k',"filled",'DisplayName','interpolation');
legend('Fontsize', 20);
xlabel('z_1','FontSize', 20); ylabel('z_2','FontSize', 20);
box on;
ax = gca;  ax.FontSize = 20; 

subplot(1,2,2);
scatter(Sqknots_mat(1,:), Sqknots_mat(2,:), [],colors{6},'.','DisplayName','quadrature');
hold on;
scatter(Sknots_in(1,:), Sknots_in(2,:),[],'k',"filled",'DisplayName','interpolation');
legend('Fontsize', 20);
xlabel('S_1', 'FontSize', 20); ylabel('S_2', 'FontSize',20);
box on;
ax = gca;  ax.FontSize = 20; 

set(h, 'Position', [100 100 1400 650]);
