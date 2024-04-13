clear, close all
[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();

%% GH test
qlevel = 5:16;
nb_qlevel = size(qlevel,2);
reerr = zeros(1,nb_qlevel);
nb_quad = zeros(1,nb_qlevel);

for i = 1:nb_qlevel
    qMaxDepth = qlevel(i);
    [reerr(i), nb_quad(i)] = GH(qMaxDepth);
end

plot(nb_quad, reerr, 'Color',colors{1}, 'Marker' ,markers{1},...
    'markersize',ms,'MarkerFaceColor',colors{1},'MarkerEdgeColor',colors{1},'LineWidth',lw,'DisplayName','Gauss Hermite');
hold on;

%% normal leja test
qlevel = 3:7;
nb_qlevel = size(qlevel,2);
reerr = zeros(1,nb_qlevel);
nb_quad = zeros(1,nb_qlevel);

for i = 1:nb_qlevel
    qMaxDepth = qlevel(i);
    [reerr(i), nb_quad(i)] = normal_leja(qMaxDepth);
end

plot(nb_quad, reerr, 'Color',colors{3}, 'Marker' ,markers{3},...
    'markersize',ms,'MarkerFaceColor',colors{3},'MarkerEdgeColor',colors{3},'LineWidth',lw,'DisplayName','normal Leja');
hold on;

%% GK test
qlevel = 2:4;
nb_qlevel = size(qlevel,2);
reerr = zeros(1,nb_qlevel);
nb_quad = zeros(1,nb_qlevel);

for i = 1:nb_qlevel
    qMaxDepth = qlevel(i);
    [reerr(i), nb_quad(i)] = GK(qMaxDepth);
end

plot(nb_quad, reerr, 'Color',colors{5}, 'Marker' ,markers{5},...
    'markersize',ms,'MarkerFaceColor',colors{5},'MarkerEdgeColor',colors{5},'LineWidth',lw,'DisplayName','Genz-Keister');
hold on;

%% figure setting
xlabel('$M$','Interpreter','latex');
ylabel('relative $L^\infty$ error', 'Interpreter','latex');

legend('Gauss Hermite','normal Leja','Genz-Keister','Interpreter','latex','Location','northeast','FontSize',fs);

box on;
axis square;

ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.FontSize = fs;
ax.LineWidth = lw;
ax.FontWeight = 'bold';
ax.FontName = 'Times New Roman';
exportgraphics(ax,'fig_err_quad_1step.eps','Resolution',300)


%%
function [err, num_qknots] = GH(qMaxDepth)
    % geometric basket put %%%
    p.S0 = 100; p.strike = 100; p.rate = 0.03; p.dividend = 0;
    p.volatility = 0.2; p.expiration = 0.25;
    p.dim = 2;
    p.correlation = 0.5*eye(p.dim) + 0.5*ones(p.dim);
    p.numTimeStep = 50;
    dt = p.expiration./p.numTimeStep;
    disc = exp(-p.rate*dt);
    type = 'geoBaskPut';
    % compute transformation
    [Q,Lambda] = eig(p.volatility^2*p.correlation);                  % eigen pairs of \Sigma*P*\Sigma^\top
    % The order of eigenvalues do not impact the sparse grid quadrature, but
    % will impact QMC.
    qmean = Q'*(p.rate - p.dividend - p.volatility^2/2)*ones(p.dim,1)*dt;
    qstd = sqrt(dt)*sqrt(Lambda);
    
    %%% generate some random knots at t_{K-1} %%%
    t = (p.numTimeStep-1)*dt;
    num_pts = 1000;
    xx = Q'*(p.rate - p.dividend - p.volatility^2/2)*ones(p.dim,1)*t + sqrt(t)*sqrt(Lambda)*randn(p.dim,num_pts);
    
    %%% quadrature knots %%%
    yknots = cell(1,p.dim);
    for k = 1:p.dim
        yknots{k} = @(n) knots_normal(n,qmean(k),qstd(k,k));
    end
    yS = smolyak_grid(p.dim, qMaxDepth,yknots,@lev2knots_lin);
    ySr = reduce_sparse_grid(yS); 
    qknots = ySr.knots;
    num_qknots = ySr.size;  
    % pts next step
    xqknots = zeros(p.dim, num_qknots, num_pts);
    for n = 1:num_pts
        xqknots(:,:,n) = xx(:,n) + qknots;
    end
    xqknots_mat = reshape(xqknots, [p.dim, num_qknots*num_pts]);
    
    %%% 1d
    xknots_1d = mean(Q*xx);
    sig = 0.2./p.dim*sqrt(0.5*p.dim^2 + 0.5*p.dim); 
    DivYield = 0.2^2/2 - sig^2/2;
    [~,C_bs] = blsprice(p.S0*exp(xknots_1d),p.strike,p.rate,dt,sig, DivYield);
    
    %%% first step (equivalent to European option)
    payoff = payoff_func(Q*xqknots_mat, p.S0, p.strike, p.dim, type);
    payoff = reshape(payoff, [num_qknots, num_pts]);
    C_ap = disc*(ySr.weights*payoff);
    
    err = max( abs(C_ap - C_bs) )/p.strike;
end

function [err, num_qknots] = normal_leja(qMaxDepth)
    % geometric basket put %%%
    p.S0 = 100; p.strike = 100; p.rate = 0.03; p.dividend = 0;
    p.volatility = 0.2; p.expiration = 0.25;
    p.dim = 2;
    p.correlation = 0.5*eye(p.dim) + 0.5*ones(p.dim);
    p.numTimeStep = 50;
    dt = p.expiration./p.numTimeStep;
    disc = exp(-p.rate*dt);
    type = 'geoBaskPut';
    % compute transformation
    [Q,Lambda] = eig(p.volatility^2*p.correlation);                  % eigen pairs of \Sigma*P*\Sigma^\top
    % The order of eigenvalues do not impact the sparse grid quadrature, but
    % will impact QMC.
    qmean = Q'*(p.rate - p.dividend - p.volatility^2/2)*ones(p.dim,1)*dt;
    qstd = sqrt(dt)*sqrt(Lambda);
    
    %%% generate some random knots at t_{K-1} %%%
    t = (p.numTimeStep-1)*dt;
    num_pts = 1000;
    xx = Q'*(p.rate - p.dividend - p.volatility^2/2)*ones(p.dim,1)*t + sqrt(t)*sqrt(Lambda)*randn(p.dim,num_pts);
    
    %%% quadrature knots %%%
    yknots = cell(1,p.dim);
    for k = 1:p.dim
        yknots{k} = @(n) knots_normal_leja(n,qmean(k),qstd(k,k),'line');
    end
    yS = smolyak_grid(p.dim, qMaxDepth,yknots,@lev2knots_doubling);
    ySr = reduce_sparse_grid(yS); 
    qknots = ySr.knots;
    num_qknots = ySr.size;  
    % pts next step
    xqknots = zeros(p.dim, num_qknots, num_pts);
    for n = 1:num_pts
        xqknots(:,:,n) = xx(:,n) + qknots;
    end
    xqknots_mat = reshape(xqknots, [p.dim, num_qknots*num_pts]);
    
    %%% 1d
    xknots_1d = mean(Q*xx);
    sig = 0.2./p.dim*sqrt(0.5*p.dim^2 + 0.5*p.dim); 
    DivYield = 0.2^2/2 - sig^2/2;
    [~,C_bs] = blsprice(p.S0*exp(xknots_1d),p.strike,p.rate,dt,sig, DivYield);
    
    %%% first step (equivalent to European option)
    payoff = payoff_func(Q*xqknots_mat, p.S0, p.strike, p.dim, type);
    payoff = reshape(payoff, [num_qknots, num_pts]);
    C_ap = disc*(ySr.weights*payoff);
    
    err = max( abs(C_ap - C_bs) )/p.strike;
end

function [err, num_qknots] = GK(qMaxDepth)
    % geometric basket put %%%
    p.S0 = 100; p.strike = 100; p.rate = 0.03; p.dividend = 0;
    p.volatility = 0.2; p.expiration = 0.25;
    p.dim = 2;
    p.correlation = 0.5*eye(p.dim) + 0.5*ones(p.dim);
    p.numTimeStep = 50;
    dt = p.expiration./p.numTimeStep;
    disc = exp(-p.rate*dt);
    type = 'geoBaskPut';
    % compute transformation
    [Q,Lambda] = eig(p.volatility^2*p.correlation);                  % eigen pairs of \Sigma*P*\Sigma^\top
    % The order of eigenvalues do not impact the sparse grid quadrature, but
    % will impact QMC.
    qmean = Q'*(p.rate - p.dividend - p.volatility^2/2)*ones(p.dim,1)*dt;
    qstd = sqrt(dt)*sqrt(Lambda);
    
    %%% generate some random knots at t_{K-1} %%%
    t = (p.numTimeStep-1)*dt;
    num_pts = 1000;
    xx = Q'*(p.rate - p.dividend - p.volatility^2/2)*ones(p.dim,1)*t + sqrt(t)*sqrt(Lambda)*randn(p.dim,num_pts);
    
    %%% quadrature knots %%%
    yknots = cell(1,p.dim);
    for k = 1:p.dim
        yknots{k} = @(n) knots_GK(n,qmean(k),qstd(k,k));
    end
    yS = smolyak_grid(p.dim, qMaxDepth,yknots,@lev2knots_GK);
    ySr = reduce_sparse_grid(yS); 
    qknots = ySr.knots;
    num_qknots = ySr.size;  
    % pts next step
    xqknots = zeros(p.dim, num_qknots, num_pts);
    for n = 1:num_pts
        xqknots(:,:,n) = xx(:,n) + qknots;
    end
    xqknots_mat = reshape(xqknots, [p.dim, num_qknots*num_pts]);
    
    %%% 1d
    xknots_1d = mean(Q*xx);
    sig = 0.2./p.dim*sqrt(0.5*p.dim^2 + 0.5*p.dim); 
    DivYield = 0.2^2/2 - sig^2/2;
    [~,C_bs] = blsprice(p.S0*exp(xknots_1d),p.strike,p.rate,dt,sig, DivYield);
    
    %%% first step (equivalent to European option)
    payoff = payoff_func(Q*xqknots_mat, p.S0, p.strike, p.dim, type);
    payoff = reshape(payoff, [num_qknots, num_pts]);
    C_ap = disc*(ySr.weights*payoff);
    
    err = max( abs(C_ap - C_bs) )/p.strike;
end