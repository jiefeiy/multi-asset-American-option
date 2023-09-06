function V0 = pricing_leja(p, type, ilevel, qlevel, beta, L)
% pricing using normal leja sparse grid quadrature
disp('--------------------------------------------------')
disp('pricing using normal leja sparse grid quadrature');
dt = p.expiration./p.numTimeStep;
disc = exp(-p.rate*dt);
MaxDepth = ilevel;                                         % interpolation level
qMaxDepth = qlevel;                                        % GK level <= 4
% bubble function 
bubble = @(x) prod(1-x.^2).^beta;
% compute transformation                                                           % scale parameter
[Q,Lambda] = eig(p.volatility^2*p.correlation);                    % eigen pairs of \Sigma*P*\Sigma^\top
% The order of eigenvalues do not impact the sparse grid quadrature, but will impact QMC.
qmean = Q'*(p.rate - p.dividend - p.volatility^2/2)*ones(p.dim,1)*dt;
qstd = sqrt(dt)*sqrt(Lambda);

%% interpolation knots %%%
knots = @(n) knots_CC(n,-1,1);                                     % knots type: CC
S = mysmolyak_grid(p.dim, MaxDepth, knots, @lev2knots_doubling);   % generate grids 
Sr = myreduce_sparse_grid(S);                                      % reduced grids
num_iknots = Sr.size; 
num_iknots_in = Sr.size_in;
zknots_in = Sr.knots_in; 
xknots_in = atanh(zknots_in)./L;
zbubble = bubble(zknots_in);
b0 = bubble(zeros(p.dim,1));

%% quadrature knots %%%
yknots = cell(1,p.dim);
for k = 1:p.dim
    yknots{k} = @(n) knots_normal_leja(n,qmean(k),qstd(k,k),'line');
end
yS = smolyak_grid(p.dim, qMaxDepth,yknots,@lev2knots_doubling);
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

%% dynamic programming %%%
payoff = payoff_func(Q*xqknots_mat, p.S0, p.strike, p.dim, type);
payoff = reshape(payoff, [num_qknots, num_iknots_in]);

EV = payoff; 
tic
for k = p.numTimeStep-1:-1:0 
    udata_in = disc*(ySr.weights*payoff).*zbubble; 
    
    if k == 0
        u0 = myinterpolate_on_sparse_grid( S,Sr, udata_in, zeros(p.dim, 1) );
        break;
    end

    u_temp = zeros(num_qknots, num_iknots_in);
    parfor n = 1:num_iknots_in
        u_temp(:,n) = myinterpolate_on_sparse_grid( S,Sr, udata_in, zqknots(:,:,n) )';
    end
    payoff = max( EV, u_temp./zqbubble );

    disp(['Computing k = ' num2str(k) '......']);
end
V0 = u0./b0
toc

disp(['The number of inner sparse grids N = ' num2str(num_iknots_in)]);
disp(['The number of quadrature points M = ' num2str(num_qknots)]);
disp(['The number of sparse grids (include boundary points) is ' num2str(num_iknots)]);
