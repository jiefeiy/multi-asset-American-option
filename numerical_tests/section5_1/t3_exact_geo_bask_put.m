clear, clc, close all
% => For multi-dimensional American geometric baskey put benchmark. 
% discounted to present first
% CC quadrature + barycentric formula at Chebyshev-Lobatto pts

d_vals = 2:16;
V0_vals = zeros(size(d_vals));
i = 1;
for d = d_vals
    V0_vals(i) = run_1d_geo_bask_put(d);
    fprintf("The reference Bermudan price of d = %d is %1.5f \n", d, V0_vals(i));
    fprintf("--------------------------------------------------\n");
    i = i+1;
end


function V0 = run_1d_geo_bask_put(d)
S0 = 100; r = 0.03; sig = 0.2./d*sqrt(0.5*d^2 + 0.5*d); 
DivYield = 0.2^2/2 - sig^2/2;
T = 0.25; kappa = 100; K = 50;

M = 400;                                 % number of interpolation pts
n = 400;                                 % number of quadrature pts
dt = T/K;
NumInterp = 5;                           % NumInterp standard error of interpolation interval
NumSig = 6;                              % NumSig standard error of quadrature interval
lb = log(S0) + (r-DivYield-sig*sig/2)*T - NumInterp*sig*sqrt(T);
ub = log(S0) + (r-DivYield-sig*sig/2)*T + NumInterp*sig*sqrt(T);

x = (ub-lb)/2*cos( pi*(0:M)/M ) + (ub+lb)/2; % interpolation pts x
C = zeros(size(x)); s = zeros(n+1, M+1);

% compute quadrature pts s(:)
for m = 1:M+1
    li = x(m) + (r-DivYield-sig*sig/2)*dt - NumSig*sig*sqrt(dt);
    ui = x(m) + (r-DivYield-sig*sig/2)*dt + NumSig*sig*sqrt(dt);
    t = pi*(0:n)'/n;
    s(:,m) = (ui-li)/2*cos(t) + (ui+li)/2;
end
V = max(kappa - exp(s(:)), 0);
V = exp(-r*T)*reshape(V, n+1, M+1);

% Backward induction
for k = K-1:-1:1
    for m = 1:M+1
        li = x(m) + (r-DivYield-sig*sig/2)*dt - NumSig*sig*sqrt(dt);
        ui = x(m) + (r-DivYield-sig*sig/2)*dt + NumSig*sig*sqrt(dt);
        t = pi*(0:n)'/n;
        ftemp = 1/sig/sqrt(2*pi*dt)*exp( -(s(:,m) - x(m) - (r-DivYield-sig*sig/2)*dt).^2/(2*dt*sig*sig) );
        C(m) = clenshaw_curtis(s(:,m), V(:,m).*ftemp, li, ui);
    end
    f = barycentric_1(x', C', lb, ub, s(:));
    g = exp(-r*k*dt)*(kappa - exp(s(:)));
    V = max([g, f], [], 2);
    V = reshape(V, n+1, M+1);
end

% compute the initial value
for m = 1:M+1
    li = x(m) + (r-DivYield-sig*sig/2)*dt - NumSig*sig*sqrt(dt);
    ui = x(m) + (r-DivYield-sig*sig/2)*dt + NumSig*sig*sqrt(dt);
    t = pi*(0:n)'/n;
    ftemp = 1/sig/sqrt(2*pi*dt)*exp( -(s(:,m) - x(m) - (r-DivYield-sig*sig/2)*dt).^2/(2*dt*sig*sig) );
    C(m) = clenshaw_curtis(s(:,m), V(:,m).*ftemp, li, ui);
end
V0 = barycentric_1(x', C', lb, ub, log(S0));
end

