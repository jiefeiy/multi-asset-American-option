clear, clc, close all
% pricing American option with d assets 
% poolobj = parpool('local',32);                    % decide the number of cpu cores

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
disp('The parameters are ');
disp(p);

%% pricing using Genz-Keister sparse grid quadrature
qlevel = 4;                                         % GK level <= 4
V0_GK = pricing_GK(p, type, ilevel, qlevel, beta, L);

%% pricing using Gauss-Hermite sparse grid quadrature
qlevel = 12;
V0_GH = pricing_GH(p, type, ilevel, qlevel, beta, L);

%% pricing using normal Leja sparse grid quadrature
qlevel = 5;
V0_leja = pricing_leja(p, type, ilevel, qlevel, beta, L);

%% pricing using QMC with scramble Sobol sequence
M = pow2(11);
V0_qmc = pricing_qmc(p, type, ilevel, M, beta, L);

%% pricing using preintegration strategy and QMC with scramble Sobol sequence
M_preint = pow2(5);
M_qmc = pow2(6);
V0_qmc_preint = pricing_qmc_preint(p, type, ilevel, M_preint, M_qmc, beta, L);

%% display pricing results
disp('=========================================================================');
disp('The results of pricing values are ')
disp(['1. The value using Genz-Keister sparse grid quadrature is ' num2str(V0_GK)]);
disp(['2. The value using Gauss-Hermite sparse grid quadrature is ' num2str(V0_GH)]);
disp(['3. The value using normal Leja sparse grid quadrature is ' num2str(V0_leja)]);
disp(['4. The value using QMC with scramble Sobol sequence is ' num2str(V0_qmc)]);
disp(['5. The value using preintegration strategy and QMC with scramble Sobol sequence is ' num2str(V0_qmc_preint)]);




