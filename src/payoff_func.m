function y = payoff_func(x, S0, strike, dim, type)
%%% define payoff function
% x = dim-by-numPts matrix
switch type
    case 'geoBaskPut'
        y = max( strike - S0*exp(sum(x)./dim), 0); 
    case 'arithmBaskPut'
        y = max( strike - sum( S0*exp(x) )/dim, 0);
    case 'maxCall'
        y = max([S0*exp(x) - strike; zeros(1,size(x,2))]);
    case 'minPut'
        y = max( strike - min(S0*exp(x)), 0); 
end


%% test example
% plot
% [Xtest, Ytest] = meshgrid(linspace(-1, 1, 100), linspace(-1, 1, 100));
% Ztest = payoff_func([Xtest(:),Ytest(:)]', p.S0, p.strike, p.dim, type);
% Ztest = reshape(Ztest, size(Xtest));
% mesh(Xtest, Ytest, Ztest);