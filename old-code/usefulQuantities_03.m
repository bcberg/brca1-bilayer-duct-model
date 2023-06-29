function output = usefulQuantities_03(sol,params)
% USEFULQUANTITIES_03 computes quantities that help describe system 
% behavior from the state history, corresponds to ODE given by
% bilayerDuct_03.m
%   Inputs:
%       sol (n x 3 double): vector of cell quantities over time for system
%       params (1x1 struct): parameter values; has fields pbar, Abar,
%       rBbar, rLbar, rLDbar, rBDbar, gammaA, gammap, gammaB, gammaL,
%       tauSA, aE, k, tCutoff, afterCutoff, pPile
%   Outputs:
%       output (1x1 struct): time courses; has fields S, SA, A, p, SE, rL,
%       rB
n = size(sol,1);
output.SI = max((sol(:,2) - sol(:,1)) / params.k, zeros(n,1));
output.SA = max(output.SI - params.tauSA, zeros(size(output.SI)));
output.A = params.Abar ./ (1 + params.gammaA * output.SA);
output.p = params.pbar ./ (1 + params.gammap * output.A);
output.SE = (sol(:,1) + sol(:,2)).^(1/3);
output.rL = params.rLbar ./ ( (1 + params.gammaL * output.SI).*(1 + ...
    params.aE * output.SE) );
output.rB = params.rBbar ./ ( (1 + params.gammaB * output.A).*(1 + ...
    params.aE * output.SE) );
end