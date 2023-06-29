function output = usefulQuantities_02(sol,params)
% USEFULQUANTITIES_02 computes quantities that help describe system 
% behavior from the state history, assuming stress equilibrates faster
% than the growth dynamics
%   Inputs:
%       sol (n x 2 double): vector of cell quantities over time for system
%       params (1x1 struct): parameter values; has fields pbar, Abar,
%       rBbar, rLbar, rLDbar, rBDbar, gammaA, gammap, gammaB, gammaL,
%       tauSA, aE, k, tCutoff, afterCutoff
%   Outputs:
%       output (1x1 struct): time courses; has fields S, SA, A, p, SE, rL,
%       rB
output.S = (sol(:,2) - sol(:,1)) / params.k;
output.SA = max(output.S - params.tauSA, zeros(size(output.S)));
output.A = params.Abar ./ (1 + params.gammaA * output.SA);
output.p = params.pbar ./ (1 + params.gammap * output.A);
output.SE = sum(sol,2) ./ (sol(1,1) + sol(1,2));
output.rL = params.rLbar ./ ( (1 + params.gammaL * output.S).*(1 + ...
    params.aE * output.SE) );
output.rB = params.rBbar ./ ( (1 + params.gammaB * output.A).*(1 + ...
    params.aE * output.SE) );
end