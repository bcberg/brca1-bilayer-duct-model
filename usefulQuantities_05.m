function output = usefulQuantities_05(sol,params)
% USEFULQUANTITIES_05 computes quantities that help describe system 
% behavior from the state history, corresponds to ODE given by
% bilayerDuct_05.m
%   Inputs:
%       sol (n x 3 double): vector of cell quantities over time for system
%       params (1x1 struct): parameter values; has fields pbar, Abar,
%       rBbar, rLbar, rBDbar, rLDbar, gammap, gammaA, gammaB, gammaL,
%       pPile, gammaPile, gammaAP, gammaFST, nA, nFST, tauSA, aE, k, 
%       tcutoff, afterCutoff
%   Outputs:
%       output (1x1 struct): time courses; has fields SI, SA, SE, A, p,
%       rB, rL, rPile
n = size(sol,1);
B = sol(:,1);
L = sol(:,2);
Lpile = sol(:,3);
output.SI = max((L - B) / params.k, zeros(n,1));
% output.SI = abs(L - B) / params.k;
output.SA = max(output.SI - params.tauSA, zeros(n,1));
% SE is no longer used to feed back on rB or rL
output.SE = B + L;
% output.rPD = params.rPDbar * (1 + params.alpha * params.gammaPD * sol(:,3)./ ...
%     (1 + params.gammaPD * sol(:,3)));
output.A = (params.Abar + params.gammaAP * Lpile) ./ ...
    ((1 + (params.gammaA * output.SA).^params.nA) .* ...
    (1 + (params.gammaFST * Lpile).^params.nFST));
output.p = params.pbar ./ (1 + params.gammap * output.A);
output.rB = params.rBbar ./ ( (1 + params.gammaB * output.A).*(1 + ...
    params.aE * B) );
output.rL = params.rLbar ./ ( (1 + params.gammaL * output.A).*(1 + ...
    params.aE * L) );
output.rPile = params.rLbar ./ (1 + params.gammaPile * output.A);
end