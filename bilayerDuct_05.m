function dxdt = bilayerDuct_05(t,x,params,paramToChange)
% BILAYERDUCT_05 implements an ODE model of bilayer duct 
% growth with chemomechanical feedback, assuming stress reaches equilibrium
% faster than the growth dynamics
%   Inputs:
%       t (double): time at which to calculate dxdt
%       x (3x1 double): values for [B(t); L(t); Lpile(t)]
%       params (1x1 struct): parameter values; has fields pbar, Abar,
%       rBbar, rLbar, rBDbar, rLDbar, gammap, gammaA, gammaB, gammaL,
%       pPile, gammaPile, gammaAP, gammaFST, nA, nFST, tauSA, aE, k, 
%       tcutoff, afterCutoff
%       paramToChange (string): name of parameter to change to the value
%       afterCutoff at time tcutoff
%   Outputs:
%       dxdt (3x1 double): time derivative of x
pbar = params.pbar;
Abar = params.Abar;
rBbar = params.rBbar;
rLbar = params.rLbar;
rBDbar = params.rBDbar;
rLDbar = params.rLDbar;
rPDbar = params.rPDbar; % r_IDbar
gammap = params.gammap;
gammaA = params.gammaA;
gammaB = params.gammaB;
gammaL = params.gammaL;
pPile = params.pPile;   % p_I
gammaPile = params.gammaPile; % gamma_I
gammaAP = params.gammaAP; % gamma_AI
gammaFST = params.gammaFST;
nA = params.nA;
nFST = params.nFST;
tauSA = params.tauSA;
aE = params.aE;
k = params.k;

if t > params.tCutoff
    switch paramToChange
        case 'gammaL'
            gammaL = params.afterCutoff;
        case 'rLbar'
            rLbar = params.afterCutoff;
        case 'aE'
            aE = params.afterCutoff;
        case 'k'
            k = params.afterCutoff;
        case 'nochange'
            % do nothing
    end
end

B = x(1);
L = x(2);
Lpile = x(3);

SI = max((L - B) / k, 0);
% SI = abs(L - B) / k;
SA = max(SI - tauSA, 0);
% SE = B + L;
% A = (Abar + gammaAP*Lpile/(1+Lpile)) / ...
A = (Abar + gammaAP*Lpile) / ...
    ((1 + (gammaA*SA)^nA)*(1 + (gammaFST*Lpile)^nFST));
p = pbar / (1 + gammap * A);
% rB = rBbar / ( (1 + gammaB * A) * (1 + aE * SE) );
% rL = rLbar / ( (1 + gammaL * A) * (1 + aE * SE) );
rB = rBbar / ( (1 + gammaB * A) * (1 + aE * B) );
rL = rLbar / ( (1 + gammaL * A) * (1 + aE * L) );
rPile = rLbar / (1 + gammaPile * A);

% dxdt = [(2*p - 1) * rB * B - rBDbar * B;
%     2*(1-p) * rB * B + (rL - rLDbar) * L - pPile * L;
%     pPile * L + (rPile - rPDbar) * Lpile];

dxdt = [(2*p - 1) * rB * B - rBDbar * B - pPile * B;
    2*(1-p) * rB * B + (rL - rLDbar) * L;
    pPile * B + (rPile - rPDbar) * Lpile];

end