function dxdt = bilayerDuct_04(t,x,params,paramToChange)
% BILAYERDUCT_04 implements an ODE model of bilayer duct 
% growth with chemomechanical feedback, assuming stress reaches equilibrium
% faster than the growth dynamics
%   Inputs:
%       t (double): time at which to calculate dxdt
%       x (3x1 double): values for [B(t); L(t); Lpile(t)]
%       params (1x1 struct): parameter values; has fields pbar, Abar,
%       rBbar, rLbar, rLDbar, rBDbar, gammaA, gammap, gammaB, gammaL,
%       tauSA, aE, k, tcutoff, afterCutoff, pPile, gammaPile, gammaRank,
%       nRank, gammaAtilde, alpha, rPDbar
%       paramToChange (string): name of parameter to change to the value
%       afterCutoff at time tcutoff
%   Outputs:
%       dxdt (3x1 double): time derivative of x
pbar = params.pbar;
Abar = params.Abar;
rBbar = params.rBbar;
rLbar = params.rLbar;
rLDbar = params.rLDbar;
rBDbar = params.rBDbar;
gammaA = params.gammaA;
nA = params.nA;     % new parameter
gammap = params.gammap;
gammaB = params.gammaB;
gammaL = params.gammaL;
tauSA = params.tauSA;
aE = params.aE;
k = params.k;
pPile = params.pPile;

% qbar = params.qbar;     % new parameter
gammaPile = params.gammaPile;   % new parameter
gammaRank = params.gammaRank;   % new parameter
nRank = params.nRank;
gammaAtilde = params.gammaAtilde;   % new parameter
alpha = params.alpha;       % new parameter
% gammaq = params.gammaq;     % new parameter
% deltaq = params.deltaq;     % new parameter
rPDbar = params.rPDbar;

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

S = max([(L - B) / k, 0]);
SA = max(S - tauSA, 0);
SE = (B + L);
% SE = log(B + L + 1);
A = Abar / (1 + (gammaA * SA)^nA);
p = pbar / (1 + gammap * A);
rL = rLbar / ( (1 + gammaL * S) * (1 + aE * SE) );
rB = rBbar / ( (1 + gammaB * A) * (1 + aE * SE) );
Atilde = A * (1 + alpha * (gammaAtilde * Lpile) / (1 + gammaAtilde * Lpile));
rPile = rLbar / (1 + gammaPile * Atilde / (1 + (gammaRank * Lpile)^nRank));
% q = qbar * gammaq * Lpile / ((1 + gammaq * Lpile) * (1 + deltaq * Lpile));

dxdt = [(2*p - 1) * rB * B - rBDbar * B;
    2*(1-p) * rB * B + (rL - rLDbar) * L - pPile * L;
    pPile * L + (rPile - rPDbar) * Lpile];
end