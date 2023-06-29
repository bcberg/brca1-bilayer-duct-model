function dxdt = bilayerDuct_03(t,x,params,paramToChange)
% BILAYERDUCT_03 implements an ODE model of bilayer duct 
% growth with chemomechanical feedback, assuming stress reaches equilibrium
% faster than the growth dynamics
%   Inputs:
%       t (double): time at which to calculate dxdt
%       x (3x1 double): values for [B(t); L(t); Lpile(t)]
%       params (1x1 struct): parameter values; has fields pbar, Abar,
%       rBbar, rLbar, rLDbar, rBDbar, gammaA, gammap, gammaB, gammaL,
%       tauSA, aE, k, tcutoff, afterCutoff, pPile
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
gammap = params.gammap;
gammaB = params.gammaB;
gammaL = params.gammaL;
tauSA = params.tauSA;
aE = params.aE;
k = params.k;
pPile = params.pPile;

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
S = max([(L - B) / k, 0]);

SA = max(S - tauSA, 0);
SE = (B + L)^(1/3);
% SE = log(B + L + 1);
A = Abar / (1 + gammaA * SA);
p = pbar / (1 + gammap * A);
rL = rLbar / ( (1 + gammaL * S) * (1 + aE * SE) );
rB = rBbar / ( (1 + gammaB * A) * (1 + aE * SE) );

Lgrowth = 2*(1-p) * rB * B + (rL - rLDbar) * L;
% dLpile = pPile * max([Lgrowth,0]);
dLpile = pPile * rL * L;
dxdt = [(2*p - 1) * rB * B - rBDbar * B;
    Lgrowth - dLpile;
    dLpile];
end