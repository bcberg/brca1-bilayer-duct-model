% model01.m
% Brady Berg
clear; close all;
format long; format compact;
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')

%% System definition

pbar = 1;
Abar = 0.1;
rBbar = 1;
rLbar = 1.1;
rLDbar = 1;
rBDbar = 0;
scaleParams = [pbar; Abar; rBbar; rLbar; rLDbar; rBDbar];

gammaA = 10;
gammap = 10;
gammaB = 40;
gammaL = 0; % currently varying this parameter; this is just a placeholder
gammaLbase = 2.5;
tauSA = 0.1;
aE = 0.05;
k = 2;
feedbackParams = [gammaA; gammap; gammaB; gammaL; gammaLbase; tauSA; aE; k];
tcutoff = 30;

%% Recreating preliminary results

tspan = [0,60];
init = [1;1;0];

feedback0_25 = feedbackParams;
feedback0_25(4) = 0.25;
sys0_25 = @(t,x) prelimModel(t,x,tcutoff,init,scaleParams,feedback0_25);
feedback1 = feedbackParams;
feedback1(4) = 1;
sys1 = @(t,x) prelimModel(t,x,tcutoff,init,scaleParams,feedback1);
feedback2_5 = feedbackParams;
feedback2_5(4) = 2.5;
sys2_5 = @(t,x) prelimModel(t,x,tcutoff,init,scaleParams,feedback2_5);

[t0_25, sol0_25] = ode45(sys0_25, tspan, init);
[t1, sol1] = ode45(sys1, tspan, init);
[t2_5, sol2_5] = ode45(sys2_5, tspan, init);

% extracting basal cell self-renewal p
sA0_25 = max(sol0_25(:,3) - tauSA,zeros(size(sol0_25(:,3))));
A0_25 = Abar ./ (1 + gammaA * sA0_25);
p0_25 = pbar ./ (1 + gammap * A0_25);
sA1 = max(sol1(:,3) - tauSA,zeros(size(sol1(:,3))));
A1 = Abar ./ (1 + gammaA * sA1);
p1 = pbar ./ (1 + gammap * A1);
sA2_5 = max(sol2_5(:,3) - tauSA,zeros(size(sol2_5(:,3))));
A2_5 = Abar ./ (1 + gammaA * sA2_5);
p2_5 = pbar ./ (1 + gammap * A2_5);

fig1 = figure(1);
set(fig1,'Position',[200, 200, 900, 500])
fig1 = tiledlayout(1,2,'TileSpacing','compact');
nexttile(1)
hold on;
plot(t0_25,sol0_25(:,1),'--','Color',"#77AC30",'DisplayName','Basal0.25')
plot(t0_25,sol0_25(:,2),'-','Color',"#77AC30",'DisplayName','Lum0.25')
plot(t1,sol1(:,1),'--','Color',"#EDB120",'DisplayName','Basal1')
plot(t1,sol1(:,2),'-','Color',"#EDB120",'DisplayName','Lum1')
plot(t2_5,sol2_5(:,1),'--','Color',"#0072BD",'DisplayName','Basal2.5')
plot(t2_5,sol2_5(:,2),'-','Color',"#0072BD",'DisplayName','Lum2.5')
hold off;
lg = legend;
lg.NumColumns = 3;
lg.Layout.Tile = 'south';
title('A')
xlabel('Time')
ylabel('Cell Number')
nexttile(2)
hold on;
plot(t0_25,p0_25,'-','Color',"#77AC30")
plot(t1,p1,'-','Color',"#EDB120")
plot(t2_5,p2_5,'-','Color',"#0072BD")
hold off;
title('B')
xlabel('Time')
ylabel('Basal Cell Self-renewal')

%% Functions

function dxdt = prelimModel(t,x,tcutoff,x0,scaleParams,feedbackParams)
% PRELIMMODEL implements the preliminary ODE model of bilayer duct growth with
% chemomechanical feedback
%   Inputs:
%       t (double): time at which to calculate dxdt
%       x (3x1 double): values for [B(t); L(t); S(t)]
%       tcutoff (double): time at which to switch to gammaL = gammaLbase
%       x0 (3x1 double): values for [B(0); L(0); S(0)]
%       scaleParams (6x1 double): parameters [pbar; Abar; rBbar; rLbar;
%       rLDbar; rBDbar]
%       feedbackParams (8x1 double): parameters [gammaA; gammap; gammaB; 
%       gammaL; gammaLbase; tauSA; aE; k];
%   Outputs:
%       dxdt (3x1 double): time derivative of x
pbar = scaleParams(1);
Abar = scaleParams(2);
rBbar = scaleParams(3);
rLbar = scaleParams(4);
rLDbar = scaleParams(5);
rBDbar = scaleParams(6);
gammaA = feedbackParams(1);
gammap = feedbackParams(2);
gammaB = feedbackParams(3);
if t < tcutoff
    gammaL = feedbackParams(4);
else
    gammaL = feedbackParams(5);
end
tauSA = feedbackParams(6);
aE = feedbackParams(7);
k = feedbackParams(8);

B = x(1);
L = x(2);
S = x(3);

sA = max(S - tauSA, 0);
sE = (B + L) / (x0(1) + x0(2));
A = Abar / (1 + gammaA * sA);
p = pbar / (1 + gammap * A);
rL = rLbar / ( (1 + gammaL * S) * (1 + aE * sE) );
rB = rBbar / ( (1 + gammaB * A) * (1 + aE * sE) );
dxdt = [(2*p - 1) * rB * B - rBDbar * B;
    2*(1-p) * rB * B + (rL - rLDbar) * L;
    L - B - k*S];
end