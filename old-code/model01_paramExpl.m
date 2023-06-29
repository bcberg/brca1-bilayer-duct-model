% model01_paramExpl.m
% Brady Berg
clear; close all;
format long; format compact;
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')

%% Base parameters

pbar = 1;
Abar = 0.1;
rBbar = 1;
rLbar = 1.1;
rLDbar = 1;
rBDbar = 0;
gammaA = 10;
gammap = 10;
gammaB = 40;
gammaL = 1;
tauSA = 0.1;
aE = 0.05;
k = 2;

afterCutoff = [];
params = struct('pbar',pbar,'Abar',Abar,'rBbar',rBbar,'rLbar',rLbar,...
    'rLDbar',rLDbar,'rBDbar',rBDbar,'gammaA',gammaA,'gammap',gammap,...
    'gammaB',gammaB,'gammaL',gammaL,'tauSA',tauSA,'aE',aE,'k',k,...
    'afterCutoff',afterCutoff);
tcutoff = 30;
tspan = [0,60];
init = [1;1;0];

%% gammaL

paramToChange = 'gammaL';
params.afterCutoff = 2.5;

params0_25 = params;
params0_25.gammaL = 0.25;
sys0_25 = @(t,x) prelimModel_param(t,x,tcutoff,init,params0_25,paramToChange);

params1 = params;
params1.gammaL = 1;
sys1 = @(t,x) prelimModel_param(t,x,tcutoff,init,params1,paramToChange);

params2_5 = params;
params2_5.gammaL = 2.5;
sys2_5 = @(t,x) prelimModel_param(t,x,tcutoff,init,params2_5,paramToChange);

[t0_25, sol0_25] = ode45(sys0_25, tspan, init);
p0_25 = basalSelfRenew(sol0_25(:,3),params0_25);
[t1, sol1] = ode45(sys1, tspan, init);
p1 = basalSelfRenew(sol1(:,3),params1);
[t2_5, sol2_5] = ode45(sys2_5, tspan, init);
p2_5 = basalSelfRenew(sol2_5(:,3),params2_5);

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
title(fig1,'Changing $\gamma_L$','Interpreter','latex')
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
clear params0_25 params1 params2_5 sys0_25 sys1 sys2_5 t0_25 t1 t2_5 ...
    sol0_25 sol1 sol2_5 p0_25 p1 p2_5

%% rLbar

paramToChange = 'rLbar';
params.afterCutoff = 0.25;

paramsA = params;
paramsA.rLbar = 0.25;
sysA = @(t,x) prelimModel_param(t,x,tcutoff,init,paramsA,paramToChange);

paramsB = params;
paramsB.rLbar = 0.75;
sysB = @(t,x) prelimModel_param(t,x,tcutoff,init,paramsB,paramToChange);

paramsC = params;
paramsC.rLbar = 1.1;
sysC = @(t,x) prelimModel_param(t,x,tcutoff,init,paramsC,paramToChange);

paramsD = params;
paramsD.rLbar = 1.5;
sysD = @(t,x) prelimModel_param(t,x,tcutoff,init,paramsD,paramToChange);

[tA, solA] = ode45(sysA, tspan, init);
pA = basalSelfRenew(solA(:,3),paramsA);
[tB, solB] = ode45(sysB, tspan, init);
pB = basalSelfRenew(solB(:,3),paramsB);
[tC, solC] = ode45(sysC, tspan, init);
pC = basalSelfRenew(solC(:,3),paramsC);
[tD, solD] = ode45(sysD, tspan, init);
pD = basalSelfRenew(solD(:,3),paramsD);

fig2 = figure(2);
set(fig2,'Position',[210, 210, 900, 500])
fig2 = tiledlayout(1,2,'TileSpacing','compact');
nexttile(1)
hold on;
plot(tA,solA(:,1),'--','Color',"#77AC30",'DisplayName', ...
    ['Basal',num2str(paramsA.rLbar)])
plot(tA,solA(:,2),'-','Color',"#77AC30",'DisplayName', ...
    ['Lum',num2str(paramsA.rLbar)])
plot(tB,solB(:,1),'--','Color',"#EDB120",'DisplayName', ...
    ['Basal',num2str(paramsB.rLbar)])
plot(tB,solB(:,2),'-','Color',"#EDB120",'DisplayName', ...
    ['Lum',num2str(paramsB.rLbar)])
plot(tC,solC(:,1),'--','Color',"#0072BD",'DisplayName', ...
    ['Basal',num2str(paramsC.rLbar)])
plot(tC,solC(:,2),'-','Color',"#0072BD",'DisplayName', ...
    ['Lum',num2str(paramsC.rLbar)])
plot(tD,solD(:,1),'--','Color',"#7E2F8E",'DisplayName', ...
    ['Basal',num2str(paramsD.rLbar)])
plot(tD,solD(:,2),'-','Color',"#7E2F8E",'DisplayName', ...
    ['Lum',num2str(paramsD.rLbar)])
hold off;
title(fig2,'Changing $\bar{r}_L$','Interpreter','latex')
lg = legend;
lg.NumColumns = 4;
lg.Layout.Tile = 'south';
title('A')
xlabel('Time')
ylabel('Cell Number')
nexttile(2)
hold on;
plot(tA,pA,'-','Color',"#77AC30")
plot(tB,pB,'-','Color',"#EDB120")
plot(tC,pC,'-','Color',"#0072BD")
plot(tD,pD,'-','Color',"#7E2F8E")
hold off;
title('B')
xlabel('Time')
ylabel('Basal Cell Self-renewal')
clear paramsA paramsB paramsC paramsD sysA sysB sysC sysD tA tB tC tD ...
    solA solB solC solD pA pB pC pD

%% aE

paramToChange = 'aE';
params.afterCutoff = 0.05;

paramsA = params;
paramsA.aE = 1e-3;
sysA = @(t,x) prelimModel_param(t,x,tcutoff,init,paramsA,paramToChange);

paramsB = params;
paramsB.aE = 5e-3;
sysB = @(t,x) prelimModel_param(t,x,tcutoff,init,paramsB,paramToChange);

paramsC = params;
paramsC.aE = 2.5e-2;
sysC = @(t,x) prelimModel_param(t,x,tcutoff,init,paramsC,paramToChange);

paramsD = params;
paramsD.aE = 0.05;
sysD = @(t,x) prelimModel_param(t,x,tcutoff,init,paramsD,paramToChange);

[tA, solA] = ode45(sysA, tspan, init);
pA = basalSelfRenew(solA(:,3),paramsA);
[tB, solB] = ode45(sysB, tspan, init);
pB = basalSelfRenew(solB(:,3),paramsB);
[tC, solC] = ode45(sysC, tspan, init);
pC = basalSelfRenew(solC(:,3),paramsC);
[tD, solD] = ode45(sysD, tspan, init);
pD = basalSelfRenew(solD(:,3),paramsD);

fig3 = figure(3);
set(fig3,'Position',[220, 220, 900, 500])
fig3 = tiledlayout(1,2,'TileSpacing','compact');
nexttile(1)
hold on;
plot(tA,solA(:,1),'--','Color',"#77AC30",'DisplayName', ...
    ['Basal',num2str(paramsA.aE)])
plot(tA,solA(:,2),'-','Color',"#77AC30",'DisplayName', ...
    ['Lum',num2str(paramsA.aE)])
plot(tB,solB(:,1),'--','Color',"#EDB120",'DisplayName', ...
    ['Basal',num2str(paramsB.aE)])
plot(tB,solB(:,2),'-','Color',"#EDB120",'DisplayName', ...
    ['Lum',num2str(paramsB.aE)])
plot(tC,solC(:,1),'--','Color',"#0072BD",'DisplayName', ...
    ['Basal',num2str(paramsC.aE)])
plot(tC,solC(:,2),'-','Color',"#0072BD",'DisplayName', ...
    ['Lum',num2str(paramsC.aE)])
plot(tD,solD(:,1),'--','Color',"#7E2F8E",'DisplayName', ...
    ['Basal',num2str(paramsD.aE)])
plot(tD,solD(:,2),'-','Color',"#7E2F8E",'DisplayName', ...
    ['Lum',num2str(paramsD.aE)])
hold off;
title(fig3,'Changing $a_E$','Interpreter','latex')
lg = legend;
lg.NumColumns = 4;
lg.Layout.Tile = 'south';
title('A')
xlabel('Time')
ylabel('Cell Number')
nexttile(2)
hold on;
plot(tA,pA,'-','Color',"#77AC30")
plot(tB,pB,'-','Color',"#EDB120")
plot(tC,pC,'-','Color',"#0072BD")
plot(tD,pD,'-','Color',"#7E2F8E")
hold off;
title('B')
xlabel('Time')
ylabel('Basal Cell Self-renewal')
clear paramsA paramsB paramsC paramsD sysA sysB sysC sysD tA tB tC tD ...
    solA solB solC solD pA pB pC pD

%% k

paramToChange = 'k';
params.afterCutoff = 2;

paramsA = params;
paramsA.k = 1;
sysA = @(t,x) prelimModel_param(t,x,tcutoff,init,paramsA,paramToChange);

paramsB = params;
paramsB.k = 2;
sysB = @(t,x) prelimModel_param(t,x,tcutoff,init,paramsB,paramToChange);

paramsC = params;
paramsC.k = 3;
sysC = @(t,x) prelimModel_param(t,x,tcutoff,init,paramsC,paramToChange);

paramsD = params;
paramsD.k = 4;
sysD = @(t,x) prelimModel_param(t,x,tcutoff,init,paramsD,paramToChange);

[tA, solA] = ode45(sysA, tspan, init);
pA = basalSelfRenew(solA(:,3),paramsA);
[tB, solB] = ode45(sysB, tspan, init);
pB = basalSelfRenew(solB(:,3),paramsB);
[tC, solC] = ode45(sysC, tspan, init);
pC = basalSelfRenew(solC(:,3),paramsC);
[tD, solD] = ode45(sysD, tspan, init);
pD = basalSelfRenew(solD(:,3),paramsD);

fig4 = figure(4);
set(fig4,'Position',[230, 230, 900, 500])
fig4 = tiledlayout(1,2,'TileSpacing','compact');
nexttile(1)
hold on;
plot(tA,solA(:,1),'--','Color',"#77AC30",'DisplayName', ...
    ['Basal',num2str(paramsA.k)])
plot(tA,solA(:,2),'-','Color',"#77AC30",'DisplayName', ...
    ['Lum',num2str(paramsA.k)])
plot(tB,solB(:,1),'--','Color',"#EDB120",'DisplayName', ...
    ['Basal',num2str(paramsB.k)])
plot(tB,solB(:,2),'-','Color',"#EDB120",'DisplayName', ...
    ['Lum',num2str(paramsB.k)])
plot(tC,solC(:,1),'--','Color',"#0072BD",'DisplayName', ...
    ['Basal',num2str(paramsC.k)])
plot(tC,solC(:,2),'-','Color',"#0072BD",'DisplayName', ...
    ['Lum',num2str(paramsC.k)])
plot(tD,solD(:,1),'--','Color',"#7E2F8E",'DisplayName', ...
    ['Basal',num2str(paramsD.k)])
plot(tD,solD(:,2),'-','Color',"#7E2F8E",'DisplayName', ...
    ['Lum',num2str(paramsD.k)])
hold off;
title(fig4,'Changing $k$','Interpreter','latex')
lg = legend;
lg.NumColumns = 4;
lg.Layout.Tile = 'south';
title('A')
xlabel('Time')
ylabel('Cell Number')
nexttile(2)
hold on;
plot(tA,pA,'-','Color',"#77AC30")
plot(tB,pB,'-','Color',"#EDB120")
plot(tC,pC,'-','Color',"#0072BD")
plot(tD,pD,'-','Color',"#7E2F8E")
hold off;
title('B')
xlabel('Time')
ylabel('Basal Cell Self-renewal')
clear paramsA paramsB paramsC paramsD sysA sysB sysC sysD tA tB tC tD ...
    solA solB solC solD pA pB pC pD

%% Modified prelimModel

function dxdt = prelimModel_param(t,x,tcutoff,x0,params,paramToChange)
% PRELIMMODEL_PARAM implements the preliminary ODE model of bilayer duct 
% growth with chemomechanical feedback; flexible version where stimulation
% parameter can be varied
%   Inputs:
%       t (double): time at which to calculate dxdt
%       x (3x1 double): values for [B(t); L(t); S(t)]
%       tcutoff (double): time of switching paramToChange = afterCutoff
%       x0 (3x1 double): values for [B(0); L(0); S(0)]
%       params (1x1 struct): parameter values; has fields pbar, Abar,
%       rBbar, rLbar, rLDbar, rBDbar, gammaA, gammap, gammaB, gammaL,
%       tauSA, aE, k, afterCutoff
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

if t > tcutoff
    switch paramToChange
        case 'gammaL'
            gammaL = params.afterCutoff;
        case 'rLbar'
            rLbar = params.afterCutoff;
        case 'aE'
            aE = params.afterCutoff;
        case 'k'
            k = params.afterCutoff;
    end
end

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

function p = basalSelfRenew(S,params)
% BASALSELFRENEW computes the basal cell self-renewal probability from the
% stress history of the system
%   Inputs:
%       S (n x 1 double): vector of stress over time for system
%       params (1x1 struct): parameter values; has fields pbar, Abar,
%       rBbar, rLbar, rLDbar, rBDbar, gammaA, gammap, gammaB, gammaL,
%       tauSA, aE, k, afterCutoff
%   Outputs:
%       p (n x 1 double): vector of basal cell self-renewal probability
%       over time
sA = max(S - params.tauSA, zeros(size(S)));
A = params.Abar ./ (1 + params.gammaA * sA);
p = params.pbar ./ (1 + params.gammap * A);
end