% model02.m
% Brady Berg
clear; close all;
format long; format compact;
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')

green = "#77AC30";
yellow = "#EDB120";
blue = "#0072BD";
purple = "#7E2F8E";

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
tCutoff = 40;
afterCutoff = [];
params = struct('pbar',pbar,'Abar',Abar,'rBbar',rBbar,'rLbar',rLbar,...
    'rLDbar',rLDbar,'rBDbar',rBDbar,'gammaA',gammaA,'gammap',gammap,...
    'gammaB',gammaB,'gammaL',gammaL,'tauSA',tauSA,'aE',aE,'k',k,...
    'tCutoff',tCutoff,'afterCutoff',afterCutoff);
tspan = [0,2*tCutoff];
init = [1;1];

%% Assuming elastic equilibrium, alter gammaL

paramToChange = 'gammaL';
params.afterCutoff = 2.5;

paramsA = params;
paramsA.gammaL = 0.25;
sysA = @(t,x) bilayerDuct_02(t,x,init,paramsA,paramToChange);

paramsB = params;
paramsB.gammaL = 1;
sysB = @(t,x) bilayerDuct_02(t,x,init,paramsB,paramToChange);

paramsC = params;
paramsC.gammaL = 2.5;
sysC = @(t,x) bilayerDuct_02(t,x,init,paramsC,paramToChange);

opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
[tA, solA] = ode45(sysA, tspan, init, opts);
qA = usefulQuantities_02(solA,paramsA);
[tB, solB] = ode45(sysB, tspan, init, opts);
qB = usefulQuantities_02(solB,paramsB);
[tC, solC] = ode45(sysC, tspan, init, opts);
qC = usefulQuantities_02(solC,paramsC);
cutoffIdxs = [find(tA>paramsA.tCutoff,1);find(tB>paramsB.tCutoff,1)];

fig1 = figure(1);
set(fig1,'Position',[200, 200, 500, 500])
fig1 = tiledlayout(2,2,'TileSpacing','compact');
nexttile(1)
hold on;
l1 = plot(tA(1:cutoffIdxs(1)), solA(1:cutoffIdxs(1),1), '--', ...
    'Color',green,'DisplayName','Basal0.25');
plot(tA(cutoffIdxs(1):end), solA(cutoffIdxs(1):end,1), '--', ...
    'Color',blue)
l2 = plot(tA(1:cutoffIdxs(1)), solA(1:cutoffIdxs(1),2), '-', ...
    'Color',green,'DisplayName','Lum0.25');
plot(tA(cutoffIdxs(1):end), solA(cutoffIdxs(1):end,2), '-', ...
    'Color',blue)
l3 = plot(tB(1:cutoffIdxs(2)), solB(1:cutoffIdxs(2),1), '--', ...
    'Color',yellow,'DisplayName','Basal1');
plot(tB(cutoffIdxs(2):end), solB(cutoffIdxs(2):end,1), '--', ...
    'Color',blue)
l4 = plot(tB(1:cutoffIdxs(2)), solB(1:cutoffIdxs(2),2), '-', ...
    'Color',yellow,'DisplayName','Lum1');
plot(tB(cutoffIdxs(2):end), solB(cutoffIdxs(2):end,2), '-', ...
    'Color',blue)
l5 = plot(tC,solC(:,1),'--','Color',blue,'DisplayName','Basal2.5');
l6 = plot(tC,solC(:,2),'-','Color',blue,'DisplayName','Lum2.5');
hold off;
title('\bf A')
xlabel('Time')
ylabel('Cell Number')
xlim(tspan)
nexttile(2)
hold on;
plot(tA,qA.p,'-','Color',green)
plot(tB,qB.p,'-','Color',yellow)
plot(tC,qC.p,'-','Color',blue)
hold off;
title('\bf B')
xlabel('Time')
ylabel('Basal Cell Self-renewal')
xlim(tspan)
nexttile(3)
hold on;
plot(tA,qA.rL,'-','Color',green)
plot(tB,qB.rL,'-','Color',yellow)
plot(tC,qC.rL,'-','Color',blue)
hold off;
title('\bf C')
xlabel('Time')
ylabel('Luminal Cell Proliferation')
nexttile(4)
hold on;
plot(tA,qA.rB,'-','Color',green)
plot(tB,qB.rB,'-','Color',yellow)
plot(tC,qC.rB,'-','Color',blue)
hold off;
title('\bf D')
xlabel('Time')
ylabel('Basal Cell Proliferation')
title(fig1,'Changing $\gamma_L$','Interpreter','latex')
lg = legend([l1,l2,l3,l4,l5,l6]);
lg.NumColumns = 3;
lg.Layout.Tile = 'south';

%% Assuming elastic equilibrium, alter rLbar

paramToChange = 'rLbar';
params.afterCutoff = 0.9;

paramsA = params;
paramsA.rLbar = 1.5;
sysA = @(t,x) bilayerDuct_02(t,x,init,paramsA,paramToChange);

paramsB = params;
paramsB.rLbar = 1.2;
sysB = @(t,x) bilayerDuct_02(t,x,init,paramsB,paramToChange);

paramsC = params;
paramsC.rLbar = 0.9;
sysC = @(t,x) bilayerDuct_02(t,x,init,paramsC,paramToChange);

opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
[tA, solA] = ode45(sysA, tspan, init, opts);
qA = usefulQuantities_02(solA,paramsA);
[tB, solB] = ode45(sysB, tspan, init, opts);
qB = usefulQuantities_02(solB,paramsB);
[tC, solC] = ode45(sysC, tspan, init, opts);
qC = usefulQuantities_02(solC,paramsC);
cutoffIdxs = [find(tA>paramsA.tCutoff,1);find(tB>paramsB.tCutoff,1)];

fig2 = figure(2);
set(fig2,'Position',[210, 210, 500, 500])
fig2 = tiledlayout(2,2,'TileSpacing','compact');
nexttile(1)
hold on;
l1 = plot(tA(1:cutoffIdxs(1)), solA(1:cutoffIdxs(1),1), '--', ...
    'Color',green,'DisplayName',['Basal',num2str(paramsA.rLbar)]);
plot(tA(cutoffIdxs(1):end), solA(cutoffIdxs(1):end,1), '--', ...
    'Color',blue)
l2 = plot(tA(1:cutoffIdxs(1)), solA(1:cutoffIdxs(1),2), '-', ...
    'Color',green,'DisplayName',['Lum',num2str(paramsA.rLbar)]);
plot(tA(cutoffIdxs(1):end), solA(cutoffIdxs(1):end,2), '-', ...
    'Color',blue)
l3 = plot(tB(1:cutoffIdxs(2)), solB(1:cutoffIdxs(2),1), '--', ...
    'Color',yellow,'DisplayName',['Basal',num2str(paramsB.rLbar)]);
plot(tB(cutoffIdxs(2):end), solB(cutoffIdxs(2):end,1), '--', ...
    'Color',blue)
l4 = plot(tB(1:cutoffIdxs(2)), solB(1:cutoffIdxs(2),2), '-', ...
    'Color',yellow,'DisplayName',['Lum',num2str(paramsB.rLbar)]);
plot(tB(cutoffIdxs(2):end), solB(cutoffIdxs(2):end,2), '-', ...
    'Color',blue)
l5 = plot(tC,solC(:,1),'--','Color',blue,'DisplayName', ...
    ['Basal',num2str(paramsC.rLbar)]);
l6 = plot(tC,solC(:,2),'-','Color',blue,'DisplayName', ...
    ['Lum',num2str(paramsC.rLbar)]);
hold off;
title('\bf A')
xlabel('Time')
ylabel('Cell Number')
xlim(tspan)
nexttile(2)
hold on;
plot(tA,qA.p,'-','Color',green)
plot(tB,qB.p,'-','Color',yellow)
plot(tC,qC.p,'-','Color',blue)
hold off;
title('\bf B')
xlabel('Time')
ylabel('Basal Cell Self-renewal')
xlim(tspan)
nexttile(3)
hold on;
plot(tA,qA.rL,'-','Color',green)
plot(tB,qB.rL,'-','Color',yellow)
plot(tC,qC.rL,'-','Color',blue)
hold off;
title('\bf C')
xlabel('Time')
ylabel('Luminal Cell Proliferation')
nexttile(4)
hold on;
plot(tA,qA.rB,'-','Color',green)
plot(tB,qB.rB,'-','Color',yellow)
plot(tC,qC.rB,'-','Color',blue)
hold off;
title('\bf D')
xlabel('Time')
ylabel('Basal Cell Proliferation')
title(fig2,'Changing $\bar{r}_L$','Interpreter','latex')
lg = legend([l1,l2,l3,l4,l5,l6]);
lg.NumColumns = 3;
lg.Layout.Tile = 'south';

%% Functions (now in separate files)

% function dxdt = bilayerDuct_02(t,x,x0,params,paramToChange)
% % BILAYERDUCT_02 implements an ODE model of bilayer duct 
% % growth with chemomechanical feedback, assuming stress reaches equilibrium
% % faster than the growth dynamics
% %   Inputs:
% %       t (double): time at which to calculate dxdt
% %       x (2x1 double): values for [B(t); L(t)]
% %       x0 (2x1 double): values for [B(0); L(0)]
% %       params (1x1 struct): parameter values; has fields pbar, Abar,
% %       rBbar, rLbar, rLDbar, rBDbar, gammaA, gammap, gammaB, gammaL,
% %       tauSA, aE, k, tcutoff, afterCutoff
% %       paramToChange (string): name of parameter to change to the value
% %       afterCutoff at time tcutoff
% %   Outputs:
% %       dxdt (2x1 double): time derivative of x
% pbar = params.pbar;
% Abar = params.Abar;
% rBbar = params.rBbar;
% rLbar = params.rLbar;
% rLDbar = params.rLDbar;
% rBDbar = params.rBDbar;
% gammaA = params.gammaA;
% gammap = params.gammap;
% gammaB = params.gammaB;
% gammaL = params.gammaL;
% tauSA = params.tauSA;
% aE = params.aE;
% k = params.k;
% 
% if t > params.tCutoff
%     switch paramToChange
%         case 'gammaL'
%             gammaL = params.afterCutoff;
%         case 'rLbar'
%             rLbar = params.afterCutoff;
%         case 'aE'
%             aE = params.afterCutoff;
%         case 'k'
%             k = params.afterCutoff;
%     end
% end
% 
% B = x(1);
% L = x(2);
% S = (L - B) / k;
% 
% SA = max(S - tauSA, 0);
% SE = (B + L) / (x0(1) + x0(2));
% A = Abar / (1 + gammaA * SA);
% p = pbar / (1 + gammap * A);
% rL = rLbar / ( (1 + gammaL * S) * (1 + aE * SE) );
% rB = rBbar / ( (1 + gammaB * A) * (1 + aE * SE) );
% dxdt = [(2*p - 1) * rB * B - rBDbar * B;
%     2*(1-p) * rB * B + (rL - rLDbar) * L];
% end
% 
% function output = usefulQuantities_02(sol,params)
% % USEFULQUANTITIES_02 computes quantities that help describe system 
% % behavior from the state history, assuming stress equilibrates faster
% % than the growth dynamics
% %   Inputs:
% %       sol (n x 2 double): vector of cell quantities over time for system
% %       params (1x1 struct): parameter values; has fields pbar, Abar,
% %       rBbar, rLbar, rLDbar, rBDbar, gammaA, gammap, gammaB, gammaL,
% %       tauSA, aE, k, tCutoff, afterCutoff
% %   Outputs:
% %       output (1x1 struct): time courses; has fields S, SA, A, p, SE, rL,
% %       rB
% output.S = (sol(:,2) - sol(:,1)) / params.k;
% output.SA = max(output.S - params.tauSA, zeros(size(output.S)));
% output.A = params.Abar ./ (1 + params.gammaA * output.SA);
% output.p = params.pbar ./ (1 + params.gammap * output.A);
% output.SE = sum(sol,2) ./ (sol(1,1) + sol(1,2));
% output.rL = params.rLbar ./ ( (1 + params.gammaL * output.S).*(1 + ...
%     params.aE * output.SE) );
% output.rB = params.rBbar ./ ( (1 + params.gammaB * output.A).*(1 + ...
%     params.aE * output.SE) );
% end