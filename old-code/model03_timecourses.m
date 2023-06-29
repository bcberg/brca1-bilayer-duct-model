% model03_timecourses.m
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
orange = "#D95319";

%% Base parameters

pbar = 1;
Abar = 0.1;
rBbar = 1;
rLbar = 1.1;
rLDbar = 1;
rBDbar = 0.05;  % nonzero basal death rate
gammaA = 10;
gammap = 10;
gammaB = 40;
gammaL = 1;
tauSA = 0.1;
aE = 0.07;
k = 4;
tCutoff = Inf;
afterCutoff = [];
pPile = 0;
params = struct('pbar',pbar,'Abar',Abar,'rBbar',rBbar,'rLbar',rLbar,...
    'rLDbar',rLDbar,'rBDbar',rBDbar,'gammaA',gammaA,'gammap',gammap,...
    'gammaB',gammaB,'gammaL',gammaL,'tauSA',tauSA,'aE',aE,'k',k,...
    'tCutoff',tCutoff,'afterCutoff',afterCutoff,'pPile',pPile);
tspan = [0,100];
opts = odeset('RelTol',1e-5,'AbsTol',1e-6);

%% gammaL = 0.25

paramsA = params;
paramsA.pPile = 0.05;
paramsA.gammaL = 0.25;

sysA = @(t,x) bilayerDuct_03(t,x,paramsA,'nochange');
initA = [5;10;0];
[tA, solA] = ode45(sysA,tspan,initA,opts);
qA = usefulQuantities_03(solA,paramsA);

fig1 = figure(1);
set(fig1,'Position',[100,100,600,600])
fig1 = tiledlayout(2,2,'TileSpacing','compact');
title(fig1,'$\gamma_L = 0.25$', 'Interpreter', 'latex')
nexttile(1)
hold on;
plot(tA, solA(:,1), 'DisplayName', 'Basal')
plot(tA, solA(:,2), 'DisplayName', 'Luminal')
plot(tA, solA(:,3), 'DisplayName', '$L_{\rm pile}$')
hold off;
xlabel('Time');     ylabel('Cell number');
legend('Location', 'northwest')
nexttile(2)
hold on;
plot(tA, qA.SI, 'DisplayName', '$S_I$')
plot(tA, qA.SE, 'DisplayName', '$S_E$')
hold off;
xlabel('Time');     ylabel('Stress level');
legend('Location','east')
nexttile(3)
hold on;
plot(tA, qA.rB, 'DisplayName', '$r_B$')
plot(tA, qA.rL, 'DisplayName', '$r_L$')
hold off;
xlabel('Time');     ylabel('Proliferation rate');
legend('Location','east')
nexttile(4)
yyaxis left
plot(tA, qA.p, 'DisplayName', 'p')
ylabel('Basal self-renewal prob.')
yyaxis right
plot(tA, qA.A, 'DisplayName', 'A')
ylabel('Activin level')
xlabel('Time')
exportgraphics(fig1,'model03_time_gammaL025.tiff')

%% gammaL = 1

paramsB = params;
paramsB.pPile = 0.05;
paramsB.gammaL = 1;
tspan = [0,60];
sysB = @(t,x) bilayerDuct_03(t,x,paramsB,'nochange');
initB = [5;10;0];
[tB, solB] = ode45(sysB,tspan,initB,opts);
qB = usefulQuantities_03(solB,paramsB);

fig2 = figure(2);
set(fig2,'Position',[100,100,600,600])
fig2 = tiledlayout(2,2,'TileSpacing','compact');
title(fig2,'$\gamma_L = 1$', 'Interpreter', 'latex')
nexttile(1)
hold on;
plot(tB, solB(:,1), 'DisplayName', 'Basal')
plot(tB, solB(:,2), 'DisplayName', 'Luminal')
plot(tB, solB(:,3), 'DisplayName', '$L_{\rm pile}$')
hold off;
xlabel('Time');     ylabel('Cell number');
legend('Location', 'northwest')
nexttile(2)
hold on;
plot(tB, qB.SI, 'DisplayName', '$S_I$')
plot(tB, qB.SE, 'DisplayName', '$S_E$')
hold off;
xlabel('Time');     ylabel('Stress level');
legend('Location','east')
nexttile(3)
hold on;
plot(tB, qB.rB, 'DisplayName', '$r_B$')
plot(tB, qB.rL, 'DisplayName', '$r_L$')
hold off;
xlabel('Time');     ylabel('Proliferation rate');
ylim([0,1.2]) % to match t200 plot
legend('Location','east')
nexttile(4)
yyaxis left
plot(tB, qB.p, 'DisplayName', 'p')
ylabel('Basal self-renewal prob.')
yyaxis right
plot(tB, qB.A, 'DisplayName', 'A')
ylabel('Activin level')
xlabel('Time')
exportgraphics(fig2,'model03_time_gammaL1_t60.tiff')

% run for longer time
tspan = [0,200];
sysB = @(t,x) bilayerDuct_03(t,x,paramsB,'nochange');
initB = [5;10;0];
[tB, solB] = ode45(sysB,tspan,initB,opts);
qB = usefulQuantities_03(solB,paramsB);

fig3 = figure(3);
set(fig3,'Position',[100,100,600,600])
fig3 = tiledlayout(2,2,'TileSpacing','compact');
title(fig3,'$\gamma_L = 1$', 'Interpreter', 'latex')
nexttile(1)
hold on;
plot(tB, solB(:,1), 'DisplayName', 'Basal')
plot(tB, solB(:,2), 'DisplayName', 'Luminal')
plot(tB, solB(:,3), 'DisplayName', '$L_{\rm pile}$')
hold off;
xlabel('Time');     ylabel('Cell number');
legend('Location', 'northwest')
nexttile(2)
hold on;
plot(tB, qB.SI, 'DisplayName', '$S_I$')
plot(tB, qB.SE, 'DisplayName', '$S_E$')
hold off;
xlabel('Time');     ylabel('Stress level');
legend('Location','east')
nexttile(3)
hold on;
plot(tB, qB.rB, 'DisplayName', '$r_B$')
plot(tB, qB.rL, 'DisplayName', '$r_L$')
hold off;
xlabel('Time');     ylabel('Proliferation rate');
legend('Location','east')
nexttile(4)
yyaxis left
plot(tB, qB.p, 'DisplayName', 'p')
ylabel('Basal self-renewal prob.')
yyaxis right
plot(tB, qB.A, 'DisplayName', 'A')
ylabel('Activin level')
xlabel('Time')
exportgraphics(fig3,'model03_time_gammaL1_t200.tiff')