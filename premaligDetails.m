% premaligDetails.m
% Brady Berg
clear; close all;
format long; format compact;
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')
set(0,'defaultLineLineWidth',0.8)
blue = "#0072BD";   % Basal
orange = "#D95319"; % Luminal
% BLI curves in black
%% Base parameters

params.pbar = 1;
params.Abar = 0.1;
params.rBbar = 1.2; % higher drive for growth
params.rLbar = 1.2; % higher drive for growth
params.rBDbar = 0.1;
params.rLDbar = 0.25;
params.rPDbar = 1;
params.gammap = 10;
params.gammaA = 20;
params.gammaB = 20;
params.gammaL = 20;  %10
params.pPile = 1e-4;    %1e-4 % 1e-3
params.gammaPile = 10;  % 20
params.gammaAP = 1;
params.gammaFST = 2;
params.nA = 1;
params.nFST = 2;
params.tauSA = 0.1;
params.aE = 0.05;     %0.075;
params.k = 5;

% optional parameter step-discontinuity
params.tCutoff = Inf;
params.afterCutoff = [];

opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'Events',@detectOvergrowth);
% opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

%% Plotting solution detail

sysA = @(t,x) bilayerDuct_05(t,x,params,'nochange');
tspan = [0,1.2*2e2];
initA = [20;20;8.911459];
% initA = [66.252600602161010; 68.252787226175542; 8];
labelA = ['$[B(0);L(0);BLI(0)] = ', mat2str(initA), '$'];
[tA, solA] = ode45(sysA, tspan, initA, opts);
tAscaled = tA / 1.2;
qA = usefulQuantities_05(solA, params);

% finding when stress threshold is crossed
% SA starts at identically 0, but becomes positive at some value and then
% zero again at some value
SA_A_critind = find(diff(sign(qA.SA)));

fig1 = figure(1);
set(fig1,'Position',[100,100,600,600])
fig1 = tiledlayout(2,2,'TileSpacing','compact');
title(fig1, labelA, 'Interpreter', 'latex')
nexttile(1)
hold on;
plot(tAscaled, solA(:,2), 'DisplayName', 'Luminal', 'LineStyle', '-', ...
    'Color', orange)
plot(tAscaled, solA(:,1), 'DisplayName', 'Basal', 'LineStyle', '--', ...
    'Color', blue)
plot(tAscaled, solA(:,3), 'DisplayName', 'BLI', 'LineStyle', '-', ...
    'Color', 'k')
hold off;
% set(gca, 'YScale', 'log')
xlabel('Time');     ylabel('Cell number');
legend('Location','east')
nexttile(2)
hold on;
plot(tAscaled, qA.SI, 'DisplayName', '$S_I$')
plot(tAscaled(SA_A_critind), qA.SI(SA_A_critind), 'ob', 'DisplayName', ...
    '$S_I = \tau_{SA}$')
hold off;
xlabel('Time');     ylabel('Stress level');
legend('Location','east')
nexttile(3)
hold on;
plot(tAscaled, qA.rL, 'DisplayName', '$r_L$', 'LineStyle', '-', ...
    'Color', orange)
plot(tAscaled, qA.rB, 'DisplayName', '$r_B$', 'LineStyle', '--', ...
    'Color', blue)
plot(tAscaled, qA.rPile, 'DisplayName', '$r_{BLI}$', 'LineStyle', '-', ...
    'Color', 'k')
hold off;
xlabel('Time');     ylabel('Proliferation rate');
legend('Location','east')
nexttile(4)
yyaxis left
plot(tAscaled, qA.p, 'DisplayName', '$p$')
ylabel('Basal self-renewal prob.')
yyaxis right
plot(tAscaled, qA.A, 'DisplayName', '$A$')
ylabel('Activin level')
xlabel('Time')

% plotting solution and first derivative
dBdt_A = ((2 * qA.p - 1) .* qA.rB - params.rBDbar - params.pPile) .* solA(:,1);
dBdt_A_critind = find(diff(sign(dBdt_A)));
dLdt_A = 2 * (1 - qA.p) .* qA.rB .* solA(:,1) + (qA.rL - params.rLDbar) .* solA(:,2);
dLdt_A_critind = find(diff(sign(dLdt_A)));
dLpiledt_A = params.pPile * solA(:,1) + (qA.rPile - params.rPDbar) .* solA(:,3);
dLpiledt_A_critind = find(diff(sign(dLpiledt_A)));

fig2 = figure(2);
set(fig2,'Position',[110,110,600,600])
fig2 = tiledlayout(2,1,'TileSpacing','compact');
title(fig2,labelA,'Interpreter','latex')
nexttile(1)
hold on;
l1 = plot(tAscaled, solA(:,2), 'DisplayName', 'Luminal', 'LineStyle', '-', ...
    'Color', orange);
plot(tAscaled(dLdt_A_critind), solA(dLdt_A_critind,2), 'or')
l2 = plot(tAscaled, solA(:,1), 'DisplayName', 'Basal', 'LineStyle', '--', ...
    'Color', blue);
plot(tAscaled(dBdt_A_critind), solA(dBdt_A_critind,1), 'ob')
l3 = plot(tAscaled, solA(:,3), 'DisplayName', 'BLI', 'LineStyle', '-', ...
    'Color', 'k');
plot(tAscaled(dLpiledt_A_critind), solA(dLpiledt_A_critind,3), 'ok')
hold off;
xlabel('Time');     ylabel('Cell number');
legend([l1,l2,l3],'Location','southeast')
nexttile(2)
hold on;
l4 = plot(tAscaled, dLdt_A, 'DisplayName', '$dL/dt$', 'LineStyle', '-', ...
    'Color', orange);
plot(tAscaled(dLdt_A_critind), dLdt_A(dLdt_A_critind), 'or')
l5 = plot(tAscaled, dBdt_A, 'DisplayName', '$dB/dt$', 'LineStyle', '--', ...
    'Color', blue);
plot(tAscaled(dBdt_A_critind), dBdt_A(dBdt_A_critind), 'ob')
l6 = plot(tAscaled, dLpiledt_A, 'DisplayName', '$d(BLI)/dt$', 'LineStyle', '-', ...
    'Color', 'k');
plot(tAscaled(dLpiledt_A_critind), dLpiledt_A(dLpiledt_A_critind), 'ok')
hold off;
ylim('tight')
xlabel('Time');     ylabel('Time derivatives');
legend([l4,l5,l6],'Location','northeast')

%% Other initial conditions

% params.tauSA = 0.1;
sysB = @(t,x) bilayerDuct_05(t,x,params,'nochange');

% initB = [20;20;8.911459];
% % initB = [solA(end,1:2).';initA(end)];
% initB = [20;20;8.65];
% initB = [66.252600602161010; 68.252787226175542; 8];
initB = [20;20;8];
labelB = ['$[B(0);L(0);BLI(0)] = ', mat2str(initB), '$'];
[tB, solB] = ode45(sysB, tspan, initB, opts);
tBscaled = tB / 1.2;
qB = usefulQuantities_05(solB, params);

% finding when stress threshold is crossed
% SA starts at identically 0, but becomes positive at some value and then
% zero again at some value
SA_B_critind = find(diff(sign(qB.SA)));

fig3 = figure(3);
set(fig3,'Position',[120,120,600,600])
fig3 = tiledlayout(2,2,'TileSpacing','compact');
title(fig3, labelB, 'Interpreter', 'latex')
nexttile(1)
hold on;
plot(tBscaled, solB(:,2), 'DisplayName', 'Luminal', 'LineStyle', '-', ...
    'Color', orange)
plot(tBscaled, solB(:,1), 'DisplayName', 'Basal', 'LineStyle', '--', ...
    'Color', blue)
plot(tBscaled, solB(:,3), 'DisplayName', 'BLI', 'LineStyle', '-', ...
    'Color', 'k')
hold off;
% set(gca, 'YScale', 'log')
xlabel('Time');     ylabel('Cell number');
legend('Location','east')
nexttile(2)
hold on;
plot(tBscaled, qB.SI, 'DisplayName', '$S_I$')
plot(tBscaled(SA_B_critind), qB.SI(SA_B_critind), 'ob', 'DisplayName', ...
    '$S_I = \tau_{SA}$')
hold off;
xlabel('Time');     ylabel('Stress level');
legend('Location','east')
nexttile(3)
hold on;
plot(tBscaled, qB.rL, 'DisplayName', '$r_L$', 'LineStyle', '-', ...
    'Color', orange)
plot(tBscaled, qB.rB, 'DisplayName', '$r_B$', 'LineStyle', '--', ...
    'Color', blue)
plot(tBscaled, qB.rPile, 'DisplayName', '$r_{BLI}$', 'LineStyle', '-', ...
    'Color', 'k')
hold off;
xlabel('Time');     ylabel('Proliferation rate');
legend('Location','east')
nexttile(4)
yyaxis left
plot(tBscaled, qB.p, 'DisplayName', '$p$')
ylabel('Basal self-renewal prob.')
yyaxis right
plot(tBscaled, qB.A, 'DisplayName', '$A$')
ylabel('Activin level')
xlabel('Time')

% plotting solution and first derivative
dBdt_B = ((2 * qB.p - 1) .* qB.rB - params.rBDbar - params.pPile) .* solB(:,1);
dBdt_B_critind = find(diff(sign(dBdt_B)));
dLdt_B = 2 * (1 - qB.p) .* qB.rB .* solB(:,1) + (qB.rL - params.rLDbar) .* solB(:,2);
dLdt_B_critind = find(diff(sign(dLdt_B)));
dLpiledt_B = params.pPile * solB(:,1) + (qB.rPile - params.rPDbar) .* solB(:,3);
dLpiledt_B_critind = find(diff(sign(dLpiledt_B)));

fig4 = figure(4);
set(fig4,'Position',[130,130,600,600])
fig4 = tiledlayout(2,1,'TileSpacing','compact');
title(fig4,labelB,'Interpreter','latex')
nexttile(1)
hold on;
l1 = plot(tBscaled, solB(:,2), 'DisplayName', 'Luminal', 'LineStyle', '-', ...
    'Color', orange);
plot(tBscaled(dLdt_B_critind), solB(dLdt_B_critind,2), 'or')
l2 = plot(tBscaled, solB(:,1), 'DisplayName', 'Basal', 'LineStyle', '--', ...
    'Color', blue);
plot(tBscaled(dBdt_B_critind), solB(dBdt_B_critind,1), 'ob')
l3 = plot(tBscaled, solB(:,3), 'DisplayName', 'BLI', 'LineStyle', '-', ...
    'Color', 'k');
plot(tBscaled(dLpiledt_B_critind), solB(dLpiledt_B_critind,3), 'ok')
hold off;
xlabel('Time');     ylabel('Cell number');
legend([l1,l2,l3],'Location','southeast')
nexttile(2)
hold on;
l4 = plot(tBscaled, dLdt_B, 'DisplayName', '$dL/dt$', 'LineStyle', '-', ...
    'Color', orange);
plot(tBscaled(dLdt_B_critind), dLdt_B(dLdt_B_critind), 'or')
l5 = plot(tBscaled, dBdt_B, 'DisplayName', '$dB/dt$', 'LineStyle', '--', ...
    'Color', blue);
plot(tBscaled(dBdt_B_critind), dBdt_B(dBdt_B_critind), 'ob')
l6 = plot(tBscaled, dLpiledt_B, 'DisplayName', '$d(BLI)/dt$', 'LineStyle', '-', ...
    'Color', 'k');
plot(tBscaled(dLpiledt_B_critind), dLpiledt_B(dLpiledt_B_critind), 'ok')
hold off;
ylim('tight')
xlabel('Time');     ylabel('Time derivatives');
legend([l4,l5,l6],'Location','northeast')

%% Functions

function [value, isterminal, direction] = detectOvergrowth(~,y)
value = y(3) - 1e9;     % when value = 0 (Lpile reaches 1e9), ...
isterminal = 1;         % terminate integration
direction = 0;
end