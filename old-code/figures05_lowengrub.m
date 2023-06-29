% figures05.m
% Brady Berg
clear; close all;
format long; format compact;
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')

green = "#078A0D";
red = "#CC0000";
yellow = "#EDB120";

%% Parameters

params.pbar = 1;
params.Abar = 0.1;
params.rBbar = 1;
params.rLbar = 1.2; % 1.1, 1.3, 1.5
params.rBDbar = 0.1;    %0.1
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

% opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'Events',@detectOvergrowth);
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

%% Low steady regime

paramsA = params;
paramsA.rBbar = 1.2;
paramsA.rLbar = 1.2;
initA1 = [20;20;0];
%initA1 = [66;68;0];
tspan = [0,2e2];
sysA = @(t,x) bilayerDuct_05(t,x,paramsA,'nochange');
[tA1, solA1] = ode45(sysA, tspan, initA1, opts);

fig1 = figure(1);
tiledlayout('flow')
nexttile
hold on
set(gca,'FontSize',24)
set(gca,'DefaultLineLineWidth', 2)
plot(tA1,solA1(:,1),'DisplayName','Basal','Color',green)
plot(tA1,solA1(:,2),'DisplayName','Luminal','Color',red)
plot(tA1,solA1(:,3),'DisplayName','BLI','Color',yellow)
hold off
xlabel('Time (arb. units)')
ylabel('Cell number (arb. units)')
ylim([0,100])
lg1 = legend;
lg1.Layout.Tile = 'east';

%exportgraphics(fig1,"C:\Users\bcber\OneDrive\Documents\1-UCI\Q3\Lowengrub Rotation\fig_noBLI-draft.pdf")

%% Med quasi-steady regime

%initA2 = [20;20;9];
%initA2 = [20;20;6.731];
initA2 = [20;20;8.911459];
tspan = [0,2e2];
sysA = @(t,x) bilayerDuct_05(t,x,paramsA,'nochange');
[tA2, solA2] = ode45(sysA, tspan, initA2, opts);

fig2 = figure(2);
tiledlayout('flow')
nexttile
hold on
set(gca,'FontSize',24)
set(gca,'DefaultLineLineWidth', 2)
plot(tA2,solA2(:,1),'DisplayName','Basal','Color',green)
plot(tA2,solA2(:,2),'DisplayName','Luminal','Color',red)
plot(tA2,solA2(:,3),'DisplayName','BLI','Color',yellow)
hold off
xlabel('Time (arb. units)')
ylabel('Cell number (arb. units)')
ylim([0,100])
lg2 = legend;
lg2.Layout.Tile = 'east';

%% Destabilizing the quasi-steady

paramsB = paramsA;
paramsB.rLbar = 1.200;
initA3 = [20;20;15];
tspan = [0,2e2];
sysB = @(t,x) bilayerDuct_05(t,x,paramsB,'nochange');
[tB, solB] = ode45(sysB, tspan, initA3, opts);

fig3 = figure(3);
tiledlayout('flow')
nexttile
hold on
set(gca,'FontSize',24)
set(gca,'DefaultLineLineWidth', 2)
plot(tB,solB(:,1),'DisplayName','Basal','Color',green)
plot(tB,solB(:,2),'DisplayName','Luminal','Color',red)
plot(tB,solB(:,3),'DisplayName','BLI','Color',yellow)
hold off
xlabel('Time (arb. units)')
ylabel('Cell number (arb. units)')
ylim([0,250])
lg3 = legend;
lg3.Layout.Tile = 'east';


%% Destabilizing the quasi-steady with larger rL

paramsB = paramsA;
paramsB.rLbar = 1.300;
initA4 = [20;20;8.911459];
tspan = [0,2e2];
sysB = @(t,x) bilayerDuct_05(t,x,paramsB,'nochange');
[tB1, solB1] = ode45(sysB, tspan, initA4, opts);

fig4 = figure(4);
tiledlayout('flow')
nexttile
hold on
set(gca,'FontSize',24)
set(gca,'DefaultLineLineWidth', 2)
plot(tB1,solB1(:,1),'DisplayName','Basal','Color',green)
plot(tB1,solB1(:,2),'DisplayName','Luminal','Color',red)
plot(tB1,solB1(:,3),'DisplayName','BLI','Color',yellow)
hold off
xlabel('Time (arb. units)')
ylabel('Cell number (arb. units)')
ylim([0,250])
lg3 = legend;
lg3.Layout.Tile = 'east';


%% Showing movement of separatrix

fig5 = figure(5);
tiledlayout('flow')
nexttile
hold on
set(gca,'FontSize',16)
set(gca,'DefaultLineLineWidth', 2)
plot(tA1,solA1(:,1)+solA1(:,2)+solA1(:,3),'-b','DisplayName','$B+L$, stable', ...
    'LineWidth',2)
% plot(tA1,solA1(:,3),'--k','DisplayName','$B+L$, stable', ...
%     'LineWidth',1)
plot(tA2,solA2(:,1)+solA2(:,2)+solA2(:,3),'-r','DisplayName','$B+L$, stable', ...
    'LineWidth',2)
% plot(tA2,solA2(:,3),'--','DisplayName','BLI, stable','Color',yellow, ...
%     'LineWidth',1)
plot(tB,solB(:,1)+solB(:,2)+solB(:,3),'-k','DisplayName','$B+L$, unstable', ...
    'LineWidth',2)
plot(tB1,solB1(:,1)+solB1(:,2)+solB1(:,3),'--r','DisplayName','$B+L$, unstable', ...
    'LineWidth',2)
% plot(tB,solB(:,3),'-','DisplayName','BLI, unstable','Color',yellow, ...
%     'LineWidth',0.75)
hold off
xlabel('Time (arb. units)')
ylabel('Total Cell number (arb. units)')
ylim([0,400])
lg4 = legend;
lg4.Layout.Tile = 'south';
lg4.NumColumns = 2;
%exportgraphics(fig4,"C:\Users\bcber\OneDrive\Documents\1-UCI\Q3\Lowengrub Rotation\fig_withBLI-draft.pdf")
%% Showing movement of separatrix 2

fig6 = figure(6);
tiledlayout('flow')
nexttile
hold on
set(gca,'FontSize',24)
set(gca,'DefaultLineLineWidth', 2)
plot(tA1,solA1(:,3),'Color',yellow,'DisplayName','$B+L$, stable', ...
    'LineWidth',2)
% plot(tA1,solA1(:,3),'--k','DisplayName','$B+L$, stable', ...
%     'LineWidth',1)
plot(tA2,solA2(:,3),'Color',yellow,'DisplayName','$B+L$, stable', ...
    'LineWidth',2)
% plot(tA2,solA2(:,3),'--','DisplayName','BLI, stable','Color',yellow, ...
%     'LineWidth',1)
plot(tB,solB(:,3),'Color',yellow,'DisplayName','$B+L$, unstable', ...
    'LineWidth',2)
plot(tB1,solB1(:,3),'Color',yellow,'DisplayName','$B+L$, unstable', ...
    'LineWidth',2)
% plot(tB,solB(:,3),'-','DisplayName','BLI, unstable','Color',yellow, ...
%     'LineWidth',0.75)
hold off
xlabel('Time (arb. units)')
ylabel('Cell number (arb. units)')
xlim([0,100])
ylim([0,100])
lg4 = legend;
lg4.Layout.Tile = 'south';
lg4.NumColumns = 2;
%exportgraphics(fig4,"C:\Users\bcber\OneDrive\Documents\1-UCI\Q3\Lowengrub Rotation\fig_withBLI-draft.pdf")

%% Functions

function [value, isterminal, direction] = detectOvergrowth(~,y)
value = y(3) - 1e9;     % when value = 0 (Lpile reaches 1e9), ...
isterminal = 1;         % terminate integration
direction = 0;
end