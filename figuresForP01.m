% figuresForP01.m
% Brady Berg
clear; close all;
format long; format compact;
set(0,'defaultTextInterpreter','none')
set(0,'defaultAxesTickLabelInterpreter','none')
set(0,'defaultLegendInterpreter','none')
set(0,'defaultAxesFontName','Arial')
set(0,'defaultLineLineWidth',3)

yellow = "#EDB120"; % this is actually yellow
blue = "#0072BD";   % Basal
green = "#77AC30";  
orange = "#D95319"; % Luminal
% BLI are black (make thicker!)
BLIthickness = 4;
fontsize = 20;

%% Base parameters

params.pbar = 1;
params.Abar = 0.1;
params.rBbar = 1.2;
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

%% Homeostasis plot

paramsA1 = params;
paramsA1.rBbar = .6;
paramsA1.rLbar = .6;
initA1 = [20;20;0.0];
tspan = [0,1.2*2e2];
sysA1 = @(t,x) bilayerDuct_05(t,x,paramsA1,'nochange');
[tA1, solA1] = ode45(sysA1, tspan, initA1, opts);

fig1 = figure(1);
hold on
plot(tA1/1.2,solA1(:,2),'Color',orange,'DisplayName','Luminal')
plot(tA1/1.2,solA1(:,1),'Color',blue,'DisplayName','Basal','LineStyle','--')
plot(tA1/1.2,solA1(:,3),'Color','k','DisplayName','BLI','LineWidth',BLIthickness)
hold off
ylim([10^(-3),10^4])
set(gca,'YScale','log','FontSize',fontsize)
% xlabel('Time')
% ylabel('Cell number')
% exportgraphics(fig1,"C:\Users\bcber\OneDrive\Documents\MATLAB\LowengrubRotation\fig-homeostasis.pdf")

%% Demand production plot

paramsA1 = params;
paramsA1.rBbar = .6;
paramsA1.rLbar = .6;
initA1 = [20;20;0.0];
tspan = [0,1.2*5e1];
sysA1 = @(t,x) bilayerDuct_05(t,x,paramsA1,'nochange');
[tA1, solA1] = ode45(sysA1, tspan, initA1, opts);

% initA2 = [22.67;23.83;0.0043];
initA2 = solA1(end,:).';
paramsA2 = params;
paramsA2.rBbar = 1.2;
paramsA1.rLbar = 1.2;
tspan = [0,1.2*1e2];
sysA2 = @(t,x) bilayerDuct_05(t,x,paramsA2,'nochange');
[tA2, solA2] = ode45(sysA2, tspan, initA2, opts);

% initA3 = [66;68;0.105];
initA3 = solA2(end,:).';
paramsA3 = paramsA1;
tspan = [0,1.2*5e1];
[tA3, solA3] = ode45(sysA1, tspan, initA3, opts);

fig2 = figure(2);
hold on
plot(tA1/1.2,solA1(:,2),'DisplayName','Luminal','Color',orange)
plot(tA2/1.2 + 50,solA2(:,2),'DisplayName','Luminal','Color',orange)
plot(tA3/1.2 + 150,solA3(:,2),'DisplayName','Luminal','Color',orange)
plot(tA1/1.2,solA1(:,1),'DisplayName','Basal','Color',blue,'LineStyle','--')
plot(tA2/1.2 + 50,solA2(:,1),'DisplayName','Basal','Color',blue,'LineStyle','--')
plot(tA3/1.2 + 150,solA3(:,1),'DisplayName','Basal','Color',blue,'LineStyle','--')
plot(tA1/1.2,solA1(:,3),'DisplayName','BLI','Color','k','LineWidth',BLIthickness)
plot(tA2/1.2 + 50,solA2(:,3),'DisplayName','BLI','Color','k','LineWidth',BLIthickness)
plot(tA3/1.2 + 150,solA3(:,3),'DisplayName','BLI','Color','k','LineWidth',BLIthickness)
hold off
xlim('tight')
ylim([10^(-3),10^4])
set(gca,'YScale','log','FontSize',fontsize)
% exportgraphics(fig2,"C:\Users\bcber\OneDrive\Documents\MATLAB\LowengrubRotation\fig-demandprod.pdf")

%% Premalignant plot

initB = [20;20;8.911459];
% initB = [20;20;8];
tspan = [0,1.2*2e2];
paramsB = params;
paramsB.rLbar = 1.2;
paramsB.rBbar = 1.2;
sysB = @(t,x) bilayerDuct_05(t,x,paramsB,'nochange');
[tB, solB] = ode45(sysB, tspan, initB, opts);

fig3 = figure(3);
hold on
plot(tB/1.2,solB(:,2),'DisplayName','Luminal','Color',orange)
plot(tB/1.2,solB(:,1),'DisplayName','Basal','Color',blue,'LineStyle','--')
plot(tB/1.2,solB(:,3),'DisplayName','BLI','Color','k','LineWidth',BLIthickness)
hold off
ylim([10^(-3),10^4])
set(gca,'YScale','log','FontSize',fontsize)

% exportgraphics(fig3,"C:\Users\bcber\OneDrive\Documents\MATLAB\LowengrubRotation\fig-premalignant.pdf")

%% Malignant plot

initC = [20;20;10];
tspan = [0,1.2*2e2];
[tC, solC] = ode45(sysB, tspan, initC, opts);

fig4 = figure(4);
hold on
plot(tC/1.2,solC(:,2),'DisplayName','Luminal','Color',orange)
plot(tC/1.2,solC(:,1),'DisplayName','Basal','Color',blue,'LineStyle','--')
plot(tC/1.2,solC(:,3),'DisplayName','BLI','Color','k','LineWidth',BLIthickness)
hold off
ylim([10^(-3),10^4])
set(gca,'YScale','log','FontSize',fontsize)
% exportgraphics(fig4,"C:\Users\bcber\OneDrive\Documents\MATLAB\LowengrubRotation\fig-malignant.pdf")