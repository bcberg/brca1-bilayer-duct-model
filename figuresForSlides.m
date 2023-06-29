% figuresForSlides.m
% Brady Berg
clear; close all;
format long; format compact;
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')
% set(0,'defaultAxesFontName','Arial')
set(0,'defaultLineLineWidth',2)

yellow = "#EDB120";
blue = "#0072BD";   % Basal
green = "#77AC30";  
orange = "#D95319"; % Luminal
% BLI are black (make thicker!)
BLIthickness = 2.5;
fontsize = 16;
C = colororder('default');

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
set(fig1,'Position',[1000,818,680,420])
hold on
plot(tA1/1.2,solA1(:,2),'Color',orange,'DisplayName','Luminal')
plot(tA1/1.2,solA1(:,1),'Color',blue,'DisplayName','Basal','LineStyle','--')
plot(tA1/1.2,solA1(:,3),'Color','k','DisplayName','BLI','LineWidth',BLIthickness)
text(70,1e3,"$\bar{r}_B=\bar{r}_L=0.6$","FontSize",fontsize)
hold off
ylim([10^(-3),10^4])
set(gca,'YScale','log','FontSize',fontsize)
xlabel('Time')
ylabel('Cell number')
legend('Location','eastoutside')
exportgraphics(fig1,"C:\Users\bcber\OneDrive\Documents\1-UCI\Q3\Lowengrub Rotation\fig-homeostasis-slides.pdf")

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
set(fig2,'Position',[1000,818,680,420])
hold on
l1 = plot(tA1/1.2,solA1(:,2),'DisplayName','Luminal','Color',orange);
plot(tA2/1.2 + 50,solA2(:,2),'DisplayName','Luminal','Color',orange)
plot(tA3/1.2 + 150,solA3(:,2),'DisplayName','Luminal','Color',orange)
l2 = plot(tA1/1.2,solA1(:,1),'DisplayName','Basal','Color',blue,'LineStyle','--');
plot(tA2/1.2 + 50,solA2(:,1),'DisplayName','Basal','Color',blue,'LineStyle','--')
plot(tA3/1.2 + 150,solA3(:,1),'DisplayName','Basal','Color',blue,'LineStyle','--')
l3 = plot(tA1/1.2,solA1(:,3),'DisplayName','BLI','Color','k','LineWidth',BLIthickness);
plot(tA2/1.2 + 50,solA2(:,3),'DisplayName','BLI','Color','k','LineWidth',BLIthickness)
plot(tA3/1.2 + 150,solA3(:,3),'DisplayName','BLI','Color','k','LineWidth',BLIthickness)
% gridlines to highlight parameter changes
text(10,1e3,["$\bar{r}_B=\bar{r}_L$", "$=0.6$"],"FontSize",fontsize)
plot([50,50],[10^(-3),10^4],'LineStyle','--','Color',[0.5,0.5,0.5],'LineWidth',1)
text(70,1e3,"$\bar{r}_B=\bar{r}_L=1.2$","FontSize",fontsize)
plot([150,150],[10^(-3),10^4],'LineStyle','--','Color',[0.5,0.5,0.5],'LineWidth',1)
text(160,1e3,["$\bar{r}_B=\bar{r}_L$", "$=0.6$"],"FontSize",fontsize)
hold off
ylim([10^(-3),10^4])
xlabel('Time')
ylabel('Cell number')
legend([l1,l2,l3],'Location','eastoutside')
set(gca,'YScale','log','FontSize',fontsize)
exportgraphics(fig2,"C:\Users\bcber\OneDrive\Documents\1-UCI\Q3\Lowengrub Rotation\fig-demandprod-slides.pdf")

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
set(fig3,'Position',[1000,818,680,420])
hold on
plot(tB/1.2,solB(:,2),'DisplayName','Luminal','Color',orange)
plot(tB/1.2,solB(:,1),'DisplayName','Basal','Color',blue,'LineStyle','--')
plot(tB/1.2,solB(:,3),'DisplayName','BLI','Color','k','LineWidth',BLIthickness)
text(70,1e3,"$\bar{r}_B=\bar{r}_L=1.2$","FontSize",fontsize)
hold off
ylim([10^(-3),10^4])
xlabel('Time')
ylabel('Cell number')
legend('Location','eastoutside')
set(gca,'YScale','log','FontSize',fontsize)

exportgraphics(fig3,"C:\Users\bcber\OneDrive\Documents\1-UCI\Q3\Lowengrub Rotation\fig-premalignant-slides.pdf")

%% Malignant plot

initC = [20;20;10];
tspan = [0,1.2*2e2];
[tC, solC] = ode45(sysB, tspan, initC, opts);

fig4 = figure(4);
set(fig4,'Position',[1000,818,680,420])
hold on
plot(tC/1.2,solC(:,2),'DisplayName','Luminal','Color',orange)
plot(tC/1.2,solC(:,1),'DisplayName','Basal','Color',blue,'LineStyle','--')
plot(tC/1.2,solC(:,3),'DisplayName','BLI','Color','k','LineWidth',BLIthickness)
text(70,1e3,"$\bar{r}_B=\bar{r}_L=1.2$","FontSize",fontsize)
hold off
ylim([10^(-3),10^4])
xlabel('Time')
ylabel('Cell number')
legend('Location','eastoutside')
set(gca,'YScale','log','FontSize',fontsize)
exportgraphics(fig4,"C:\Users\bcber\OneDrive\Documents\1-UCI\Q3\Lowengrub Rotation\fig-malignant-slides.pdf")

%% Partial phase portrait: dBLI/dt

sysA = @(t,x) bilayerDuct_05(t,x,params,'nochange');
tspan = [0,1.2*2e2];
initA = [20;20;8.911459];
% initA = [66.252600602161010; 68.252787226175542; 8];
labelA = ['$[B(0);L(0);BLI(0)] = ', mat2str(initA), '$'];
[tA, solA] = ode45(sysA, tspan, initA, opts);
tAscaled = tA / 1.2;
qA = usefulQuantities_05(solA, params);

t1 = find(tAscaled>1,1) - 1; % index where tAscaled approx. 1
t50 = find(tAscaled>50,1) - 1; % index where tAscaled approx. 50
t90 = find(tAscaled>90,1) - 1; % index where tAscaled approx. 90

BLI = linspace(0,20,100);
A = @(bli,SA) (params.Abar + params.gammaAP * bli) ./ ((1 + params.gammaA * SA) ...
    .* (1 + (params.gammaFST * bli).^2));
rI = @(bli,SA) params.rLbar ./ (1 + params.gammaPile * A(bli,SA));
dBLIdt = @(bli,B,SA) params.pPile * B + (rI(bli,SA) - params.rPDbar) .* bli;

fig5 = figure(5);
set(fig5,'Position',[100,100,850,600])
fig5 = tiledlayout(2,3,"TileSpacing","compact");
nexttile(1)
hold on;
plot(tAscaled, solA(:,2), 'DisplayName', 'Luminal', 'LineStyle', '-', ...
    'Color', orange)
plot(tAscaled, solA(:,1), 'DisplayName', 'Basal', 'LineStyle', '--', ...
    'Color', blue)
plot(tAscaled, solA(:,3), 'DisplayName', 'BLI', 'LineStyle', '-', ...
    'Color', 'k')
plot(tAscaled(t1)*ones(1,2),[0,100],'-','LineWidth',1,'Color',[0.5,0.5,0.5])
plot(tAscaled(t50)*ones(1,2),[0,100],'--','LineWidth',1,'Color',[0.5,0.5,0.5])
plot(tAscaled(t90)*ones(1,2),[0,100],'-.','LineWidth',1,'Color',[0.5,0.5,0.5])
hold off;
xlabel('Time','FontSize',10);     ylabel('Cell number','FontSize',10);
nexttile(2,[2,2])
hold on;
l1 = plot(BLI,dBLIdt(BLI,solA(t1,1),qA.SA(t1)),'-k','DisplayName','$t\approx 1$');
l2 = plot(BLI,dBLIdt(BLI,solA(t50,1),qA.SA(t50)),'--k','DisplayName','$t\approx 50$');
l3 = plot(BLI,dBLIdt(BLI,solA(t90,1),qA.SA(t90)),'-.k','DisplayName','$t\approx 90$');
plot(solA(t1,3),dBLIdt(solA(t1,3),solA(t1,1),qA.SA(t1)),'or')
plot(solA(t50,3),dBLIdt(solA(t50,3),solA(t50,1),qA.SA(t50)),'or')
plot(solA(t90,3),dBLIdt(solA(t90,3),solA(t90,1),qA.SA(t90)),'or')
plot([min(BLI),max(BLI)],[0,0],'--','LineWidth',1,'Color',[0.5,0.5,0.5])
hold off;
xlabel('$I$','FontSize',fontsize)
ylabel('$dI/dt$','FontSize',fontsize)
nexttile(4)
hold on;
plot(tAscaled, qA.SA, 'DisplayName', '$S_A$', 'Color', C(3,:))
plot(tAscaled(t1)*ones(1,2),[0,0.5],'-','LineWidth',1,'Color',[0.5,0.5,0.5])
plot(tAscaled(t50)*ones(1,2),[0,0.5],'--','LineWidth',1,'Color',[0.5,0.5,0.5])
plot(tAscaled(t90)*ones(1,2),[0,0.5],'-.','LineWidth',1,'Color',[0.5,0.5,0.5])
hold off;
xlabel('Time','FontSize',10);     ylabel('$S_A(t)$ (gap formation)','FontSize',10);
lg = legend([l1,l2,l3],"FontSize",fontsize);
lg.Layout.Tile = 'east';
exportgraphics(fig5,"C:\Users\bcber\OneDrive\Documents\1-UCI\Q3\Lowengrub Rotation\fig-premaligDetail-slides.pdf")