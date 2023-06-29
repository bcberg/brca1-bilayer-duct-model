% model04.m
% Brady Berg
clear; close all;
format long; format compact;
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')

blue = "#0072BD";
orange = "#D95319";
yellow = "#EDB120";
green = "#77AC30";
purple = "#7E2F8E";

%% Base parameters

pbar = 1;
Abar = 0.1;
rBbar = 1;
rLbar = 3;
rLDbar = 1;
rBDbar = 0.05;
gammaA = 10;
nA = 1;
gammap = 10;
gammaB = 40;
gammaL = 1;
tauSA = 0.1;
aE = 0.075;
k = 4;
tCutoff = Inf;
afterCutoff = [];
pPile = 0.05;
gammaPile = 20;
gammaRank = 1/2;
nRank = 2;
gammaAtilde = 1;
alpha = 10;
rPDbar = 1;

params = struct('pbar',pbar,'Abar',Abar,'rBbar',rBbar,'rLbar',rLbar,...
    'rLDbar',rLDbar,'rBDbar',rBDbar,'gammaA',gammaA,'nA',nA, ...
    'gammap',gammap,'gammaB',gammaB,'gammaL',gammaL,'tauSA',tauSA, ...
    'aE',aE,'k',k,'tCutoff',tCutoff,'afterCutoff',afterCutoff, ...
    'pPile',pPile,'gammaPile',gammaPile,'gammaRank',gammaRank, ...
    'nRank',nRank,'gammaAtilde',gammaAtilde,'alpha',alpha,'rPDbar',rPDbar);

opts = odeset('RelTol',1e-5,'AbsTol',1e-6,'Events',@detectOvergrowth);

%% Looking at time courses

paramsA = params;
tspan = [0,40];
sysA = @(t,x) bilayerDuct_04(t,x,paramsA,'nochange');
[tA1, solA1] = ode45(sysA, tspan, [4;3;0], opts);
[tA2, solA2] = ode45(sysA, tspan, [5;4;0], opts);

fig1 = figure(1);
set(fig1,'Position',[100,100,900,450])
tiledlayout(1,2,'TileSpacing','compact')
nexttile(1)
hold on;
plot(tA1, solA1(:,1), 'DisplayName', 'Basal')
plot(tA1, solA1(:,2), 'DisplayName', 'Luminal')
plot(tA1, solA1(:,3), 'DisplayName', '$L_{\rm pile}$')
plot(tA2, solA2(:,1), '--', 'Color', [0.5 0.5 0.5])
plot(tA2, solA2(:,2), '--', 'Color', [0.5 0.5 0.5])
plot(tA2, solA2(:,3), '--', 'Color', [0.5 0.5 0.5])
hold off;
xlabel('Time')
ylabel('Cell number')
xlim([0,40])
ylim([0,12])
title('$[B,L,L_{\rm pile}](0) = [4,3,0]$')
nexttile(2)
hold on;
plot(tA2, solA2(:,1), 'DisplayName', 'Basal')
plot(tA2, solA2(:,2), 'DisplayName', 'Luminal')
plot(tA2, solA2(:,3), 'DisplayName', '$L_{\rm pile}$')
hold off;
lg = legend;
lg.Layout.Tile = 'south';
xlabel('Time')
ylabel('Cell number')
title('$[B,L,L_{\rm pile}](0) = [5,4,0]$')
exportgraphics(fig1, ...
    "C:\Users\bcber\OneDrive\Documents\1-UCI\Q3\Lowengrub Rotation\model04_time-1.tiff")

%% Phase planes

tspanPhase = [0,5e2];
n = 15;
initB = linspace(0.1,20,n);
initL = initB;
initVec = cat(3,meshgrid(initB,initL).',meshgrid(initB,initL),zeros(n,n));
paramsB = params;
sysB = @(t,x) bilayerDuct_04(t,x,paramsB,'nochange');

solCurvesB = cell(n,n);
overgrown = zeros(n,n);
parfor idx = 1:n^2
    [row,col] = ind2sub([n,n],idx);
    init = reshape(initVec(row,col,:),[3,1]);
    [~, sol] = ode45(sysB, tspanPhase, init, opts);
    if sol(end,3) > 1e3
        overgrown(idx) = 1;
    end
    solCurvesB{idx} = sol;
end
overgrown = logical(overgrown);

fig2 = figure(2);
set(fig2,'Position',[110,110,900,450])
tiledlayout(1,2,'TileSpacing','compact');
nexttile(1)
hold on;
for idx = 1:numel(solCurvesB)
    sol = solCurvesB{idx};
    plot(sol(:,1), sol(:,2) + sol(:,3), 'Color', blue)
    plot(sol(end,1), sol(end,2) + sol(end,3), '.r')
end
hold off;
xlim([0,20])
ylim([0,20])
xlabel('Basal cell number')
ylabel('Total luminal cells')
title('Trajectories')
nexttile(2)
initX = initVec(:,:,1);
initY = initVec(:,:,2);
hold on;
plot(initX(overgrown),initY(overgrown),'*r')
plot(initX(~overgrown),initY(~overgrown),'.k')
hold off;
xlim([0,20])
ylim([0,20])
xlabel('Initial basal cell number')
ylabel('Initial luminal cell number')
title('Marking overgrowth')
exportgraphics(fig2, ...
    "C:\Users\bcber\OneDrive\Documents\1-UCI\Q3\Lowengrub Rotation\model04_phase-1.tiff")

%% Functions

function [value, isterminal, direction] = detectOvergrowth(~,y)
value = y(3) - 1e6;     % when value = 0 (Lpile reaches 1e6), ...
isterminal = 1;         % terminate integration
direction = 0;
end