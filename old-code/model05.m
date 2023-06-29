% model05.m
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
% saved parameter set 1
% params.pbar = 1;
% params.Abar = 0.1;
% params.rBbar = 1;
% params.rLbar = 1;
% params.rBDbar = 0.04;
% params.rLDbar = 0.05;
% params.rPDbar = 1;
% params.gammap = 10;
% params.gammaA = 20;
% params.gammaB = 2;
% params.gammaL = 2;
% params.pPile = 1e-3;    % 1e-3
% params.gammaPile = 20;
% params.gammaAP = 5;
% params.gammaFST = 290;
% params.nA = 1;
% params.nFST = 1;
% params.tauSA = 0.1;
% params.aE = 0.075;
% params.k = 4;
% params.alpha = 0.5;
% params.gammaPD = 1e-2;

% params.pbar = 1;
% params.Abar = 0.1;
% params.rBbar = 1.5;
% params.rLbar = 1.5;
% params.rBDbar = 0.04;
% params.rLDbar = 0.05;
% params.rPDbar = 1;
% params.gammap = 10;
% params.gammaA = 20;
% params.gammaB = 2;
% params.gammaL = 2;
% params.pPile = 1e-3;    % 1e-3
% params.gammaPile = 20;
% params.gammaAP = 5;
% params.gammaFST = 0.1;
% params.nA = 1;
% params.nFST = 1;
% params.tauSA = 0.1;
% params.aE = 0.075;
% params.k = 4;
% params.alpha = 0.5;
% params.gammaPD = 1e-2;

% saved parameters set 2
% params.pbar = 1;
% params.Abar = 0.1;
% params.rBbar = 1;
% params.rLbar = 2;
% params.rBDbar = 0.05;
% params.rLDbar = 0.5;
% params.rPDbar = 1;
% params.gammap = 10;
% params.gammaA = 20;
% params.gammaB = 40;
% params.gammaL = 2;
% params.pPile = 1e-4;    % 1e-3
% params.gammaPile = 20;
% params.gammaAP = 5;
% params.gammaFST = 3;
% params.nA = 1;
% params.nFST = 1;
% params.tauSA = 0.1;
% params.aE = 0.075;
% params.k = 4;
% params.alpha = 0.5;
% params.gammaPD = 1e-2;

params.pbar = 1;
params.Abar = 0.1;
params.rBbar = 1;
params.rLbar = 1.3; % 1.1, 1.3, 1.5
params.rBDbar = 0.1;    %0.1
params.rLDbar = 0.25;
params.rPDbar = 1;
params.gammap = 10;
params.gammaA = 1;
params.gammaB = 20;
params.gammaL = 10;  %10
params.pPile = 0.9e-3;    % 1e-3
params.gammaPile = 10;  % 20
params.gammaAP = 1;
params.gammaFST = 1;
% params.gammaFSTstar = 1;
params.nA = 1;
params.nFST = 2;
params.tauSA = 0.1;
params.aE = 0.05;     %0.075;
params.k = 5;
params.alpha = 0;
% params.gammaPD = 1e-2;

% optional parameter step-discontinuity
params.tCutoff = Inf;
params.afterCutoff = [];

opts = odeset('RelTol',1e-5,'AbsTol',1e-6,'Events',@detectOvergrowth);

%% Looking at time courses
close all;
paramsA1 = params;
% paramsA1.rBbar = 1.1;
paramsA1.rLbar = 1.2;
% paramsA1.gammaFST = 60;
paramsA2 = params;
% paramsA2.rBbar = 1.3;
paramsA2.rLbar = 1.2;
% paramsA2.gammaFST = 120;
paramsA3 = params;
% paramsA3.rBbar = 1.5;
paramsA3.rLbar = 1.3;
% paramsA3.gammaFST = 240;
tspan = [0,2e2];
sysA1 = @(t,x) bilayerDuct_05(t,x,paramsA1,'nochange');
sysA2 = @(t,x) bilayerDuct_05(t,x,paramsA2,'nochange');
sysA3 = @(t,x) bilayerDuct_05(t,x,paramsA3,'nochange');
initA1 = [20;20;0];
initA2 = [20;20;15.5];
initA3 = [20;20;15.5];
initlabel = ['$[B,L,L_{\rm pile}](0)=[', num2str(initA1(1)), ',', ...
    num2str(initA1(2)), ',', num2str(initA1(3)), ']$, $\alpha = ', ...
    num2str(params.alpha), '$'];
A1label = ['$\bar{r}_L =', num2str(paramsA1.rLbar),'$'];
A2label = ['$\bar{r}_L =', num2str(paramsA2.rLbar),'$'];
A3label = ['$\bar{r}_L =', num2str(paramsA3.rLbar),'$'];
% initA1 = [5;5;0];
% A1label = ['$[B,L,L_{\rm pile}](0) = [', num2str(initA1(1)), ',', ...
%     num2str(initA1(2)), ',', num2str(initA1(3)), ']$',', $\alpha = ', ...
%     num2str(paramsA1.alpha), '$'];
% initA2 = [5;5;0];
% A2label = ['$[B,L,L_{\rm pile}](0) = [', num2str(initA2(1)), ',', ...
%     num2str(initA2(2)), ',', num2str(initA2(3)), ']$',', $\alpha = ', ...
%     num2str(paramsA2.alpha), '$'];
[tA1, solA1] = ode45(sysA1, tspan, initA1, opts);
[tA2, solA2] = ode45(sysA2, tspan, initA2, opts);
[tA3, solA3] = ode45(sysA3, tspan, initA3, opts);

fig1 = figure(1);
set(fig1,'Position',[100,100,900,400])
fig1 = tiledlayout(1,3,'TileSpacing','compact');
nexttile(1)
hold on;
plot(tA1, solA1(:,1), 'DisplayName', 'Basal')
plot(tA1, solA1(:,2), 'DisplayName', 'Luminal')
plot(tA1, solA1(:,3), 'DisplayName', '$L_{\rm pile}$')
% plot(tA2, solA2(:,1), '--', 'Color', [0.5 0.5 0.5])
% plot(tA2, solA2(:,2), '--', 'Color', [0.5 0.5 0.5])
% plot(tA2, solA2(:,3), '--', 'Color', [0.5 0.5 0.5])
hold off;
% set(gca,'YScale','log')
ylim([0,200])
xlabel('Time')
ylabel('Cell number')
title(A1label)
nexttile(2)
hold on;
plot(tA2, solA2(:,1), 'DisplayName', 'Basal')
plot(tA2, solA2(:,2), 'DisplayName', 'Luminal')
plot(tA2, solA2(:,3), 'DisplayName', '$L_{\rm pile}$')
hold off;
% set(gca,'YScale','log')
ylim([0,200])
xlabel('Time')
ylabel('Cell number')
title(A2label)
nexttile(3)
hold on;
plot(tA3, solA3(:,1), 'DisplayName', 'Basal')
plot(tA3, solA3(:,2), 'DisplayName', 'Luminal')
plot(tA3, solA3(:,3), 'DisplayName', '$L_{\rm pile}$')
hold off;
% set(gca,'YScale','log')
ylim([0,200])
xlabel('Time')
ylabel('Cell number')
title(A3label)
lg = legend;
lg.Layout.Tile = 'south';
% title(fig1,initlabel,'Interpreter', 'latex')
exportgraphics(fig1,...
    "C:\Users\bcber\OneDrive\Documents\1-UCI\Q3\Lowengrub Rotation\model05_time-1.tiff")

q1 = usefulQuantities_05(solA1,paramsA1);
q2 = usefulQuantities_05(solA2,paramsA2);
q3 = usefulQuantities_05(solA3,paramsA3);
% solA1 detail
fig3 = figure(3);
set(fig3,'Position',[120,120,600,600])
fig3 = tiledlayout(2,2,'TileSpacing','compact');
title(fig3,A1label,'Interpreter','latex')
nexttile(1)
hold on;
plot(tA1, solA1(:,1), 'DisplayName', 'Basal')
plot(tA1, solA1(:,2), 'DisplayName', 'Luminal')
plot(tA1, solA1(:,3), 'DisplayName', '$L_{\rm pile}$')
hold off;
% set(gca,'YScale','log')
xlabel('Time');     ylabel('Cell number');
legend('Location','east')
nexttile(2)
hold on;
plot(tA1, q1.SI, 'DisplayName', '$S_I$')
plot(tA1, q1.SE, 'DisplayName', '$S_E$')
hold off;
ylabel('Stress level');
xlabel('Time');     
legend('Location','east')
nexttile(3)
hold on;
plot(tA1, q1.rB, 'DisplayName', '$r_B$')
plot(tA1, q1.rL, 'DisplayName', '$r_L$')
plot(tA1, q1.rPile, 'DisplayName', '$r_{\rm pile}$')
hold off;
xlabel('Time');     ylabel('Proliferation rate');
legend('Location','east')
nexttile(4)
yyaxis left
plot(tA1, q1.p, 'DisplayName', 'p')
ylabel('Basal self-renewal prob.')
yyaxis right
plot(tA1, q1.A, 'DisplayName', 'A')
ylabel('Activin level')
xlabel('Time')
exportgraphics(fig3,...
    "C:\Users\bcber\OneDrive\Documents\1-UCI\Q3\Lowengrub Rotation\model05_time-detail-1.tiff")

% solA2 detail
fig4 = figure(4);
set(fig4,'Position',[130,130,600,600])
fig4 = tiledlayout(2,2,'TileSpacing','compact');
title(fig4,A2label,'Interpreter','latex')
nexttile(1)
hold on;
plot(tA2, solA2(:,1), 'DisplayName', 'Basal')
plot(tA2, solA2(:,2), 'DisplayName', 'Luminal')
plot(tA2, solA2(:,3), 'DisplayName', '$L_{\rm pile}$')
hold off;
% set(gca,'YScale','log')
xlabel('Time');     ylabel('Cell number');
legend('Location','east')
nexttile(2)
hold on;
plot(tA2, q2.SI, 'DisplayName', '$S_I$')
plot(tA2, q2.SE, 'DisplayName', '$S_E$')
hold off;
ylabel('Stress level');
xlabel('Time');     
legend('Location','east')
nexttile(3)
hold on;
plot(tA2, q2.rB, 'DisplayName', '$r_B$')
plot(tA2, q2.rL, 'DisplayName', '$r_L$')
plot(tA2, q2.rPile, 'DisplayName', '$r_{\rm pile}$')
hold off;
xlabel('Time');     ylabel('Proliferation rate');
legend('Location','northeast')
nexttile(4)
yyaxis left
plot(tA2, q2.p, 'DisplayName', 'p')
ylabel('Basal self-renewal prob.')
yyaxis right
plot(tA2, q2.A, 'DisplayName', 'A')
ylabel('Activin level')
xlabel('Time')
exportgraphics(fig4,...
    "C:\Users\bcber\OneDrive\Documents\1-UCI\Q3\Lowengrub Rotation\model05_time-detail-2.tiff")

% solA3 detail
fig5 = figure(5);
set(fig5,'Position',[140,140,600,600])
fig5 = tiledlayout(2,2,'TileSpacing','compact');
title(fig5,A3label,'Interpreter','latex')
nexttile(1)
hold on;
plot(tA3, solA3(:,1), 'DisplayName', 'Basal')
plot(tA3, solA3(:,2), 'DisplayName', 'Luminal')
plot(tA3, solA3(:,3), 'DisplayName', '$L_{\rm pile}$')
hold off;
% set(gca,'YScale','log')
xlabel('Time');     ylabel('Cell number');
legend('Location','northwest')
nexttile(2)
hold on;
plot(tA3, q3.SI, 'DisplayName', '$S_I$')
plot(tA3, q3.SE, 'DisplayName', '$S_E$')
hold off;
ylabel('Stress level');
xlabel('Time');     
legend('Location','east')
nexttile(3)
hold on;
plot(tA3, q3.rB, 'DisplayName', '$r_B$')
plot(tA3, q3.rL, 'DisplayName', '$r_L$')
plot(tA3, q3.rPile, 'DisplayName', '$r_{\rm pile}$')
hold off;
xlabel('Time');     ylabel('Proliferation rate');
legend('Location','northeast')
nexttile(4)
yyaxis left
plot(tA3, q3.p, 'DisplayName', 'p')
ylabel('Basal self-renewal prob.')
yyaxis right
plot(tA3, q3.A, 'DisplayName', 'A')
ylabel('Activin level')
xlabel('Time')
exportgraphics(fig5,...
    "C:\Users\bcber\OneDrive\Documents\1-UCI\Q3\Lowengrub Rotation\model05_time-detail-3.tiff")

%% Phase portrait

% tspanPhase = [0,5e2];
% n = 20;
% initB = linspace(0.1,40,n);
% initL = initB;
% initVec = cat(3,meshgrid(initB,initL).',meshgrid(initB,initL),zeros(n,n));
% paramsB = params;
% sysB = @(t,x) bilayerDuct_05(t,x,paramsB,'nochange');
% 
% solCurvesB = cell(n,n);
% overgrown = zeros(n,n);
% parfor idx = 1:n^2
%     [row,col] = ind2sub([n,n],idx);
%     initA1 = reshape(initVec(row,col,:),[3,1]);
%     [~, sol] = ode45(sysB, tspanPhase, initA1, opts);
%     if sol(end,3) > 1e8
%         overgrown(idx) = 1;
%     end
%     solCurvesB{idx} = sol;
% end
% overgrown = logical(overgrown);
% 
% fig2 = figure(2);
% set(fig2,'Position',[110,110,900,450])
% tiledlayout(1,2,'TileSpacing','compact');
% nexttile(1)
% hold on;
% for idx = 1:numel(solCurvesB)
%     sol = solCurvesB{idx};
%     plot(sol(:,1), sol(:,2) + sol(:,3), 'Color', blue)
%     plot(sol(end,1), sol(end,2) + sol(end,3), '.r')
% end
% hold off;
% xlim([0,40])
% ylim([0,40])
% xlabel('Basal cell number')
% ylabel('Total luminal cells')
% title('Trajectories')
% nexttile(2)
% initX = initVec(:,:,1);
% initY = initVec(:,:,2);
% hold on;
% plot(initX(overgrown),initY(overgrown),'*r')
% plot(initX(~overgrown),initY(~overgrown),'.k')
% hold off;
% xlim([0,40])
% ylim([0,40])
% xlabel('Initial basal cell number')
% ylabel('Initial luminal cell number')
% title('Marking overgrowth')
% % exportgraphics(fig2,...
% %     "C:\Users\bcber\OneDrive\Documents\1-UCI\Q3\Lowengrub Rotation\model05_phase-1.tiff")

%% Functions

function [value, isterminal, direction] = detectOvergrowth(~,y)
value = y(3) - 1e9;     % when value = 0 (Lpile reaches 1e9), ...
isterminal = 1;         % terminate integration
direction = 0;
end