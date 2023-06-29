% model02_phaseplane.m
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
tCutoff = Inf;
afterCutoff = [];
params = struct('pbar',pbar,'Abar',Abar,'rBbar',rBbar,'rLbar',rLbar,...
    'rLDbar',rLDbar,'rBDbar',rBDbar,'gammaA',gammaA,'gammap',gammap,...
    'gammaB',gammaB,'gammaL',gammaL,'tauSA',tauSA,'aE',aE,'k',k,...
    'tCutoff',tCutoff,'afterCutoff',afterCutoff);
tspan = [0,1e4];
init = [1;1];

%% Phase planes, varying gammaL

nB = 10;
nL = 10;
initB = linspace(0.1, 10, nB);
initL = linspace(0.1, 10, nL);
initVec = cat(3,meshgrid(initB,initL),meshgrid(initB,initL).');

% gammaL = 0.25
paramsA1 = params;
paramsA1.gammaL = 0.25;
solCurvesA1 = cell(nB,nL);
opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
parfor idx = 1:(nB*nL)
    [row,col] = ind2sub([nB,nL],idx);
    init = reshape(initVec(row,col,:),[2,1]);
    sys = @(t,x) bilayerDuct_02(t,x,init,paramsA1,'nochange');
    [~,sol] = ode45(sys, tspan, init, opts);
    solCurvesA1{idx} = sol;
end

% gammaL = 1
paramsA2 = params;
paramsA2.gammaL = 1;
solCurvesA2 = cell(nB,nL);
opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
parfor idx = 1:(nB*nL)
    [row,col] = ind2sub([nB,nL],idx);
    init = reshape(initVec(row,col,:),[2,1]);
    sys = @(t,x) bilayerDuct_02(t,x,init,paramsA2,'nochange');
    [~,sol] = ode45(sys, tspan, init, opts);
    solCurvesA2{idx} = sol;
end

% gammaL = 2.5
paramsA3 = params;
paramsA3.gammaL = 2.5;
solCurvesA3 = cell(nB,nL);
opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
parfor idx = 1:(nB*nL)
    [row,col] = ind2sub([nB,nL],idx);
    init = reshape(initVec(row,col,:),[2,1]);
    sys = @(t,x) bilayerDuct_02(t,x,init,paramsA3,'nochange');
    [~,sol] = ode45(sys, tspan, init, opts);
    solCurvesA3{idx} = sol;
end

% zoomed out version
fig1 = figure(1);
set(fig1,'Position',[100,100,700,1000])
fig1 = tiledlayout(3,2,"TileSpacing","compact");
nexttile(1)
hold on;
for idx = 1:numel(solCurvesA1)
    sol = solCurvesA1{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsA1.gammaL),'$: Trajectories'])
nexttile(2)
hold on;
for idx = 1:numel(solCurvesA1)
    sol = solCurvesA1{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')

nexttile(3)
hold on;
for idx = 1:numel(solCurvesA2)
    sol = solCurvesA2{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsA2.gammaL),'$: Trajectories'])
nexttile(4)
hold on;
for idx = 1:numel(solCurvesA2)
    sol = solCurvesA2{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')

nexttile(5)
hold on;
for idx = 1:numel(solCurvesA3)
    sol = solCurvesA3{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsA3.gammaL),'$: Trajectories'])
nexttile(6)
hold on;
for idx = 1:numel(solCurvesA3)
    sol = solCurvesA3{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')
saveas(fig1,'model02_phase_gammaL_zoomout.svg')
close(gcf)
% zoomed in version
fig2 = figure(2);
set(fig2,'Position',[100,100,700,1000])
fig2 = tiledlayout(3,2,"TileSpacing","compact");
nexttile(1)
hold on;
for idx = 1:numel(solCurvesA1)
    sol = solCurvesA1{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlim([0,10])
ylim([0,10])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsA1.gammaL),'$: Trajectories'])
nexttile(2)
hold on;
for idx = 1:numel(solCurvesA1)
    sol = solCurvesA1{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlim([0,10])
ylim([0,10])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')

nexttile(3)
hold on;
for idx = 1:numel(solCurvesA2)
    sol = solCurvesA2{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlim([0,10])
ylim([0,10])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsA2.gammaL),'$: Trajectories'])
nexttile(4)
hold on;
for idx = 1:numel(solCurvesA2)
    sol = solCurvesA2{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlim([0,10])
ylim([0,10])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')

nexttile(5)
hold on;
for idx = 1:numel(solCurvesA3)
    sol = solCurvesA3{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlim([0,10])
ylim([0,10])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsA3.gammaL),'$: Trajectories'])
nexttile(6)
hold on;
for idx = 1:numel(solCurvesA3)
    sol = solCurvesA3{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlim([0,10])
ylim([0,10])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')
saveas(fig2,'model02_phase_gammaL_zoomin.svg')
close(gcf)

%% Phase planes, varying rLbar

nB = 10;
nL = 10;
initB = linspace(0.1, 10, nB);
initL = linspace(0.1, 10, nL);
initVec = cat(3,meshgrid(initB,initL),meshgrid(initB,initL).');

% rLbar = 1.5
paramsA1 = params;
paramsA1.rLbar = 1.5;
solCurvesA1 = cell(nB,nL);
opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
parfor idx = 1:(nB*nL)
    [row,col] = ind2sub([nB,nL],idx);
    init = reshape(initVec(row,col,:),[2,1]);
    sys = @(t,x) bilayerDuct_02(t,x,init,paramsA1,'nochange');
    [~,sol] = ode45(sys, tspan, init, opts);
    solCurvesA1{idx} = sol;
end

% rLbar = 1.2
paramsA2 = params;
paramsA2.rLbar = 1.2;
solCurvesA2 = cell(nB,nL);
opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
parfor idx = 1:(nB*nL)
    [row,col] = ind2sub([nB,nL],idx);
    init = reshape(initVec(row,col,:),[2,1]);
    sys = @(t,x) bilayerDuct_02(t,x,init,paramsA2,'nochange');
    [~,sol] = ode45(sys, tspan, init, opts);
    solCurvesA2{idx} = sol;
end

% rLbar = 0.9
paramsA3 = params;
paramsA3.rLbar = 0.9;
solCurvesA3 = cell(nB,nL);
opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
parfor idx = 1:(nB*nL)
    [row,col] = ind2sub([nB,nL],idx);
    init = reshape(initVec(row,col,:),[2,1]);
    sys = @(t,x) bilayerDuct_02(t,x,init,paramsA3,'nochange');
    [~,sol] = ode45(sys, tspan, init, opts);
    solCurvesA3{idx} = sol;
end

% zoomed out version
fig3 = figure(3);
set(fig3,'Position',[100,100,700,1000])
fig3 = tiledlayout(3,2,"TileSpacing","compact");
nexttile(1)
hold on;
for idx = 1:numel(solCurvesA1)
    sol = solCurvesA1{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\bar{r}_L = ',num2str(paramsA1.rLbar),'$: Trajectories'])
nexttile(2)
hold on;
for idx = 1:numel(solCurvesA1)
    sol = solCurvesA1{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')

nexttile(3)
hold on;
for idx = 1:numel(solCurvesA2)
    sol = solCurvesA2{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\bar{r}_L = ',num2str(paramsA2.rLbar),'$: Trajectories'])
nexttile(4)
hold on;
for idx = 1:numel(solCurvesA2)
    sol = solCurvesA2{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')

nexttile(5)
hold on;
for idx = 1:numel(solCurvesA3)
    sol = solCurvesA3{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\bar{r}_L = ',num2str(paramsA3.rLbar),'$: Trajectories'])
nexttile(6)
hold on;
for idx = 1:numel(solCurvesA3)
    sol = solCurvesA3{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')
saveas(fig3,'model02_phase_rLbar_zoomout.svg')
close(gcf)
% zoomed in version
fig4 = figure(4);
set(fig4,'Position',[100,100,700,1000])
fig4 = tiledlayout(3,2,"TileSpacing","compact");
nexttile(1)
hold on;
for idx = 1:numel(solCurvesA1)
    sol = solCurvesA1{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlim([0,10])
ylim([0,10])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\bar{r}_L = ',num2str(paramsA1.rLbar),'$: Trajectories'])
nexttile(2)
hold on;
for idx = 1:numel(solCurvesA1)
    sol = solCurvesA1{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlim([0,10])
ylim([0,10])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')

nexttile(3)
hold on;
for idx = 1:numel(solCurvesA2)
    sol = solCurvesA2{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlim([0,10])
ylim([0,10])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\bar{r}_L = ',num2str(paramsA2.rLbar),'$: Trajectories'])
nexttile(4)
hold on;
for idx = 1:numel(solCurvesA2)
    sol = solCurvesA2{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlim([0,10])
ylim([0,10])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')

nexttile(5)
hold on;
for idx = 1:numel(solCurvesA3)
    sol = solCurvesA3{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlim([0,10])
ylim([0,10])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\bar{r}_L = ',num2str(paramsA3.rLbar),'$: Trajectories'])
nexttile(6)
hold on;
for idx = 1:numel(solCurvesA3)
    sol = solCurvesA3{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlim([0,10])
ylim([0,10])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')
saveas(fig4,'model02_phase_rLbar_zoomin.svg')
close(gcf)

%% Phase planes, varying gammaL -- nonzero basal death

nB = 10;
nL = 10;
initB = linspace(0.1, 20, nB);
initL = linspace(0.1, 20, nL);
initVec = cat(3,meshgrid(initB,initL),meshgrid(initB,initL).');

% gammaL = 0.25
paramsA1 = params;
paramsA1.gammaL = 0.25;
paramsA1.rBDbar = 0.05;
solCurvesA1 = cell(nB,nL);
opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
parfor idx = 1:(nB*nL)
    [row,col] = ind2sub([nB,nL],idx);
    init = reshape(initVec(row,col,:),[2,1]);
    sys = @(t,x) bilayerDuct_02(t,x,init,paramsA1,'nochange');
    [~,sol] = ode45(sys, tspan, init, opts);
    solCurvesA1{idx} = sol;
end

% gammaL = 1
paramsA2 = params;
paramsA2.gammaL = 1;
paramsA2.rBDbar = 0.05;
solCurvesA2 = cell(nB,nL);
opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
parfor idx = 1:(nB*nL)
    [row,col] = ind2sub([nB,nL],idx);
    init = reshape(initVec(row,col,:),[2,1]);
    sys = @(t,x) bilayerDuct_02(t,x,init,paramsA2,'nochange');
    [~,sol] = ode45(sys, tspan, init, opts);
    solCurvesA2{idx} = sol;
end

% gammaL = 2.5
paramsA3 = params;
paramsA3.gammaL = 2.5;
paramsA3.rBDbar = 0.05;
solCurvesA3 = cell(nB,nL);
opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
parfor idx = 1:(nB*nL)
    [row,col] = ind2sub([nB,nL],idx);
    init = reshape(initVec(row,col,:),[2,1]);
    sys = @(t,x) bilayerDuct_02(t,x,init,paramsA3,'nochange');
    [~,sol] = ode45(sys, tspan, init, opts);
    solCurvesA3{idx} = sol;
end

% zoomed out version
fig5 = figure(5);
set(fig5,'Position',[100,100,700,1000])
fig5 = tiledlayout(3,2,"TileSpacing","compact");
nexttile(1)
hold on;
for idx = 1:numel(solCurvesA1)
    sol = solCurvesA1{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsA1.gammaL),'$: Trajectories'])
nexttile(2)
hold on;
for idx = 1:numel(solCurvesA1)
    sol = solCurvesA1{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')

nexttile(3)
hold on;
for idx = 1:numel(solCurvesA2)
    sol = solCurvesA2{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsA2.gammaL),'$: Trajectories'])
nexttile(4)
hold on;
for idx = 1:numel(solCurvesA2)
    sol = solCurvesA2{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')

nexttile(5)
hold on;
for idx = 1:numel(solCurvesA3)
    sol = solCurvesA3{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsA3.gammaL),'$: Trajectories'])
nexttile(6)
hold on;
for idx = 1:numel(solCurvesA3)
    sol = solCurvesA3{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')
saveas(fig5,'model02_phase_gammaL+basaldeath_zoomout.svg')
close(gcf)
% zoomed in version
fig6 = figure(6);
set(fig6,'Position',[100,100,700,1000])
fig6 = tiledlayout(3,2,"TileSpacing","compact");
nexttile(1)
hold on;
for idx = 1:numel(solCurvesA1)
    sol = solCurvesA1{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlim([0,max(initB)])
ylim([0,max(initL)])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsA1.gammaL),'$: Trajectories'])
nexttile(2)
hold on;
for idx = 1:numel(solCurvesA1)
    sol = solCurvesA1{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlim([0,max(initB)])
ylim([0,max(initL)])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')

nexttile(3)
hold on;
for idx = 1:numel(solCurvesA2)
    sol = solCurvesA2{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlim([0,max(initB)])
ylim([0,max(initL)])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsA2.gammaL),'$: Trajectories'])
nexttile(4)
hold on;
for idx = 1:numel(solCurvesA2)
    sol = solCurvesA2{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlim([0,max(initB)])
ylim([0,max(initL)])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')

nexttile(5)
hold on;
for idx = 1:numel(solCurvesA3)
    sol = solCurvesA3{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlim([0,max(initB)])
ylim([0,max(initL)])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsA3.gammaL),'$: Trajectories'])
nexttile(6)
hold on;
for idx = 1:numel(solCurvesA3)
    sol = solCurvesA3{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlim([0,max(initB)])
ylim([0,max(initL)])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')
saveas(fig6,'model02_phase_gammaL+basaldeath_zoomin.svg')
close(gcf)

%% Phase planes, varying rLbar -- nonzero basal death

nB = 10;
nL = 10;
initB = linspace(0.1, 10, nB);
initL = linspace(0.1, 10, nL);
initVec = cat(3,meshgrid(initB,initL),meshgrid(initB,initL).');

% rLbar = 1.5
paramsA1 = params;
paramsA1.rLbar = 1.5;
paramsA1.rBDbar = 0.05;
solCurvesA1 = cell(nB,nL);
opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
parfor idx = 1:(nB*nL)
    [row,col] = ind2sub([nB,nL],idx);
    init = reshape(initVec(row,col,:),[2,1]);
    sys = @(t,x) bilayerDuct_02(t,x,init,paramsA1,'nochange');
    [~,sol] = ode45(sys, tspan, init, opts);
    solCurvesA1{idx} = sol;
end

% rLbar = 1.2
paramsA2 = params;
paramsA2.rLbar = 1.2;
paramsA2.rBDbar = 0.05;
solCurvesA2 = cell(nB,nL);
opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
parfor idx = 1:(nB*nL)
    [row,col] = ind2sub([nB,nL],idx);
    init = reshape(initVec(row,col,:),[2,1]);
    sys = @(t,x) bilayerDuct_02(t,x,init,paramsA2,'nochange');
    [~,sol] = ode45(sys, tspan, init, opts);
    solCurvesA2{idx} = sol;
end

% rLbar = 0.9
paramsA3 = params;
paramsA3.rLbar = 0.9;
paramsA3.rBDbar = 0.05;
solCurvesA3 = cell(nB,nL);
opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
parfor idx = 1:(nB*nL)
    [row,col] = ind2sub([nB,nL],idx);
    init = reshape(initVec(row,col,:),[2,1]);
    sys = @(t,x) bilayerDuct_02(t,x,init,paramsA3,'nochange');
    [~,sol] = ode45(sys, tspan, init, opts);
    solCurvesA3{idx} = sol;
end

% zoomed out version
fig7 = figure(7);
set(fig7,'Position',[100,100,700,1000])
fig7 = tiledlayout(3,2,"TileSpacing","compact");
nexttile(1)
hold on;
for idx = 1:numel(solCurvesA1)
    sol = solCurvesA1{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\bar{r}_L = ',num2str(paramsA1.rLbar),'$: Trajectories'])
nexttile(2)
hold on;
for idx = 1:numel(solCurvesA1)
    sol = solCurvesA1{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')

nexttile(3)
hold on;
for idx = 1:numel(solCurvesA2)
    sol = solCurvesA2{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\bar{r}_L = ',num2str(paramsA2.rLbar),'$: Trajectories'])
nexttile(4)
hold on;
for idx = 1:numel(solCurvesA2)
    sol = solCurvesA2{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')

nexttile(5)
hold on;
for idx = 1:numel(solCurvesA3)
    sol = solCurvesA3{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\bar{r}_L = ',num2str(paramsA3.rLbar),'$: Trajectories'])
nexttile(6)
hold on;
for idx = 1:numel(solCurvesA3)
    sol = solCurvesA3{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')
saveas(fig7,'model02_phase_rLbar+basaldeath_zoomout.svg')
close(gcf)
% zoomed in version
fig8 = figure(8);
set(fig8,'Position',[100,100,700,1000])
fig8 = tiledlayout(3,2,"TileSpacing","compact");
nexttile(1)
hold on;
for idx = 1:numel(solCurvesA1)
    sol = solCurvesA1{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlim([0,10])
ylim([0,10])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\bar{r}_L = ',num2str(paramsA1.rLbar),'$: Trajectories'])
nexttile(2)
hold on;
for idx = 1:numel(solCurvesA1)
    sol = solCurvesA1{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlim([0,10])
ylim([0,10])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')

nexttile(3)
hold on;
for idx = 1:numel(solCurvesA2)
    sol = solCurvesA2{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlim([0,10])
ylim([0,10])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\bar{r}_L = ',num2str(paramsA2.rLbar),'$: Trajectories'])
nexttile(4)
hold on;
for idx = 1:numel(solCurvesA2)
    sol = solCurvesA2{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlim([0,10])
ylim([0,10])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')

nexttile(5)
hold on;
for idx = 1:numel(solCurvesA3)
    sol = solCurvesA3{idx};
    plot(sol(:,1),sol(:,2),'-','Color',blue)
    plot(sol(end,1),sol(end,2),'.r')
end
hold off;
xlim([0,10])
ylim([0,10])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\bar{r}_L = ',num2str(paramsA3.rLbar),'$: Trajectories'])
nexttile(6)
hold on;
for idx = 1:numel(solCurvesA3)
    sol = solCurvesA3{idx};
    plot(sol(end-50:end,1),sol(end-50:end,2),'.r')
end
hold off;
xlim([0,10])
ylim([0,10])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Steady states only')
saveas(fig8,'model02_phase_rLbar+basaldeath_zoomin.svg')
close(gcf)