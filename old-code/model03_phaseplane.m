% model03_phaseplane.m
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
rBDbar = 0.05;
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
tspan = [0,1e3];

%% Phase planes, varying gammaL -- nonzero basal death, no piling

n = 15;
initB = linspace(0.1, 20, n);
initL = linspace(0.1, 20, n);
initVec = cat(3,meshgrid(initB,initL),meshgrid(initB,initL).');
opts = odeset('RelTol',1e-5,'AbsTol',1e-6);

% gammaL = 0.25
paramsA1 = params;
paramsA1.gammaL = 0.25;
paramsA1.rBDbar = 0.05;
solCurvesA1 = cell(n);
parfor idx = 1:n^2
    [row,col] = ind2sub([n,n],idx);
    init = [reshape(initVec(row,col,:),[2,1]);0];
    sys = @(t,x) bilayerDuct_03(t,x,paramsA1,'nochange');
    [~,sol] = ode45(sys, tspan, init, opts);
    solCurvesA1{idx} = sol;
end

% gammaL = 1
paramsA2 = params;
paramsA2.gammaL = 1;
paramsA2.rBDbar = 0.05;
solCurvesA2 = cell(n);
parfor idx = 1:n^2
    [row,col] = ind2sub([n,n],idx);
    init = [reshape(initVec(row,col,:),[2,1]);0];
    sys = @(t,x) bilayerDuct_03(t,x,paramsA2,'nochange');
    [~,sol] = ode45(sys, tspan, init, opts);
    solCurvesA2{idx} = sol;
end

% gammaL = 2.5
paramsA3 = params;
paramsA3.gammaL = 2.5;
paramsA3.rBDbar = 0.05;
solCurvesA3 = cell(n);
parfor idx = 1:n^2
    [row,col] = ind2sub([n,n],idx);
    init = [reshape(initVec(row,col,:),[2,1]);0];
    sys = @(t,x) bilayerDuct_03(t,x,paramsA3,'nochange');
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
    plot(sol(:,1),sol(:,2)+sol(:,3),'-','Color',blue)
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
axis tight
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsA1.gammaL),'$: Trajectories'])
nexttile(2)
hold on;
for idx = 1:numel(solCurvesA1)
    sol = solCurvesA1{idx};
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
axis tight
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Endpoints only')

nexttile(3)
hold on;
for idx = 1:numel(solCurvesA2)
    sol = solCurvesA2{idx};
    plot(sol(:,1),sol(:,2)+sol(:,3),'-','Color',blue)
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
axis tight
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsA2.gammaL),'$: Trajectories'])
nexttile(4)
hold on;
for idx = 1:numel(solCurvesA2)
    sol = solCurvesA2{idx};
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
axis tight
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Endpoints only')

nexttile(5)
hold on;
for idx = 1:numel(solCurvesA3)
    sol = solCurvesA3{idx};
    plot(sol(:,1),sol(:,2)+sol(:,3),'-','Color',blue)
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
axis tight
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsA3.gammaL),'$: Trajectories'])
nexttile(6)
hold on;
for idx = 1:numel(solCurvesA3)
    sol = solCurvesA3{idx};
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
axis tight
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Endpoints only')
exportgraphics(fig1,'model03_phase_gammaL+basaldeath_zoomout.pdf',...
    'ContentType','vector')
% close(gcf)
% zoomed in version
fig2 = figure(2);
set(fig2,'Position',[100,100,700,1000])
fig2 = tiledlayout(3,2,"TileSpacing","compact");
nexttile(1)
hold on;
for idx = 1:numel(solCurvesA1)
    sol = solCurvesA1{idx};
    plot(sol(:,1),sol(:,2)+sol(:,3),'-','Color',blue)
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
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
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
xlim([0,max(initB)])
ylim([0,max(initL)])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Endpoints only')

nexttile(3)
hold on;
for idx = 1:numel(solCurvesA2)
    sol = solCurvesA2{idx};
    plot(sol(:,1),sol(:,2)+sol(:,3),'-','Color',blue)
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
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
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
xlim([0,max(initB)])
ylim([0,max(initL)])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Endpoints only')

nexttile(5)
hold on;
for idx = 1:numel(solCurvesA3)
    sol = solCurvesA3{idx};
    plot(sol(:,1),sol(:,2)+sol(:,3),'-','Color',blue)
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
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
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
xlim([0,max(initB)])
ylim([0,max(initL)])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Endpoints only')
exportgraphics(fig2,'model03_phase_gammaL+basaldeath_zoomin.pdf',...
    'ContentType','vector')
% close(gcf)

%% Phase planes, varying gammaL -- nonzero basal death, some piling

n = 10;
initB = linspace(0.1, 20, n);
initL = linspace(0.1, 20, n);
initVec = cat(3,meshgrid(initB,initL),meshgrid(initB,initL).');
opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
paramsB = params;
paramsB.pPile = 0.05;

% gammaL = 0.25
paramsB1 = paramsB;
paramsB1.gammaL = 0.25;
solCurvesB1 = cell(n);
parfor idx = 1:n^2
    [row,col] = ind2sub([n,n],idx);
    init = [reshape(initVec(row,col,:),[2,1]);0];
    sys = @(t,x) bilayerDuct_03(t,x,paramsB1,'nochange');
    [~,sol] = ode45(sys, tspan, init, opts);
    solCurvesB1{idx} = sol;
end

% gammaL = 1
paramsB2 = paramsB;
paramsB2.gammaL = 1;
solCurvesB2 = cell(n);
parfor idx = 1:n^2
    [row,col] = ind2sub([n,n],idx);
    init = [reshape(initVec(row,col,:),[2,1]);0];
    sys = @(t,x) bilayerDuct_03(t,x,paramsB2,'nochange');
    [~,sol] = ode45(sys, tspan, init, opts);
    solCurvesB2{idx} = sol;
end

% gammaL = 2.5
paramsB3 = paramsB;
paramsB3.gammaL = 2.5;
solCurvesB3 = cell(n);
parfor idx = 1:n^2
    [row,col] = ind2sub([n,n],idx);
    init = [reshape(initVec(row,col,:),[2,1]);0];
    sys = @(t,x) bilayerDuct_03(t,x,paramsB3,'nochange');
    [~,sol] = ode45(sys, tspan, init, opts);
    solCurvesB3{idx} = sol;
end

% zoomed out version
fig5 = figure(5);
set(fig5,'Position',[100,100,700,1000])
fig5 = tiledlayout(3,2,"TileSpacing","compact");
nexttile(1)
hold on;
for idx = 1:numel(solCurvesB1)
    sol = solCurvesB1{idx};
    plot(sol(:,1),sol(:,2)+sol(:,3),'-','Color',blue)
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
axis tight
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsB1.gammaL),'$: Trajectories'])
nexttile(2)
hold on;
for idx = 1:numel(solCurvesB1)
    sol = solCurvesB1{idx};
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
axis tight
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Endpoints only')

nexttile(3)
hold on;
for idx = 1:numel(solCurvesB2)
    sol = solCurvesB2{idx};
    plot(sol(:,1),sol(:,2)+sol(:,3),'-','Color',blue)
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
axis tight
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsB2.gammaL),'$: Trajectories'])
nexttile(4)
hold on;
for idx = 1:numel(solCurvesB2)
    sol = solCurvesB2{idx};
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
axis tight
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Endpoints only')

nexttile(5)
hold on;
for idx = 1:numel(solCurvesB3)
    sol = solCurvesB3{idx};
    plot(sol(:,1),sol(:,2)+sol(:,3),'-','Color',blue)
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
axis tight
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsB3.gammaL),'$: Trajectories'])
nexttile(6)
hold on;
for idx = 1:numel(solCurvesB3)
    sol = solCurvesB3{idx};
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
axis tight
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Endpoints only')
exportgraphics(fig5,'model03_phase_gammaL+bd+lp_zoomout.pdf',...
    'ContentType','vector')
% close(gcf)
% zoomed in version
fig6 = figure(6);
set(fig6,'Position',[100,100,700,1000])
fig6 = tiledlayout(3,2,"TileSpacing","compact");
nexttile(1)
hold on;
for idx = 1:numel(solCurvesB1)
    sol = solCurvesB1{idx};
    plot(sol(:,1),sol(:,2)+sol(:,3),'-','Color',blue)
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
xlim([0,max(initB)])
ylim([0,max(initL)])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsB1.gammaL),'$: Trajectories'])
nexttile(2)
hold on;
for idx = 1:numel(solCurvesB1)
    sol = solCurvesB1{idx};
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
xlim([0,max(initB)])
ylim([0,max(initL)])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Endpoints only')

nexttile(3)
hold on;
for idx = 1:numel(solCurvesB2)
    sol = solCurvesB2{idx};
    plot(sol(:,1),sol(:,2)+sol(:,3),'-','Color',blue)
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
xlim([0,max(initB)])
ylim([0,max(initL)])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsB2.gammaL),'$: Trajectories'])
nexttile(4)
hold on;
for idx = 1:numel(solCurvesB2)
    sol = solCurvesB2{idx};
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
xlim([0,max(initB)])
ylim([0,max(initL)])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Endpoints only')

nexttile(5)
hold on;
for idx = 1:numel(solCurvesB3)
    sol = solCurvesB3{idx};
    plot(sol(:,1),sol(:,2)+sol(:,3),'-','Color',blue)
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
xlim([0,max(initB)])
ylim([0,max(initL)])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title(['$\gamma_L = ',num2str(paramsB3.gammaL),'$: Trajectories'])
nexttile(6)
hold on;
for idx = 1:numel(solCurvesB3)
    sol = solCurvesB3{idx};
    plot(sol(end,1),sol(end,2)+sol(end,3),'.r')
end
hold off;
xlim([0,max(initB)])
ylim([0,max(initL)])
xlabel('Basal cell number')
ylabel('Luminal cell number')
title('Endpoints only')
exportgraphics(fig6,'model03_phase_gammaL+bd+lp_zoomin.pdf',...
    'ContentType','vector')
% close(gcf)
%% Plotting sample timecourses

sys1 = @(t,x) bilayerDuct_03(t,x,paramsA2,'nochange');
[t1,sol1] = ode45(sys1,[0,5e2],[1.5214;0.1;0],opts);
% q1 = usefulQuantities_03(sol1,paramsA2);
[t2,sol2] = ode45(sys1,[0,5e2],[0.1;2.9429;0],opts);
[t3,sol3] = ode45(sys1,[0,5e2],[4.5;8.5;0],opts);
[t4,sol4] = ode45(sys1,[0,5e2],[10;5;0],opts);
% q2 = usefulQuantities_03(sol2,paramsA2);
fig3 = figure(3);
set(fig3,'Position',[100,100,650,400])
fig3 = tiledlayout(1,2,'TileSpacing','compact');
nexttile(1)
hold on;
plot(t1,sol1(:,1),'--','DisplayName','Basal \#1','Color',blue)
plot(t1,sol1(:,2),'DisplayName','Luminal \#1','Color',blue)
plot(t2,sol2(:,1),'--','DisplayName','Basal \#2','Color',orange)
plot(t2,sol2(:,2),'DisplayName','Luminal \#2','Color',orange)
hold off;
ylim([0,5])
title(['$\gamma_L = ',num2str(paramsA2.gammaL),'$: Die out'])
xlabel('Time')
ylabel('Cell number')
legend
nexttile(2)
hold on;
plot(t3,sol3(:,1),'--','DisplayName','Basal \#3','Color',blue)
plot(t3,sol3(:,2),'DisplayName','Luminal \#3','Color',blue)
plot(t4,sol4(:,1),'--','DisplayName','Basal \#4','Color',orange)
plot(t4,sol4(:,2),'DisplayName','Luminal \#4','Color',orange)
hold off;
ylim([0,10])
title(['$\gamma_L = ',num2str(paramsA2.gammaL),'$: Stable population'])
xlabel('Time')
ylabel('Cell number')
legend
exportgraphics(fig3,'model03_gammaL=1_sampTimecourses.pdf', ...
    'ContentType','vector')

%% Hysteresis?

paramToChange = 'gammaL';
params.rBDbar = 0; % currently, hysteresis seems to require rBDbar = 0
params.tCutoff = 50;
tspan = [0,2*params.tCutoff];
params.afterCutoff = 2.5;
init = [1;1;0];

paramsA = params;
paramsA.gammaL = 0.25;
sysA = @(t,x) bilayerDuct_03(t,x,paramsA,paramToChange);

paramsB = params;
paramsB.gammaL = 1;
sysB = @(t,x) bilayerDuct_03(t,x,paramsB,paramToChange);

paramsC = params;
paramsC.gammaL = 2.5;
sysC = @(t,x) bilayerDuct_03(t,x,paramsC,paramToChange);

[tA, solA] = ode45(sysA, tspan, init, opts);
qA = usefulQuantities_02(solA,paramsA);
[tB, solB] = ode45(sysB, tspan, init, opts);
qB = usefulQuantities_02(solB,paramsB);
[tC, solC] = ode45(sysC, tspan, init, opts);
qC = usefulQuantities_02(solC,paramsC);
cutoffIdxs = [find(tA>paramsA.tCutoff,1);find(tB>paramsB.tCutoff,1)];

fig4 = figure(4);
set(fig4,'Position',[100, 100, 600, 600])
fig4 = tiledlayout(2,2,'TileSpacing','compact');
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
ylim([0,max(qA.rL)])
title('\bf C')
xlabel('Time')
ylabel('Luminal Cell Proliferation')
nexttile(4)
hold on;
plot(tA,qA.rB,'-','Color',green)
plot(tB,qB.rB,'-','Color',yellow)
plot(tC,qC.rB,'-','Color',blue)
hold off;
ylim([0,0.3])
title('\bf D')
xlabel('Time')
ylabel('Basal Cell Proliferation')
title(fig4,'Changing $\gamma_L$','Interpreter','latex')
lg = legend([l1,l2,l3,l4,l5,l6]);
lg.NumColumns = 3;
lg.Layout.Tile = 'south';
exportgraphics(fig4,'model03_gammaL_hysteresis.pdf','ContentType','vector')