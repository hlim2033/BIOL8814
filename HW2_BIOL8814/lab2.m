%% NUM1
pars.a = 5;
pars.b = 1;
pars.r = 30;
x0 = [0 0];
t0 = 0;
tf = 1;
[t, x] = ode45(@proteinproductino, [t0 tf], x0, [], pars);
figure
plot(t, x, 'LineWidth', 3)
xlabel('time (hr)', 'FontSize', 20)
ylabel('nMol', 'FontSize', 20)
legend('mRNA', 'Protein')
set(gca, 'FontSize', 20)

pars.a = 5;
pars.b = 1;
pars.r = 30;
x0 = [0 0];
t0 = 0;
tf = 24;
[t, x] = ode45(@proteinproductino, [t0 tf], x0, [], pars);
figure
plot(t, x, 'LineWidth', 3)
xlabel('time (hr)', 'FontSize', 20)
ylabel('nMol', 'FontSize', 20)
legend('mRNA', 'Protein')
set(gca, 'FontSize', 20)

%% NUM2
mnullcline = pars.a/pars.b ; 
pnullcine = @(p) pars.b.*p./pars.r;

pvec = 0:20:200;
mvec = 0:1:8;

plot(mnullcline*ones(size(pvec)), pvec, 'LineWidth', 3)
hold on
plot(pnullcine(pvec),pvec,'Linewidth',3)
hold off

legend('dmdt = 0', 'dpdt = 0')
xlim([min(mvec) max(mvec)])
ylim([min(pvec) max(pvec)])
xlabel('mRNA Conc (nMol)', 'FontSize' , 20)
ylabel('Protein Conc (nMol)', 'FontSize', 20)
set(gca, 'FontSize', 20)

%%
dmdt = @(m) pars.a - pars.b.*m
dpdt = @(m,p) pars.r.*m - pars.b.*p;
[MM,PP] = meshgrid(mvec,pvec)
hold on
quiver(MM,PP,dmdt(MM), dpdt(MM,PP), 'k', 'MaxHeadSize', 0.015,'LineWidth', 2);
plot(mnullcline*ones(size(pvec)), pvec, 'LineWidth', 3)
plot(pnullcine(pvec),pvec,'Linewidth',3)
xlim([min(mvec) max(mvec)])
ylim([min(pvec) max(pvec)])
xlabel('mRNA Conc (nMol)', 'FontSize' , 20)
ylabel('Protein Conc (nMol)', 'FontSize', 20)
hold off


%%
m0_vec = [0 1 7];
p0_vec = [0 150 25];
t0 = 0;
tf =24;
for i=1:3
    x0 = [m0_vec(i) p0_vec(i)];
    [t, x] = ode45( @proteinproductino, [t0 tf], x0, [], pars);
    
    subplot(3,1,i)
    %tmph = plot(x(:,1),x(:,2),'LineWidth',3);
    hold on 
    plot(t, x(:,1), 'LineWidth', 3)
    plot(t, x(:,2), 'LineWidth', 3)
    hold off
    %set(tmph, 'color','r')
end

%%
m0_vec = [0 1 7];
p0_vec = [0 150 25];
t0 = 0;
tf =24;
for i=1:3
    x0 = [m0_vec(i) p0_vec(i)];
    [t, x] = ode45( @proteinproductino, [t0 tf], x0, [], pars);
    
    fig1 = plot(x(:,1),x(:,2),'LineWidth',3);
    hold on
    
    set(fig1, 'color','r')
    
end
quiver(MM,PP,dmdt(MM), dpdt(MM,PP), 'k', 'MaxHeadSize', 0.015,'LineWidth', 2);
xlim([min(mvec) max(mvec)])
ylim([min(pvec) max(pvec)])

%% NUM3
allnull = @(p) pars.r.*pars.a/pars.b - pars.b*p;
equil_num = fzero(allnull, 100);

%%
pstar = pars.r*pars.a/pars.b^2
mstar = pars.a/pars.b

pperturb = pstar.*.01;
mperturb = mstar.*.01;

x0 = [mstar + mperturb pstar+pperturb];
x1 = [mstar - mperturb pstar-pperturb];
t0 = 0;
tf = 10;

[t, x] = ode45(@proteinproductino, [t0 tf], x0, [], pars)
[t1,x1] = ode45(@proteinproductino, [t0 tf], x1, [], pars)
figure
hold on 
plot(t, x(:,2), 'LineWidth', 3)
plot(t1, x1(:,2), 'LineWidth', 3)
xlabel('time (hr)', 'FontSize', 20)
ylabel('Concentration (nMol)', 'FontSize', 20)
legend('Protein (starting above)', 'Protein (starting below)')
set(gca, 'FontSize', 20)
%% NUM3
avals = [5, 5, 10];
bvals = [1, 0.5 ,1];
rvals = [30, 30, 30];
for i= 1:length(avals)
    pars.a = avals(i);
    pars.b = bvals(i);
    pars.r = rvals(i);
    a = pars.a; b= pars.b; r=pars.r;
    mstar = a/b;
    pstar = r*a/b^2;
    x0 = [mstar*1.01 pstar*1.02];
    t0 = 0;
    tf= 10;
    [t, x] = ode45(@proteinproductino, [t0:0.2:tf], x0, [], pars);
    u = abs(x(:,1) - mstar);
    v = abs(x(:,2) - pstar);
    figure
    tmph = semilogy(t, [u v], 'o')
    set(tmph, 'markersize', 8, 'markerfacecolor', [0.75 0.75 0.75]);
    set(tmph(1), 'marker', 'o')
    set(tmph(2), 'marker', 'd')
    xlabel('time (hr)', 'FontSize', 20, 'interpreter', 'latex')
    ylabel('Perturbed concentration (nM)', 'FontSize', 20, 'interpreter', 'latex')
    legend('u', 'v');
    set(gca, 'fontsize', 20)
    xlim([0,5])
    hold on
    tmph = semilogy(t, exp(-b*t), 'r-');
    hold off
end

%% NUM4
bvec = 0.01:0.01:0.2;
mstarvec = zeros(1, length(bvec));
pstarvec = zeros(1, length(bvec));
for jj = 1:length(bvec)
    pars.b = bvec(jj)
    pstarvec(jj) = pars.r*pars.a/pars.b^2;
    mstarvec(jj) = pars.a/pars.b;
end
figure(1)
plot(bvec, mstarvec, '.', 'MarkerSize', 15)
xlabel('b','FontSize',20)
ylabel('mRNA Conc (nMol)', 'FontSize',20)
set(gca, 'FontSize',20)


%% NUM5
pars.a = 5;
pars.b = 1;
pars.r = 30;
x0 = [0 0];
t0 = 0;
tf = 1;
[t, x] = ode45(@proteinproductino, [t0 tf], x0, [], pars);
tmph = plot(t,x)

tmpl = legend('mRNA', 'Protein');
set(tmpl, 'Location', 'NorthWest')
legend('boxoff')

xlabel('time (hr)', 'FontSize', 20)
ylabel('Concentration (nMol)', 'FontSize', 20)

set(gca, 'FontSize', 20)
set(tmph, 'Marker', '.','MarkerSize', 20);
