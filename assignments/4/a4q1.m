%% SIPE (2014-2015): Assignment 4

clear all
close all
clc

%% Load Dataset

global t y

% Data1
load('Data1');   % output y; timevector t

% Plot data
Timefit = figure; clf;
plot(t,y,'r','linewidth',4); hold on

%% Grid search (1)

% Start Time
tic

% Upper & Lower Bounds [a,b]
lb = [0 0];
ub = [10 10];

% Resolution
ares = lb(1):0.1:ub(1);
bres = lb(2):0.1:ub(2);

% Error & Parameter Calculation
e1 = zeros(length(ares),length(bres));
for qq = 1:length(ares)
    for ss = 1:length(bres)
        par = [ares(qq) bres(ss)];
        e1(qq,ss) = errfun(par);
    end
end

% Error plot (watch out use e transposed)
errsurf = figure; clf;
surfc(ares,bres,e1.'); hold on;
alpha(0.3);
xlabel('a', 'fontsize', 16);
ylabel('b', 'fontsize', 16);
zlabel('e', 'fontsize', 16);

% Find Optimum (parest & error residual)
[err1, index] = min(e1(:));
[i, j] = ind2sub(size(e1), index);
amin = ares(i);
bmin = ares(j);
parest1 = [amin bmin];
%clear tmp

% Estimate & Plot yest
[~,yest1] = errfun(parest1);

figure(Timefit);
plot(t,yest1,'b--','linewidth',2); hold on;
legend('Measurements','Grid Search')

% Stop Time
time1 = toc;

%% Gradient Search (2)

% Start Time
tic

% Gradient Search Options (Use 'doc optimset' for more options)
options = optimset('lsqnonlin');
options = optimset(options, 'Display','iter');
options = optimset(options,'outputfcn',@storepathway);

% Set pathway variables
global pathway
pathway.x = []; pathway.y = []; pathway.z = [];

% Initial Parameter Value (First step)
par0 = [5 5];

% LSQnonlin
[parest2,~,e2,~,~,~,J]=lsqnonlin(@errfun2,par0,lb,ub,options);
err2 = e2*e2';  % sum of squares (=sum(e2.^2))

% Plot Pathway
num_pathways = 10;
cmap = hsv(num_pathways);
h = zeros(1, num_pathways);

figure(errsurf)
h(1) = plot3(pathway.x, pathway.y, pathway.z,...
    '.-', 'color', cmap(1, :), 'linewidth', 4);
plot3(pathway.x(1), pathway.y(1), pathway.z(1),...
    '.', 'color', cmap(1, :), 'markersize', 25);
plot3(pathway.x(end), pathway.y(end), pathway.z(end),...
    '.', 'color', cmap(1, :), 'markersize', 50);

% Estimate & Plot yest
[~,yest2] = errfun2(parest2);

figure(Timefit);
plot(t,yest2,'g--','linewidth',2); hold on;
legend('Measurements','Grid Search','Gradient Search')

% Stop Time
time2 = toc;

% run gradient search with different initial conditions
figure(errsurf);
lgndstr = cell([1, num_pathways]);
lgndstr{1} = sprintf('a = 5, b = 5');
pari = rand([num_pathways - 1, 2]);
mid = (ub - lb)/2;
pari = [pari(:, 1)*mid(1) + mid(1),...
        pari(:, 2)*mid(2) + mid(2)];

for i = 2:num_pathways
    lgndstr{i} = sprintf('a = %g, b = %g', pari(i - 1, :));
    pathway.x = []; pathway.y = []; pathway.z = [];
    pest = lsqnonlin(@errfun2,pari(i - 1, :),lb,ub,options);
    h(i) = plot3(pathway.x, pathway.y, pathway.z,...
        '.-', 'color', cmap(i, :), 'linewidth', 4);
    plot3(pathway.x(1), pathway.y(1), pathway.z(1),...
        '.', 'color', cmap(i, :), 'markersize', 25);
    plot3(pathway.x(end), pathway.y(end), pathway.z(end),...
        '.', 'color', cmap(i, :), 'markersize', 50);
end
%l = legend(h, lgndstr);
%set(l, 'fontsize', 14);
eps_save('optsurf')

%%
[par, e, iter, time, h] = combinedsearch(0:10, 0:10, true);
eps_save('combsearch')

%% Generic Search (3)

% Start Time
tic

% Generic Algorithm Search (Check what they do with 'help gaoptimset')
options = gaoptimset;
options = gaoptimset(options,'populationsize', 50);
options = gaoptimset(options,'elitecount', 15);
options = gaoptimset(options,'crossoverfraction', 0.8);
options = gaoptimset(options,'popinitrange', [0; 8]);
options = gaoptimset(options,'plotfcn', @customplotfcn);
options = gaoptimset(options,'plotinterval', 1);

% Error Contour plot (watch out use e transposed)
figure(99); clf;
contour(ares,bres,e1.'); hold on;
daspect([1 1 1])

% Generic Algorithm
[parest3, err3] = ga(@errfun, 2,[],[],[],[],lb,ub,[],options);

% Estimate & Plot yest
[~,yest3] = errfun(parest3);

figure(Timefit);
plot(t,yest3,'k--','linewidth',2); hold on;
legend('Measurements','Grid Search','Gradient Search','Generic Algorithm')

% Stop Time
time3 = toc;
