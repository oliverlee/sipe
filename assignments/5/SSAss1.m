clear all
close all
clc

format long;
axfs = 18;
lgndfs = 16;

hankelmatrix = @(y, i, s, N) hankel(y(i+1:i+s, :), y(i+s:i+N+s-1, :));

%% A. create system & data, B. Data matrices, C. Column space Os + SVD
Nt = 1000;                                % No. Samples of simulation
dt = 0.01;                                % Sample time
t = (0:Nt-1).'*dt;                        % Time vector

% M=0.5;B=1;K=150;                      % System parameters
M = 0.1; B = 0.5; K = 450;
sys0 = c2d(ss(tf(1,[M B K])), dt);         % discrete system representation
A0 = sys0.a; B0 = sys0.b; C0 = sys0.c; D0 = sys0.d;

s = 3; ii = 0; N = 200; % Hankel matrix starting at y(ii) size s x N
r = size(B0, 2); % number of inputs
l = size(C0, 1); % number of outputs

u0 = zeros(Nt, r);                         % zero input
x0 = [2, 1; 1, 2; 1, 1];                % vector of initial conditions
rep = size(x0, 1);
lgnd = cell(rep, 1);

fig1 = figure(1); clf;
fig2 = figure(2); clf;
h2 = zeros([rep, 1]);
colset = cool(rep);

hankel_0_3_200 = @(y) hankelmatrix(y, ii, s, 200);
for i = 1:rep % loop for different initial conditions
    x0i = x0(i, :)';
    y0 = lsim(sys0, u0, t, x0i); % clean output
    lgnd{i} = sprintf('x0 = [%d, %d]''', x0i(1), x0i(2));

    YsN = hankel_0_3_200(y0);
    UsN = hankel_0_3_200(u0);

    [U, S, V] = svd(YsN, 'econ');           % Singular value decomposition
    Us = U*S;
    disp(sprintf('singular values S for x0 = [%d, %d]', x0i(1), x0i(2)));
    disp(S);

    figure(1);
    plot(t, y0, 'color', colset(i, :), 'linewidth', 2); hold on;

    figure(2);
    h2(i) = plot3(YsN(1,:), YsN(2,:), YsN(3,:), 'color', colset(i, :),...
        'linewidth', 2); hold on;
    for j = 1:size(S, 1)
        line([0, Us(1, j)], [0, Us(2, j)], [0, Us(3, j)],...
            'color', colset(i, :), 'linestyle', '--', 'linewidth', 10);
    end

    figure(3);
    semilogy(S, 'color', colset(i, :), 'marker', '*'); hold on;
end

figure(1);
set(fig1, 'name', 'Output time domain')
xlabel('t [s]', 'fontsize', axfs)
ylabel('y [m]', 'fontsize', axfs)
l1 = legend(lgnd);
set(l1, 'fontsize', lgndfs);

figure(2)
set(gcf,'name','Data Space')
xlabel('y(k) [m]', 'fontsize', axfs);
ylabel('y(k + 1) [m]', 'fontsize', axfs);
zlabel('y(k + 2) [m]', 'fontsize', axfs);
l2 = legend(h2, lgnd);
set(l2, 'fontsize', lgndfs);

figure(3)
xlim([0 s+1])
xlabel('singular value index', 'fontsize', axfs);
ylabel('singular value', 'fontsize', axfs);
l3 = legend(lgnd);
set(l3, 'fontsize', lgndfs);

eps_save('datamatrix1a', figure(2));

%% C
figure(3); clf
figure(4); clf
% System parameters
M = [0.1, 0.5, 0.1];
B = [0.5, 1, 5];
K = [450, 150, 450];
x0_c = x0(end, :);
colset = cool(numel(M));
lgnd = cell(numel(M), 1);

for i = 1:numel(M) % loop for different system parameters
    sys = c2d(ss(tf(1, [M(i), B(i), K(i)])), dt);
    y = lsim(sys, u0, t, x0_c);
    lgnd{i} = sprintf('M = %g, B = %g, K = %d', M(i), B(i), K(i));

    YsN = hankel_0_3_200(y);

    figure(4);
    plot3(YsN(1,:), YsN(2,:), YsN(3,:), 'color', colset(i, :),...
        'linewidth', 2); hold on;

    figure(3);
    [~, S_c, ~] = svd(YsN, 'econ');           % Singular value decomposition
    semilogy(1:length(S_c), diag(S_c), '*', 'color', colset(i, :));, hold on;
    disp(['singular values for ', lgnd{i}])
    S_c(1)
    S_c(5)
    S_c(end)
end

figure(3)
xlim([0 s+1])
xlabel('singular value index', 'fontsize', axfs);
ylabel('singular value', 'fontsize', axfs);
l3 = legend(lgnd);
set(l3, 'fontsize', lgndfs);

figure(4)
set(gcf,'name', sprintf('Data Space with x0 = [%d, %d]', x0_c(1), x0_c(2)));
xlabel('y(k) [m]', 'fontsize', axfs);
ylabel('y(k + 1) [m]', 'fontsize', axfs);
zlabel('y(k + 2) [m]', 'fontsize', axfs);
l2 = legend(lgnd);
set(l2, 'fontsize', lgndfs);

eps_save('svd1d', figure(3));
eps_save('datamatrix1c', figure(4));


%% D. Retrieve system from SVD
n = 2;                                 % System order
Un = U(:, 1:n);                        % Reduced output singular vectors
Vn = V(:, 1:n);                        % Reduced input singular vectors
Sn = S(1:n, 1:n);                      % Reduced singular value matrix
SnVtn = Sn*Vn';

Aid = Un(1:(s - 1)*l, :)\Un(l + 1:s*l, :);
Bid = SnVtn(:, 1:r);
Cid = Un(1:l, :);
Did = y0(1);
sysid = ss(Aid, Bid, Cid, Did, dt);       % Identified system out of SVD
eig0 = eig(A0);
eigid = eig(Aid);
disp('original system eigenvalues');
disp(eig0);
disp('identified system eigenvalues');
disp(eigid);


%% E. plot end result
x0_e = x0(end, :);
y0 = lsim(sys0, u0, t, x0_e);
yi = lsim(sysid, u0, t, x0_e);
x0T = SnVtn(:, 1); % note this uses the same x0
yiT = lsim(sysid, u0, t, x0T);

colset = cool(3);
figure(1); clf
set(fig1, 'name', 'Output time domain')
plot(t, y0, 'color', colset(1, :), 'linewidth', 2); hold on;
plot(t, yi, '--', 'color', colset(2, :), 'linewidth', 2);
plot(t, yiT, '--', 'color', colset(3, :), 'linewidth', 2); hold off;
xlabel('t [s]', 'fontsize', axfs)
ylabel('y [m]', 'fontsize', axfs)
l1 = legend('y0', 'yi', 'yi_T');
set(l1, 'fontsize', lgndfs);
xlim([-0.1, 2.0])
ylim([-0.5, 0.5])
eps_save('output1e', figure(1));


%% F. redo 3 plots with noise and compare to original
yk = y0 + 1e-2*randn(size(y0));           % output noise added

YsN = hankel_0_3_200(y0);
YksN = hankel_0_3_200(yk);

[U, S, V] = svd(YsN, 'econ');           % Singular value decomposition
[Uk, Sk, Vk] = svd(YksN, 'econ');           % Singular value decomposition

lgnd = {'y0', 'yk', 'y_{id}', 'yk_{id}'};
colset = cool(4);

figure(2); clf
plot3(YsN(1,:), YsN(2,:), YsN(3,:), 'color', colset(1, :),...
    'linewidth', 2); hold on;
plot3(YksN(1,:), YksN(2,:), YksN(3,:), 'color', colset(3, :),...
    'linewidth', 2); hold off;
set(gcf,'name','Data Space')
xlabel('y(k) [m]', 'fontsize', axfs);
ylabel('y(k + 1) [m]', 'fontsize', axfs);
zlabel('y(k + 2) [m]', 'fontsize', axfs);
l2 = legend(lgnd(1:2));
set(l2, 'fontsize', lgndfs);

figure(3); clf
semilogy(1:length(S), diag(S), '*', 'color', colset(1, :));, hold on;
semilogy(1:length(Sk), diag(Sk), '*', 'color', colset(3, :));, hold off;
xlim([0, s + 1])
set(gcf, 'name', 'singular values');
xlabel('singular value index', 'fontsize', axfs);
ylabel('singular value', 'fontsize', axfs);
l3 = legend(lgnd(1:2));
set(l3, 'fontsize', lgndfs);

[Aid, Bid, Cid, Did, x0id] = estimatesystem(YsN, UsN, s, n);
[Akid, Bkid, Ckid, Dkid, x0kid] = estimatesystem(YksN, UsN, s, n);

yid = lsim(sysid, u0, t, x0id);
ykid = lsim(sysid, u0, t, x0kid);

figure(1); clf
set(fig1, 'name', 'Output time domain')
plot(t, y0, 'color', colset(1, :), 'linewidth', 2); hold on;
plot(t, yk, 'color', colset(3, :), 'linewidth', 2);
plot(t, yid, 'color', colset(2, :), 'linewidth', 2);
plot(t, ykid, 'color', colset(4, :), 'linewidth', 2); hold off;
xlabel('t [s]', 'fontsize', axfs)
ylabel('y [m]', 'fontsize', axfs)
l1 = legend(lgnd);
set(l1, 'fontsize', lgndfs);
xlim([-0.1, 2.0])
ylim([-0.5, 0.5])

eps_save('output1f', figure(1));
eps_save('datamatrix1f', figure(2));
eps_save('svd1f', figure(3));

%% G. increase s to reduce effect of noise
%Narray = N:100:size(yk, 1) - s;
%Sarray = zeros(length(Narray), s);
%lgnd = cell(length(Narray), 1);
%colset = cool(length(Narray));
sarray = s:1:10;
Sarray = cell(length(sarray), 1);
lgnd = cell(length(sarray), 1);
colset = cool(length(sarray));

figure(3); clf;
%for i = 1:length(Narray)
for i = 1:length(sarray)
    YsiN = hankelmatrix(yk, ii, sarray(i), N);
    [~, Si, ~] = svd(YsiN, 'econ');
    Sarray{i} = diag(Si)';
    semilogy(1:sarray(i), diag(Si), '*', 'color', colset(i, :)); hold on;
    %lgnd{i} = sprintf('N = %d', Narray(i))
    lgnd{i} = sprintf('s = %d', sarray(i))
end
hold off;
xlim([0, sarray(end) + 1])
set(gcf, 'name', 'singular values');
xlabel('singular value index', 'fontsize', axfs);
ylabel('singular value', 'fontsize', axfs);
l3 = legend(lgnd);
set(l3, 'fontsize', lgndfs);

%% H. fourth order system
% m1*x1dd = k2*(x2 - x1) - k1*x1 + c2*(x2d - x1d) - c1*x1d
% m2*x2dd = - k2*(x2 - x1) - c2*(x2d - x1d)

m1 = M(1); m2 = M(1);
c1 = B(1); c2 = B(1);
k1 = K(1); k2 = K(1);
A1 = [-(c1 + c2)/m1,  c2/m1, -(k1 + k2)/m1,  k2/m1;
              c2/m2, -c2/m2,         k2/m2, -k2/m2;
                  1,      0,              0,     0;
                  0,      1,              0,     0];

disp('original system eigenvalues');
%disp(sort(eig(A1)));
B1 = [0, 1/m2, 0, 0]';
%C1 = [0, 0, 1, 0;
%      0, 0, 0, 1];
C1 = [0, 0, 1, 0];
D1 = 0;
sys4 = c2d(ss(A1, B1, C1, D1), dt);         % discrete system representation
disp(sort(eig(sys4.a)));
u40 = zeros(Nt, 1);                         % zero input
x40 = [1, 1, 1, 1]';
y40 = lsim(sys4, u40, t, x40);
s = 5; ii = 0; N = 200; % Hankel matrix starting at yk(ii) size s x N
YsN = hankelmatrix(y40, 0, s, N);
UsN = hankelmatrix(u40, 0, s, N);
[U,S,V] = svd(YsN, 'econ');           % Singular value decomposition

n = 4;
%l = size(y40, 2);                           % number of outputs
%r = size(u40, 2);                           % number of inputs
%Un = U(:, 1:n);                        % Reduced output singular vectors
%Vn = V(:, 1:n);                        % Reduced input singular vectors
%Sn = S(1:n, 1:n);                      % Reduced singular value matrix
%SnVtn = Sn*Vn';
[A4id, B4id, C4id, D4id] = estimatesystem(YsN, UsN, s, n);
disp('identified system eigenvalues');
disp(sort(eig(A4id)));
sysid = ss(A4id, B4id, C4id, D4id, dt);       % Identified system out of SVD
