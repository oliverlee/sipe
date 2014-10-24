clear all
close all
clc

%% question 1
fs = 100; % sample frequency in Hz
dt = 1/fs; % sample time in seconds
T = 50; % total sample time in seconds
N = T/dt; % number of samples
t = (0:N - 1).'*dt; % time vector
f = (0:N - 1).'/T; % frequency vector (double sided)

H = tf(1, [0.01, 0.03, 1]);
H.InputName = 'u';
H.OutputName = 'H_out';
G = tf(1, [1, 0.1]);
G.InputName = 'y';
G.OutputName = 'G_out';
SumU = sumblk('u = r - G_out');
SumY = sumblk('y = n + H_out');
sys = connect(H, G, SumU, SumY, {'r', 'n'}, {'u', 'y'});

rng(1262) % set seed for generation of white noise input
r_var = 1;
r = normrnd(0, r_var, size(t));

rng(1465) % set seed for generation of white noise input
n_var = 0.1;
n = normrnd(0, n_var, size(t));

output = lsim(sys, [r, n], t);
u = output(:, 1);
y = output(:, 2);

dat_ur = iddata(u, r, dt);
dat_yr = iddata(y, r, dt);
n = 4; % order is 4

% two stage method % u, r, ARMAX % y, u, FIR
s11 = armax(dat_ur, [n, n, n, 0]);
up = lsim(s11, r, t);
dat_yup = iddata(y, up, dt);
s12 = oe(dat_yup, [n, n, 0]);

% coprime method % y, r, ARX % u, r, ARMAX
s21 = armax(dat_yr, [n, n, n, 0]);
s22 = armax(dat_ur, [n, n, n, 0]);

f_half = f(1:ceil(length(f)/2));
H_true = squeeze(freqresp(H, f_half, 'Hz'));
H_2stage = squeeze(freqresp(s12, f_half, 'Hz'));
H_coprime = squeeze(freqresp(s21/s22, f_half, 'Hz'));
lgd = {'H_{True}', 'H_{2stage}', 'H_{coprime}'};

cmap = hsv(3);
figure()
subplot(2, 1, 1)
    semilogx(f_half, mag2db(abs(H_true)), '--',...
             'color', cmap(1, :), 'linewidth', 2);
    hold on;
    semilogx(f_half, mag2db(abs(H_2stage)),...
             'color', cmap(2, :), 'linewidth', 2);
    semilogx(f_half, mag2db(abs(H_coprime)),...
             'color', cmap(3, :), 'linewidth', 2);
    hold off;
    ylabel('magnitude(H(f)) [dB]', 'Fontsize', 18);
    xlim([0.01 f_half(end)])
subplot(2, 1, 2)
    semilogx(f_half, 180/pi*phase(H_true), '--',...
             'color', cmap(1, :), 'linewidth', 2);
    hold on;
    semilogx(f_half, 180/pi*phase(H_2stage),...
             'color', cmap(2, :), 'linewidth', 2);
    semilogx(f_half, 180/pi*phase(H_coprime),...
             'color', cmap(3, :), 'linewidth', 2);
    hold off;
    l = legend(lgd, 'location', 'southwest');
    set(l, 'FontSize', 18);
    ylabel('phase(H(f)) [deg]', 'Fontsize', 18);
    ylim([-210  30]);
    xlim([0.01 f_half(end)]);
eps_save('question1b')
