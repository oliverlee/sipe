%% Assignment 2: Open & Closed Loop Identification
%  [WB2301] SIPE 2013-2014

clear all
close all
clc

%% Question 1
fs = 1000; % Hz
T = 10; % seconds
dt = 1/fs;
t = 0:dt:T-dt;
f_sig = 3; % Hz
u = 2*sin(2*pi*f_sig*t) + randn(size(t));

f_butter_cutoff = 25; % Hz
[B, A] = butter(3, f_butter_cutoff/(fs/2));
y1000 = filter(B, A, u);
y500 = y1000(2:2:end);
y250 = y1000(4:4:end);
y100 = y1000(10:10:end);

%% Question 1a
Y1000 = fft(y1000) / length(y1000);
f1000 = (0:length(y1000) - 1) * fs/length(y1000);
figure('Name', 'Question 1a')
subplot(2, 1, 1)
    plot(f1000, mag2db(abs(Y1000)), 'r', 'linewidth', 2)
    xlabel('frequency [Hz]', 'Fontsize', 18);
    ylabel('magnitude(Y(f)) [dB]', 'Fontsize', 18);
    set(gca, 'FontSize', 18)
subplot(2, 1, 2)
    plot(f1000, 180/pi*phase(Y1000), 'r', 'linewidth', 2)
    xlabel('frequency [Hz]', 'Fontsize', 18);
    ylabel('phase(Y(f)) [deg]', 'Fontsize', 18);
    set(gca, 'FontSize', 18)
eps_save('question1a')

%% Question 1b
Y500 = fft(y500) / length(y500);
Y250 = fft(y250) / length(y250);
Y100 = fft(y100) / length(y100);
f500 = f1000(1:end/2);
f250 = f1000(1:end/4);
f100 = f1000(1:end/10);

cmap = hsv(4);
figure('Name', 'Question 1b'); clf;
    semilogx(f1000, mag2db(abs(Y1000)), 'color', cmap(1,:), 'linewidth', 2)
    hold on; % hold after first plot or axes are not set correctly
    semilogx(f500, mag2db(abs(Y500)), 'color', cmap(2,:), 'linewidth', 2)
    semilogx(f250, mag2db(abs(Y250)), 'color', cmap(3,:), 'linewidth', 2)
    semilogx(f100, mag2db(abs(Y100)), 'color', cmap(4,:), 'linewidth', 2)
    hold off;
    xlabel('frequency [Hz]', 'Fontsize', 18);
    ylabel('magnitude(Y(f)) [dB]', 'Fontsize', 18);
    set(gca, 'FontSize', 18)
    l = legend({'1000 Hz', '500 Hz', '250 Hz', '100 Hz'});
    set(l, 'FontSize', 18);
eps_save('question1b')

%% Question 1c
y1000_5 = y1000(1:length(y1000) * 5/T);
Y1000 = fft(y1000);
Y1000_5 = fft(y1000_5);
f1000_5 = f1000(1:2:end);

cmap = hsv(2);
figure('Name', 'Question 1c'); clf;
    plot(f1000, mag2db(abs(Y1000)), 'color', cmap(1,:),'linewidth', 2)
    hold on;
    plot(f1000_5, mag2db(abs(Y1000_5)), 'color', cmap(2,:),'linewidth', 2)
    hold off;
    xlabel('frequency [Hz]', 'Fontsize', 18);
    ylabel('magnitude(Y(f)) [dB]', 'Fontsize', 18);
    set(gca, 'FontSize', 18)
    l = legend({sprintf('T = %d s', T), 'T = 5 s'});
    set(l, 'FontSize', 18);
eps_save('question1c')

%% Question 2

% Some administration
dt = 0.01;                 % Sample time
T = 50;                   % Total time
N = T/dt;                 % number of samples
t = (0:N-1).'*dt;         % Time vector
Fs = 1/dt;                 % Sample frequency
f = (0:N-1).'/T;          % Frequency vector (double sided)

%options=simset('OutputPoints','specified'); % only produce output for the time values in vector <t>
%load_system('Ass2_System_Open');
%
%%default values
%set_param('Ass2_System_Open/Input u: white noise','Variance','0.5');
%set_param('Ass2_System_Open/Noise n','Variance','0.1');
%
% Simulate model
%sim('Ass2_System_Open',t,options)
%clear tout

H = tf(1, [0.01, 0.03, 1]);
u_var = 0.5;
n_var = 0.1;
rng(1234) % set seed for generation of white noise input
u = normrnd(0, u_var, size(t));
rng(876348) % set seed for generation of additive output noise
y = lsim(H, u, t) + normrnd(0, n_var, size(t));

% Plot input and output
figure(100)
subplot(211)
    plot(t, u, 'r')
    xlabel('time [s]')
    ylabel('Input')
    title ('Input and output signals')
subplot(212)
    plot(t, y, 'r')
    ylabel('Output')
    xlabel('time [s]')

%% Question 2c

% Define the true system
H_true = squeeze(freqresp(H, f, 'Hz'));

% Estimate Open loop transfer function (H_ol) & Coherence (C_ol)
U = fft(u);
Y = fft(y);
Suu = U.*conj(U)/length(u);
Syu = Y.*conj(U)/length(u);

H_ol = Syu ./ Suu;
C_ol = Syu ./ sqrt(Suu.*Syu);

cmap = flipud(hsv(2));
plotspectrum(f, H_true, f, H_ol, C_ol, cmap, {'True', 'Raw'});
eps_save('question2c')

%% Question 2d
segment = 1:10;
wf = cell(1, length(segment));
wSuu = cell(1, length(segment));
wSyu = cell(1, length(segment));
wH_ol = cell(1, length(segment));
wC_ol = cell(1, length(segment));
for i = segment
    wf{i}= f(1:i:floor(end/i)*i);
    wSuu{i} = welchspectrum(u, u, i);
    wSyu{i} = welchspectrum(y, u, i);
    wH_ol{i} = wSyu{i} ./ wSuu{i};
    wC_ol{i} = wSyu{i} ./ sqrt(wSuu{i}.*wSyu{i});
end

cmap = jet(length(segment) + 1);
legendstr = cell(1, length(segment) + 1);
legendstr{1} = 'True';
for i = segment
    legendstr{i + 1} = sprintf('N = %d', segment(i));
end
plotspectrum(f, H_true, wf, wH_ol, wC_ol, cmap, legendstr)
eps_save('question2d')

%% Question 2f, 2g
high_var = 5;

% set variance of input signal u(t) to a high value
rng(1234) % set seed for generation of white noise input
u_f = normrnd(0, high_var, size(t));
rng(876348) % set seed for generation of additive output noise
y_f = lsim(H, u, t) + normrnd(0, n_var, size(t));

% set variance of additive output noise n(t) to a high value
rng(1234) % set seed for generation of white noise input
u_g = normrnd(0, u_var, size(t));
rng(876348) % set seed for generation of additive output noise
y_g = lsim(H, u, t) + normrnd(0, high_var, size(t));

% what do we change the variance to in part 2g? should welch averaging be used
% to calculate the frequency response H and coherence?
f_fg = {wf{end}, wf{end}};
Suu_f = welchspectrum(u_f, u_f, 10);
Syu_f = welchspectrum(y_f, u_f, 10);
Suu_g = welchspectrum(u_g, u_g, 10);
Syu_g = welchspectrum(y_g, u_g, 10);
H_ol_fg = {Syu_f ./ Suu_f, Syu_g ./ Suu_g};
C_ol_fg = {Syu_f ./ sqrt(Suu_f .* Syu_f), Syu_g ./ sqrt(Suu_g .* Syu_g)};

cmap = flipud(hsv(3));
legendstr = cell(1, 3);
legendstr = {'True', '\sigma_u^2 = 5', '\sigma_n^2 = 5'};
plotspectrum(f, H_true, f_fg, H_ol_fg, C_ol_fg, cmap, legendstr)
eps_save('question2f')
