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
    plot(f1000, phase(Y1000), 'r', 'linewidth', 2)
    xlabel('frequency [Hz]', 'Fontsize', 18);
    ylabel('phase(Y(f)) [rad]', 'Fontsize', 18);
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
    legend({'1000 Hz', '500 Hz', '250 Hz', '100 Hz'});
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
    legend({sprintf('T = %d s', T), 'T = 5 s'});
eps_save('question1c')

%% Question 2

% Some administration
dt  = 0.01;                 % Sample time
T   = 50;                   % Total time
N   = T/dt;                 % number of samples
t   = (0:N-1).'*dt;         % Time vector

Fs  = 1/dt;                 % Sample frequency
f   = (0:N-1).'/T;          % Frequency vector (double sided)
options=simset('OutputPoints','specified'); % only produce output for the time values in vector <t>
load_system('Ass2_System_Open');

%default values
set_param('Ass2_System_Open/Input u: white noise','Variance','0.5');
set_param('Ass2_System_Open/Noise n','Variance','0.1');

% Simulate model
sim('Ass2_System_Open',t,options)
clear tout

% Plot input and output
figure(100)
subplot(211)
    plot(t,u)
    xlabel('time [s]')
    ylabel('Input')
    title ('Input and output signals')
subplot(212)
    plot(t,y)
    ylabel('Output')
    xlabel('time [s]')

%% Question 2c

% Define the true system
    !! Fill in your own lines !!
H_true = [];


% Estimate Open loop transfer function (H_ol) & Coherence (C_ol)
    !! Fill in your own lines !!
H_ol = [];
C_ol = [];

% Plot
fmax = Fs/2; % Highest frequency in the plot

figure('Name','Question 1c')
subplot(311)
    loglog(f, abs(H_true), 'k--', 'linewidth', 2); hold on;
    loglog(f, abs(H_ol), 'b', 'linewidth', 2);
    title('Open-loop spectral estimation of H')
    legend({'True','Raw'})
    xlim([0.01 fmax])
subplot(312)
    semilogx(f, 180/pi*angle(H_true), 'k--', 'linewidth', 2); hold on;
    semilogx(f, 180/pi*angle(H_ol), 'b', 'linewidth', 2); 
    ylim([-180  180]); xlim([0.01 fmax]);
subplot(313)
    semilogx(f, abs(C_ol), 'b', 'linewidth', 2); hold on;
    ylim([0 1.1]); xlim([0.01 fmax]);
