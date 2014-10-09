clear all
close all
clc

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
Syy = Y.*conj(Y)/length(u);

H_ol = Syu ./ Suu;
%C_ol = Syu ./ sqrt(Suu.*Syu);
C_ol = abs(Syu).^2 ./ (Suu.*Syy);

cmap = flipud(hsv(2));
plotspectrum(f, H_true, f, H_ol, C_ol, cmap, {'True', 'Raw'});
eps_save('question2c')

%% Question 2d
segment = 1:10;
wf = cell(1, length(segment));
wSuu = cell(1, length(segment));
wSyu = cell(1, length(segment));
wSyy = cell(1, length(segment));
wH_ol = cell(1, length(segment));
wC_ol = cell(1, length(segment));
for i = segment
    wf{i}= f(1:i:floor(end/i)*i);
    wSuu{i} = welchspectrum(u, u, i);
    wSyu{i} = welchspectrum(y, u, i);
    wSyy{i} = welchspectrum(y, y, i);
    wH_ol{i} = wSyu{i} ./ wSuu{i};
    %wC_ol{i} = wSyu{i} ./ sqrt(wSuu{i}.*wSyy{i});
    wC_ol{i} = abs(wSyu{i}).^2 ./ (wSuu{i}.*wSyy{i});
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
% set variance of input signal u(t) to a high value
rng(1234) % set seed for generation of white noise input
u_f = normrnd(0, 100*u_var, size(t));
rng(876348) % set seed for generation of additive output noise
y_f = lsim(H, u, t) + normrnd(0, n_var, size(t));

% set variance of additive output noise n(t) to a high value
rng(1234) % set seed for generation of white noise input
u_g = normrnd(0, u_var, size(t));
rng(876348) % set seed for generation of additive output noise
y_g = lsim(H, u, t) + normrnd(0, 100*n_var, size(t));

% what do we change the variance to in part 2g? should welch averaging be used
% to calculate the frequency response H and coherence?
f_fg = {wf{end}, wf{end}};
Suu_f = welchspectrum(u_f, u_f, 10);
Syu_f = welchspectrum(y_f, u_f, 10);
Syy_f = welchspectrum(y_f, y_f, 10);
Suu_g = welchspectrum(u_g, u_g, 10);
Syu_g = welchspectrum(y_g, u_g, 10);
Syy_g = welchspectrum(y_g, y_g, 10);
H_ol_fg = {Syu_f ./ Suu_f, Syu_g ./ Suu_g};
C_ol_fg = {abs(Syu_f).^2 ./ (Suu_f .* Syy_f),...
           abs(Syu_g).^2 ./ (Suu_g .* Syy_g)};

cmap = flipud(hsv(3));
legendstr = cell(1, 3);
legendstr = {'True',...
    sprintf('\\sigma_u^2 = %d', 100*u_var),...
    sprintf('\\sigma_n^2 = %d', 100*n_var)};
plotspectrum(f, H_true, f_fg, H_ol_fg, C_ol_fg, cmap, legendstr)
eps_save('question2f')
