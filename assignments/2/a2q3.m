clear all
close all
clc

%% Question 3
% Some administration
dt = 0.01;                 % Sample time
T = 500;                   % Total time
N = T/dt;                 % number of samples
t = (0:N-1).'*dt;         % Time vector
Fs = 1/dt;                 % Sample frequency
f = (0:N-1).'/T;          % Frequency vector (double sided)

H = tf(1, [0.01, 0.03, 1]);
H.InputName = 'u';
H.OutputName = 'H_out';
G = tf(1, [1, 0.1]);
G.InputName = 'y';
G.OutputName = 'G_out';
SumU = sumblk('u = r - G_out');
SumY = sumblk('y = n + H_out');
sys = connect(H, G, SumU, SumY, {'r', 'n'}, {'u', 'y'});

nsegment = 30;
f_cl= f(1:nsegment:floor(end/nsegment)*nsegment);

%% 3a
rng(1262) % set seed for generation of white noise input
r_var = 1;
r = normrnd(0, r_var, size(t));

rng(1465) % set seed for generation of white noise input
n_var = 0.1;
n = normrnd(0, n_var, size(t));

output = lsim(sys, [r, n], t);
u = output(:, 1);
y = output(:, 2);

Suu = welchspectrum(u, u, nsegment);
Syu = welchspectrum(y, u, nsegment);
Syy = welchspectrum(y, y, nsegment);
H_cl = Syu ./ Suu;
C_cl = abs(Syu).^2 ./ (Suu .* Syy);
H_true = squeeze(freqresp(H, f, 'Hz'));

cmap = flipud(hsv(2));
plotspectrum(f, H_true, f_cl, H_cl, C_cl, cmap, {'True', 'Open Loop'});
eps_save('question3a')

%% 3b
rng(1465) % set seed for generation of white noise input
n_var_scale = [0.01, 1, 10, 100, 1000];
n = cell(size(n_var_scale));
H_i = cell(size(n_var_scale));
C_i = cell(size(n_var_scale));
lgd = cell(1, 2 + length(n_var_scale));
lgd{1} = 'H_{True}';
lgd{2} = '-G^{-1}_{True}';

for i = 1:length(n_var_scale)
    n{i} = normrnd(0, n_var_scale(i)*n_var, size(t));
    output = lsim(sys, [r, n{i}], t);
    u = output(:, 1);
    y = output(:, 2);
    [H_i{i}, C_i{i}, f_out] = estimatesystem(y, u, 30, f);
    lgd{i + 2} = sprintf('H_{ol} (\\sigma_n = %g)', n_var_scale(i)*n_var);
end

H_true = squeeze(freqresp(H, f, 'Hz'));
minus_Ginv_true = squeeze(freqresp(-1/G, f, 'Hz'));

cmap = flipud(hsv(length(lgd)));
fmax = f(round(end/2));
figure()
subplot(3, 1, 1)
    semilogx(f, mag2db(abs(H_true)), '--',...
             'color', cmap(end, :), 'linewidth', 2);
    hold on;
    semilogx(f, mag2db(abs(minus_Ginv_true)), '--',...
             'color', cmap(1, :), 'linewidth', 2);
    for i = 1:length(H_i)
        semilogx(f_out, mag2db(abs(H_i{i})),...
                 'color', cmap(end - i, :), 'linewidth', 2);
    end
    hold off;
    ylabel('magnitude(H(f)) [dB]', 'Fontsize', 18);
    xlim([0.01 fmax])

subplot(3, 1, 2)
    semilogx(f, 180/pi*phase(H_true), '--',...
             'color', cmap(end, :), 'linewidth', 2);
    hold on;
    semilogx(f, 180/pi*phase(minus_Ginv_true) - 180, '--',...
             'color', cmap(1, :), 'linewidth', 2);
    for i = 1:length(H_i)
        phase_offset = 180;
        if i == 1
            phase_offset = 0;
        end
        semilogx(f_out, 180/pi*phase(H_i{i}) - phase_offset,...
                 'color', cmap(end - i, :), 'linewidth', 2);
    end
    hold off;
    l = legend(lgd, 'location', 'west');
    set(l, 'FontSize', 18);
    ylabel('phase(H(f)) [deg]', 'Fontsize', 18);
    ylim([-180  180]);
    xlim([0.01 fmax]);

subplot(3, 1, 3)
    for i = 1:length(H_i)
        semilogx(f_out, abs(C_i{i}),...
                 'color', cmap(end - i, :), 'linewidth', 2);
        hold on;
    end
    hold off;
    xlabel('frequency [Hz]', 'Fontsize', 18);
    ylabel('magnitude(C_{yu}(f))', 'Fontsize', 18);
    ylim([0 1.1]);
    xlim([0.01 fmax]);
eps_save('question3b')

%% 3d
rng(1465) % set seed for generation of white noise input
n = normrnd(0, 20, size(t));
output = lsim(sys, [r, n], t);
u = output(:, 1);
y = output(:, 2);

[H_yu, C_yu, f_out] = estimatesystem(y, u, 30, f);

Srr = welchspectrum(r, r, nsegment);
Sur = welchspectrum(u, r, nsegment);
Syr = welchspectrum(y, r, nsegment);
Syy = welchspectrum(y, y, nsegment);

H_cl = Syr ./ Sur;
C_cl = abs(Syr).^2 ./ (Srr .* Syy);

H_i = {H_yu, H_cl};
C_i = {C_yu, C_cl};
f_i = {f_out, f_out};

plotspectrum(f, H_true, f_i, H_i, C_i, flipud(hsv(1 + length(f_i))),...
    {'True', 'H_{ol}', 'H_{cl}'});
eps_save('question3d')
