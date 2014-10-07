function plotspectrum(f_true, H_true, f_i, H_i, C_i, cmap, legendstr)

fmax = f_true(round(end/2));

%cmap = jet(length(f_i) + 1);
%legendstr = cell(1, length(f_i) + 1);
if ~iscell(f_i)
    f_i = {f_i};
end
if ~iscell(H_i)
    H_i = {H_i};
end
if ~iscell(C_i)
    C_i = {C_i};
end
N = length(f_i);

figure()

subplot(3, 1, 1)
    semilogx(f_true, mag2db(abs(H_true)), '--',...
             'color', cmap(end, :), 'linewidth', 2);
    hold on;
    for i = 1:N
        semilogx(f_i{i}, mag2db(abs(H_i{i})),...
                 'color', cmap(end - i, :), 'linewidth', 2);
    end
    hold off;
    ylabel('magnitude(H(f)) [dB]', 'Fontsize', 18);
    xlim([0.01 fmax])

subplot(3, 1, 2)
    semilogx(f_true, 180/pi*phase(H_true), '--',...
             'color', cmap(end, :), 'linewidth', 2);
    hold on;
    for i = 1:N
        semilogx(f_i{i}, 180/pi*phase(H_i{i}),...
                 'color', cmap(end - i, :), 'linewidth', 2);
    end
    hold off;
    l = legend(legendstr, 'location', 'east');
    set(l, 'FontSize', 18);
    ylabel('phase(H(f)) [deg]', 'Fontsize', 18);
    ylim([-180  180]);
    xlim([0.01 fmax]);

subplot(3, 1, 3)
    for i = 1:N
        semilogx(f_i{i}, abs(C_i{i}),...
                 'color', cmap(end - i, :), 'linewidth', 2);
        hold on;
    end
    hold off;
    xlabel('frequency [Hz]', 'Fontsize', 18);
    ylabel('magnitude(C_{yu}(f))', 'Fontsize', 18);
    ylim([0 1.8]);
    xlim([0.01 fmax]);
