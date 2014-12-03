function [minpar, mine, iter, time, h] = combinedsearch(ares, bres, plotfig)
if nargin < 3
    plotfig = false;
end

if plotfig
    global pathway
    pws = cell([length(ares) * length(bres), 3]);
else
    h = 0;
end
tic;

lb = [min(ares), min(bres)];
ub = [max(ares), max(bres)];
options = optimset('lsqnonlin');
options = optimset(options, 'Display', 'off');
options = optimset(options, 'outputfcn', @storepathway);
iter = 0;
mine = inf;

for qq = 1:length(ares)
    for ss = 1:length(bres)
        if plotfig
            [pathway.x pathway.y, pathway.z] = deal([]);
        end
        par = [ares(qq) bres(ss)];
        [parest, ~, e, ~, out, ~, ~] =...
            lsqnonlin(@errfun2, par, lb, ub, options);
        sse = e' * e;
        if sse < mine
            mine = sse;
            minpar = parest;
        end
        iter = iter + out.iterations;
        if plotfig
            pws((qq - 1)*length(bres) + ss, :) =...
                {pathway.x, pathway.y, pathway.z};
        end
    end
end
time = toc;

if plotfig
    h = figure(); hold on;
    for i = 1:size(pws, 1)
        px = pws{i, 1}; py = pws{i, 2}; pz = pws{i, 3};
        surface([px(:), px(:)], [py(:), py(:)], [pz(:), pz(:)],...
            'edgecolor', 'flat', 'facecolor', 'none', 'linewidth', 4);
        plot3(px(1), py(1), pz(1), 'k.', 'markersize', 15);
        plot3(px(end), py(end), pz(end), 'k.', 'markersize', 25);
    end
    xlabel('a', 'fontsize', 16);
    ylabel('b', 'fontsize', 16);
    zlabel('e', 'fontsize', 16);
end
