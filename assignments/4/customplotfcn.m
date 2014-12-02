function state = customplotfcn(options, state, flag)
% Custom function to plot intermediate populations during optimization.
% Author: Jasper Schuurmans and Winfred Mugge

persistent hgen

gen = state.Population;
if(strcmp(flag,'init')) % Set up the plot
    figure(99)
    hold on;   
    hgen = plot(gen(:,1), gen(:,2), 'r.', 'markersize', 15);
    set(gca, 'fontsize', 15)
    set(gcf, 'units', 'normalized', 'position', [0.1 0.1 0.8 0.8])
    title('This is the initial population. Press a key to start...')
    xlabel('a')
    ylabel('b')
    pause
else
    figure(99)
    set(hgen, 'xdata', gen(:,1), 'ydata', gen(:,2))
end
title(['Population after iteration ', num2str(state.Generation)])
drawnow
% pause(0.1) % Remove comment to slow down!