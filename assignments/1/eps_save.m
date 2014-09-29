function eps_save(filename, figh)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    figh = gcf;
end
position = get(gcf,'Position');
set(gcf, 'PaperPositionMode', 'auto',...
    'Position', [position(1:2) 1400 1000])
hgexport(figh, [filename, '.eps'],...
    hgexport('factorystyle'), 'Format', 'epsc2');
end

