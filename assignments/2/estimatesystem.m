function [H, C, fout] = estimatesystem(y, u, nsegment, f)
if nargin < 4
    f = 1;
end

Suu = welchspectrum(u, u, nsegment);
Syu = welchspectrum(y, u, nsegment);
Syy = welchspectrum(y, y, nsegment);
H = Syu ./ Suu;
C = abs(Syu).^2 ./ (Suu .* Syy);
fout = f(1:nsegment:floor(end/nsegment)*nsegment);

end