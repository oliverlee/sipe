function [e, Hef] = errfunMBK(P, frf, fvec, mCoh, est)
if nargin < 5
    est = false;
end

PP = num2cell(P);
[M, B, K] = deal(PP{:});

s = 2*pi*1j*fvec;
s2 = -4*pi^2*fvec.^2;

Hef = 1./(M*s2 + B.*s + K);

if ~est
    e = abs(log(frf./Hef));
else
    e = abs(log(frf./Hef)) .* mCoh .* sqrt(1./fvec);
end
