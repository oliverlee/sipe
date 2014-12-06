function [e, Hef] = errfunMBKKv(P, frf, fvec, mCoh)

PP = num2cell(P);
[M, B, K, Kv, td, w] = deal(PP{:});

s = 2*pi*1j*fvec;
s2 = -4*pi^2*fvec.^2;

Hact = w^2./(s2 + 2*0.7*w*s + w^2);
Hef = 1./(M*s2 + (B + Kv*exp(-td*s).*Hact).*s + K);

e = abs(log(frf./Hef)) .* mCoh .* sqrt(1./fvec);
