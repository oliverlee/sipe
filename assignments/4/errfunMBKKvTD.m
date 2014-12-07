function e = errfunMBKKvTD(P, ut, yt, t);

PP = num2cell(P);
[M, B, K, Kv, td, w] = deal(PP{:});

set_param('modmbkkv/invM', 'gain', num2str(1/M));
set_param('modmbkkv/B', 'gain', num2str(B));
set_param('modmbkkv/K', 'gain', num2str(K));
set_param('modmbkkv/Kv', 'gain', num2str(Kv));
set_param('modmbkkv/delay', 'delay', num2str(td));
set_param('modmbkkv/Hact',...
    'numerator', ['[ ' num2str(w.^2) ' ]']);
set_param('modmbkkv/Hact',...
    'denominator', ['[ 1 ' num2str(2*0.7*w) ' ' num2str(w.^2) ']']);
[~, ~, ysim] = sim('modmbkkv', t, [], [t ut]);

e = yt - ysim;
