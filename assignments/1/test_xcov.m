%% test out crosscov and xcov

n = randn(100000, 1);
disp(sprintf('size of input is %d', length(n)));

tic
c1 = xcov(n, n, 'biased');
t1 = toc;
disp(sprintf('time for xcov(n, n) is %g', t1))

tic
c2 = crosscov(n, n, 'biased');
t2 = toc;
disp(sprintf('time for crosscov(n, n) is %g', t2))
