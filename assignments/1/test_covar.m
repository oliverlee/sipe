%% test out covar and xcov

n = randn(100000, 1);

tic
c1 = xcov(n, n, length(n) - 1, 'biased');
t1 = toc

tic
c2 = covar(n, n, length(n) - 1, 0);
t2 = toc

t2 / t1

size(c1)
size(c2)