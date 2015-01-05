function [A, B, C, D, x0] = estimatesystem(YsN, UsN, s, n)

[U, S, V] = svd(YsN, 'econ');           % Singular value decomposition
Un = U(:, 1:n);                        % Reduced output singular vectors
Vn = V(:, 1:n);                        % Reduced input singular vectors
Sn = S(1:n, 1:n);                      % Reduced singular value matrix
SnVtn = Sn*Vn';

% n is size of state
% r is size of input
% l is size of output
r = size(UsN, 1)/s;
l = size(YsN, 1)/s;

A = Un(1:(s - 1)*l, :)\Un(l + 1:s*l, :);
B = SnVtn(:, 1:r);
C = Un(1:l, :);
D = YsN(1:l, 1:r);
x0 = SnVtn(:, 1);
