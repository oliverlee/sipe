function [c, lags] =  covar(y, u, maxlag, unbiased)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
n = length(y);
if nargin < 4
    unbiased = 0;
end
if nargin < 3
    maxlag = n - 1;
end

% if size(y, 1) == 1
%     y = y' - mean(y);
% else
%     y = y - mean(y);
% end
% 
% if size(u, 1) == 1
%     u = u - mean(u);
% else
%     u = u' - mean(u);
% end
y = y - mean(y);
u = u - mean(u);

lags = (-maxlag:maxlag)';
m = length(lags);
c = zeros(size(lags));

% bias = n*ones(size(lags));
% if ~biased
%     bias = bias - abs(lags);
% end
if unbiased ~= 1
    unbiased = 0;
end

% if isempty(gcp)
%     parpool;
% end

% parfor i = 1:maxlag + 1
for i = 1:maxlag + 1
    k = maxlag - i + 1;
    bias = n - k * unbiased;
    c(i) = dot(u(1 + k:end), y(1:end - k)) / bias;
end
% parfor i = maxlag + 2:m
for i = maxlag + 2:m
    k = i - maxlag - 1;
    bias = n - k * unbiased;
    c(i) = dot(u(1:end - k), y(1 + k:end)) / bias;
end


% c = c./bias;

% 
% y_conv = zeros(n, m);
% for i = 1:maxlag + 1
%     k = maxlag - i + 1;
%     y_conv(1:end - k, i) = y(1 + k:end);
% end
% for i = maxlag + 1:m - 1
%     k = i - maxlag;
%     y_conv(1 + k:end, i + 1) = y(1:end - k);
% end
% 
% bias = n*ones(size(lags));
% if ~biased
%     bias = bias - abs(lags);
% end
% 
% c = (u * y_conv)' ./ bias;


% % the columns of this matrix are circular shifted versions of y
% tic
% y_shift = cell2mat(...
%     arrayfun(@(i) circshift(y, i), lags, 'uniformoutput', false));
% toc
% % create map with zero values
% tic
% y_zeros = [triu(ones(n, maxlag), maxlag - n + 1), ...
%            tril(ones(n, maxlag + 1), 0)];
% toc
% tic
% y_conv = y_zeros.*y_shift;
% toc
% bias = n*ones(size(lags));
% if ~biased
%     bias = bias - [-1*(-maxlag:0), 1:maxlag];
% end
% 
% tic
% c = (u * y_conv) ./ bias;
% toc
% 
end