function [c, lags] =  crosscov(y, u, maxlag, scaleopt)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    scaleopt = 'biased';
end
if nargin < 3
    maxlag = length(y) - 1;
end
if ischar(maxlag)
    scaleopt = maxlag;
    maxlag = length(y) - 1;
end

lags = -maxlag:maxlag;
n = 2*length(y) - 1;

% xcov(y, x) = xcorr(y -mean(y), u - mean(u))
% where xcorr can be efficiently computed in the frequency domain as a
% multiplcation instead of a convolution in the time domain.
phi = ifft(conj(fft(u - mean(u), n)) .* fft(y - mean(y), n));
if size(phi, 1) == 1
    phi = [phi(end - (maxlag + 1) + 2:end), phi(1:maxlag + 1)];
else
    phi = [phi(end - (maxlag + 1) + 2:end); phi(1:maxlag + 1)];
end

bias = length(y) * ones(size(phi));
if strcmp(scaleopt, 'unbiased')
    bias = bias - abs(lags);
end
c = phi ./ bias;

end