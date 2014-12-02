function [out] = freqAvg(in,nrbands)
%
% function [out] = freqAvg(in,nrbands)
%
% Averages the input 'in' over the number of frequency bands 'nrbands'
% version: February 15, 2007

N = length(in);

[m,n]   = size(in);

% 'in' should be a vector (not an array)
if m > n
    in=in.';
end


Nmod    = floor((N/nrbands));   % number of remaining frequencies after averaging
tmp     = zeros(nrbands,Nmod);	% initialization of temporary matrix for averaging: nrband rows and nmod columns
tmp(:)  = in(1:nrbands*Nmod);   % arrange the samples of 'in' in the elements of tmp.
out     = mean(tmp,1);          % average over columns and make it a vector

% If the input was a column vector, out needs to be transposed
if m > n
    out = out.';
end

