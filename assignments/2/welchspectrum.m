function wSyu = welchspectrum(y,u,Nsegment)
% segments the data in Nsegment and calculates the spectrum averaged over
% the segments
% version: February 20, 2008

N     = length(u);
[m,n] = size(u);

% u and y should be vectors (not an array)
if m > n
    u = u.';
    y = y.';
end

% Number of samples per segment
Nd  = floor(N/Nsegment);

% data should an integer number of segments 
u   = u(1:Nd*Nsegment);
y   = y(1:Nd*Nsegment);

%reshape data, fft, and take the average
u   = reshape(u,Nd,Nsegment);
y   = reshape(y,Nd,Nsegment);
U   = fft(u);
Y   = fft(y);
Syu = 1/Nd * Y .* conj(U) ;
wSyu = mean(Syu,2).';

% If the input was a column vector, out needs to be transposed
if m > n
    wSyu = wSyu.';
end