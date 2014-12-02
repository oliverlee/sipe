%% Assignment 3: Identification
%  [WB2301] SIPE 2014-2015

%% Question 1: Time domain models in closed loop
clear all
close all
clc

% Load model
options=simset('OutputPoints','specified'); % only produce output for the time values in vector <t>
load_system('Ass3_SISOclosed')

r_in=!!!

%% Question 3: Detecting nonlinearities
clear all
close all
clc

% Load model
options=simset('OutputPoints','specified'); % only produce output for the time values in vector <t>
load_system('Ass3_SISOopen_NL');

%Default values
set_param('Ass3_SISOopen_NL/Noise n','Variance','0.1');   
set_param('Ass3_SISOopen_NL/Dead Zone','LowerValue','-0.1');
set_param('Ass3_SISOopen_NL/Dead Zone','UpperValue','0.1');


%% Question 3e: Fast method

% Some administration
T   = 200;                  % Total time
fs  = 200;                  % Sample frequency
N   = T*fs;                 % number of samples
tv   = (0:N-1).'/fs;        % Time vector
fv   = (0:N-1).'/T;         % Frequency vector (double sided)

%load input signal
load('rSignals.mat')

u_in=[tv repmat(rFastMethod,10,1)];

sim('Ass3_SISOopen_NL',tv,options)

%% DO NOT CHANGE:
M=1;P=10;
ur=reshape(u,N/P,P);
yr=reshape(y,N/P,P);
N=length(ur(:,1));
Ur=fft(ur)/sqrt(N);
Yr=fft(yr)/sqrt(N);
MeasHarm=1:N/2;
Uall=zeros(M,P,length(MeasHarm));
Yall=zeros(M,P,length(MeasHarm));
Uall(1,:,:) = Ur(2:max(MeasHarm)+1,:).';
Yall(1,:,:) = Yr(2:max(MeasHarm)+1,:).';

%Fast method
[Y, Yc, U, H, freq] = Fast_NL_Anal(Yall, Uall, rFastMethodExHarm, MeasHarm, fs, N);

% True system transfer function of the linear block
s = 2*pi*1j*freq.E';
Hlin = 1./(0.01*s.^2+0.03*s+1);

% plot
figure
loglog(freq.E,abs(Hlin),'k--',freq.E, abs(H.mean), 'k', freq.E, abs(H.stdn.E), 'g', freq.E, abs(H.stdNL), 'r')
title('FAST - ONE REAL.')
xlabel('Frequency [Hz]')
xlim([0.05 40])
legend('H_{0}','H_{BLA}','\sigma_{noise}','\sigma_{total}')
