close all
clear all
clc
%% parameters & system
Nt=1000;
dt=0.1;
t=(0:Nt-1).'*dt;
n = 6;        % system order
l = 3;        % no. of outputs
r = 3;        % no. of inputs
 
%% Generate random system
% Create your mimo system here!
sys=[];
%% Simulation
ucl=randn(Nt,r);            % clean input
uk=ucl+1e-1*randn(Nt,r);    % add noise
ycl=lsim(sys,ucl,t);        % clearn output
yk=ycl+1*randn(Nt,l);       % add noise

%% Get state space model
s=5;       % window size of data matrices
[Ai,Bi,Ci,Di,Ki,S]=getsshp(uk,yk,n,s);
figure(4)
semilogy(S,'r*')
sysi=ss(Ai,Bi,Ci,Di,dt);
yi=lsim(sysi,ucl,t);

%% plots of in- and outputs
figure(1)
for ii=1:r
    subplot(r,1,ii)
    plot(t,uk(:,ii),t,ucl(:,ii));
end
ylabel('uk')

figure(2)
for ii=1:l
    subplot(l,1,ii)
    plot(t,yk(:,ii),t,ycl(:,ii),t,yi(:,ii));
end
ylabel('yk')
legend('y0','yk','yi')


disp(['VAF with noisefree in/output = ' num2str(reshape(vaf(yi,ycl),1,l))])

%% Bode plots
figure(3);
bode(sys,sysi)
