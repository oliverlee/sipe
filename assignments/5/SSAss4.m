close all
clear all
clc

axfs = 18;
lgndfs = 16;

%% parameters & system
Nt=1000;
dt=0.1;
t=(0:Nt-1).'*dt;
n = 6;        % system order
l = 3;        % no. of outputs
r = 3;        % no. of inputs
 
%% Generate random system
% Create your mimo system here!

m1 = 0.1; m2 = 0.1; m3 = 0.1;
c1 = 0.5; c2 = 0.5; c3 = 0.5;
k1 = 450; k2 = 450; k3 = 450;
A0 = [
    -(c1 + c2)/m1,         c2/m1,      0, -(k1 + k2)/m1,         k2/m1,      0;
            c1/m2, -(c2 + c3)/m2,  c3/m2,         k1/m2, -(k2 + k3)/m2,  k3/m2;
                0,         c3/m3, -c3/m3,             0,         k3/m3, -k3/m3;
        eye(3), zeros(3)];
B0 = [diag([1/m1, 1/m2, 1/m3]); zeros(3)];
C0 = [zeros(3), eye(3)];
D0 = 0;
sys = c2d(ss(A0, B0, C0, D0), dt);         % discrete system representation

%% Simulation
ucl=randn(Nt,r);            % clean input
uk=ucl+1e-1*randn(Nt,r);    % add noise
ycl=lsim(sys,ucl,t);        % clearn output
yk=ycl+1*randn(Nt,l);       % add noise

%% Get state space model
s=20;       % window size of data matrices TODO
[Ai,Bi,Ci,Di,Ki,S]=getsshp(uk,yk,n,s);
figure(4)
semilogy(S,'r*')
set(gcf, 'name', 'singular values');
xlabel('singular value index', 'fontsize', axfs);
ylabel('singular value', 'fontsize', axfs);
%disp('singular values')
%S
%disp('singular values differences')
%abs(diff(log(S)))
sysi=ss(Ai,Bi,Ci,Di,dt);
yi=lsim(sysi,ucl,t);

%% plots of in- and outputs
tlim = [0, 4];

figure(1)
for ii=1:r
    subplot(r,1,ii)
    plot(t,ucl(:,ii),t,uk(:,ii));
    ylabel(sprintf('u%d', ii), 'fontsize', axfs)
    xlim(tlim);
end
xlabel('time [s]', 'fontsize', axfs)
l1 = legend('ucl','uk');
set(l1, 'fontsize', lgndfs);
subplot(r,1,round(r/2))
eps_save('input2h', figure(1));

figure(2)
for ii=1:l
    subplot(l,1,ii)
    plot(t,ycl(:,ii),t,yk(:,ii),t,yi(:,ii));
    ylabel(sprintf('y%d', ii), 'fontsize', axfs)
    xlim(tlim);
end
xlabel('time [s]', 'fontsize', axfs)
l2 = legend('ycl','yk','yi');
set(l2, 'fontsize', lgndfs);
subplot(r,1,round(r/2))
eps_save('output2h', figure(2));

disp(['VAF with noisefree in/output = ' num2str(reshape(vaf(yi,ycl),1,l))])

%% Bode plots
figure(3);
bode(sys,sysi)
l3 = legend('orig','iden');
set(l3, 'fontsize', lgndfs);
eps_save('bode2h', figure(3));
