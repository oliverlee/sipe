clear all
close all
clc
load('TVdatset')
dt=t(2)-t(1);
Ns=length(t);
axfs = 18;
lgndfs = 16;

%% A. Data plots
figure(1)
subplot(511)
plot(t,uk)
ylabel('uk')
subplot(512)
plot(t,yk)
ylabel('yk')
subplot(513)
plot(t,ek(:,1))
ylabel('ek1')
subplot(514)
plot(t,ek(:,2))
ylabel('ek2')
subplot(515)
plot(t,ek(:,3))
ylabel('ek3')

u0 = uk
ek1 = ek(:, 1);
ek2 = ek(:, 2);
ek3 = ek(:, 3);

%% Get LTI state space model
n=2;s=10;
u1 = u0
u2 = [u0 ek]
[As,Bs,Cs,Ds,Ks,S]=getsshp(u1,yk,n,s);
sysi1 = ss(As,Bs,Cs,Ds,dt);
yi1 = lsim(sysi1,u1,t);
[As,Bs,Cs,Ds,Ks,S]=getsshp(u2,yk,n,s);
sysi2 = ss(As,Bs,Cs,Ds,dt);
yi2 = lsim(sysi2,u2,t);

figure(3)
S0 = S
semilogy(S,'b*')
xlabel('singular value index', 'fontsize', axfs);
ylabel('singular value', 'fontsize', axfs);
eps_save('svd3i', figure(3))

figure(2)
colset = cool(4);
subplot(311)
plot(t, u0, 'color', colset(1, :), 'linewidth', 2); hold on;
plot(t, ek1, 'color', colset(2, :), 'linewidth', 2);
plot(t, ek2, 'color', colset(3, :), 'linewidth', 2);
plot(t, ek3, 'color', colset(4, :), 'linewidth', 2); hold off;
ylabel('uk [N-m]', 'fontsize', axfs)
l1 = legend({'uk', 'ek1', 'ek2', 'ek3'})
set(l1, 'fontsize', lgndfs);
subplot(312)
plot(t, yk, 'color', colset(1, :), 'linewidth', 2); hold on;
plot(t, yi1, '--', 'color', colset(2, :), 'linewidth', 2);
plot(t, yi2, '--', 'color', colset(4, :), 'linewidth', 2); hold off;
ylabel('yk [degrees]', 'fontsize', axfs)
l1 = legend({'yk', 'yi1', 'yi2'})
set(l1, 'fontsize', lgndfs);
subplot(313)
Nvaf=25;
plot(t,mvaf(yk,yi1,Nvaf), '--', 'color', colset(2, :), 'linewidth', 2), hold on;
plot(t,mvaf(yk,yi2,Nvaf), '--', 'color', colset(4, :), 'linewidth', 2), hold off;
ylabel('VAF', 'fontsize', axfs)
xlabel('time [s]', 'fontsize', axfs)
l1 = legend({'yi1', 'yi2'})
set(l1, 'fontsize', lgndfs);
eps_save('vaf3j', figure(2))

uk = u0;
[As,Bs,Cs,Ds,Ks,S]=getsshp(u2,yk,n,s);
sysi = ss(As,Bs,Cs,Ds,dt);
yi = lsim(sysi,u2,t);

%% B.
% LPV Identification
mucutoff=.25;                                % change filter cutoff here
[bmu,amu]=butter(2,mucutoff/.5*dt);
colset = cool(4);
f1 = figure()
for i = 1:3
muk=[ones(1,Ns);filtfilt(bmu,amu,ek(:,i))']'; % create scheduling function here
[Al,Bl,Cl,Dl,Kl,S]=getlpv(uk,yk,muk,n);
yid=simlpv(Al,Bl,Cl,Dl,Kl,uk,muk(:,2:end));

%
%figure()
%subplot(411)
%plot(uk)
%subplot(412)
%plot(yk), hold on
%plot(yi,'m--')
%plot(yid,'r--')
%
%subplot(413)
%plot(t,muk)
%
%subplot(414)
Nvaf=25;
plot(t,mvaf(yk,yid,Nvaf),'-', 'color', colset(i, :), 'linewidth', 2), hold on;
end
plot(t,mvaf(yk,yi,Nvaf),'-', 'color', colset(4, :), 'linewidth', 2), hold off;
ylabel('VAF', 'fontsize', axfs)
xlabel('time [s]', 'fontsize', axfs)
l1 = legend('ylpv (ek1)', 'ylpv (ek2)', 'ylpv (ek3)', 'ylti')
set(l1, 'fontsize', lgndfs);
ylim([0, 101])

%%
muk=[ones(1,Ns);filtfilt(bmu,amu,ek(:,:))']'; % create scheduling function here
[Al,Bl,Cl,Dl,Kl,S]=getlpv(uk,yk,muk,n);
Al
mesh(abs(Al))
xlabel('x', 'fontsize', axfs)
ylabel('y', 'fontsize', axfs)
zlabel('z', 'fontsize', axfs)
eps_save('mesh3m', gcf);

%% C. Bode plots
close all;
muk=[ones(1,Ns);filtfilt(bmu,amu,ek(:,2))']'; % create scheduling function here
[Al,Bl,Cl,Dl,Kl,S]=getlpv(uk,yk,muk,n);

NSteps=5;
muSteps=min(muk(:,2)):(max(muk(:,2))-min(muk(:,2)))/(NSteps-1):max(muk(:,2));
fv=0.1:0.01:5;
H=zeros(length(fv),NSteps);

for jj=1:length(muSteps)
    % create linear ss
    Ac=Al(:,1:n)+muSteps(jj)*Al(:,n+1:2*n); % Ac=combined matrix Af=found LPV 
    Bc=(Bl(:,1)+muSteps(jj)*Bl(:,2));
    Cc=Cl(1:n);Dc=Dl(1);
    sysLPV=(ss(Ac,Bc,Cc,Dc,dt));
    H(:,jj)=squeeze(freqresp(sysLPV,2*pi*fv));
end

colset = cool(5);
lgnd = cell(size(colset, 1), 1);
subplot(211)
for i = 1:size(colset, 1)
    loglog(fv, abs(H(:, i)), 'color', colset(i, :)); hold on;
    lgnd{i} = sprintf('\\mu_k = %g', muSteps(i));
end
hold off;
ylabel('gain', 'fontsize', axfs)
h2 = subplot(212)
semilogx(fv, unwrap(180/pi*(angle(H))));
for i = 1:size(colset, 1)
    semilogx(fv, 180/pi*unwrap(angle(H(:, i))), 'color', colset(i, :)); hold on;
end
hold off;
ylabel('phase [degrees]', 'fontsize', axfs)
xlabel('frequency [Hz]', 'fontsize', axfs)
l1 = legend(lgnd);
set(l1, 'fontsize', lgndfs);
eps_save('frf3n')
