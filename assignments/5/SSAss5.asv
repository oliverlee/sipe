clear all
close all
clc
load('TVdatset')
dt=t(2)-t(1);
Ns=length(t);
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

%% Get LTI state space model
n=2;s=10;
uk=[uk ek];
[As,Bs,Cs,Ds,Ks,S]=getsshp(uk,yk,n,s);
sysi=ss(As,Bs,Cs,Ds,dt);
yi=lsim(sysi,uk,t);

figure(3)
semilogy(S,'b*')

figure(2)
subplot(311)
plot(t,uk)
ylabel('uk')
subplot(312)
plot(t,yk,t,yi)
ylabel('yk')
subplot(313)
Nvaf=25;
plot(t,mvaf(yk,yi,Nvaf),'m--'), hold on;


% B.
% LPV Identification
mucutoff=.25;                                % change filter cutoff here
[bmu,amu]=butter(2,mucutoff/.5*dt);
muk=[ones(1,Ns);filtfilt(bmu,amu,ek(:,2))']'; % create scheduling function here
[Al,Bl,Cl,Dl,Kl,S]=getlpv(uk,yk,muk,n);
yid=simlpv(Al,Bl,Cl,Dl,Kl,uk,muk(:,2:end));

%
figure(2)
subplot(411)
plot(uk)
subplot(412)
plot(yk), hold on
plot(yi,'m--')
plot(yid,'r--')

subplot(413)
plot(t,muk)

subplot(414)
Nvaf=25;
plot(t,mvaf(yk,yid,Nvaf),'r--'), hold on;
plot(t,mvaf(yk,yi,Nvaf),'m--'), hold on;
legend('ylpv','ylti')


% C. Bode plots
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

figure
subplot(211)
loglog(fv,abs(H));
subplot(212)
semilogx(unwrap(angle(H)));
