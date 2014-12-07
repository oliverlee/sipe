%% SIPE (2014-2015): Assignment 4
%% PART2

clear all
close all
clc

%% PART A: Load DATA & Plot FRF's
% Load DATA
load('DatRec1.mat')
t0=matrix(1,:)';        % extract vectors from data-matrix
u0=matrix(3,:)';
y0=matrix(4,:)';
fs=500;
dt=1/fs;
T0=100;
N=T0*fs+1;
fv0=(0:1/T0:fs)';
colset=cool(4);
nBands=4;

figure(1)
set(gcf,'Name','Measured Signals')
subplot(211)
plot(t0,y0)
ylabel('angle [rad]', 'fontsize', 18);
subplot(212)
plot(t0,u0)
xlabel('time [s]', 'fontsize', 18);
ylabel('torque [N-m]', 'fontsize', 18);
% FRF's
% FRF whole dataset
U0=fft(detrend(u0));
Y0=fft(detrend(y0));

Suu0=U0.*conj(U0);
Syy0=Y0.*conj(Y0);
Syu0=Y0.*conj(U0);

mSuu0=freqAvg(Suu0,nBands);
mSyu0=freqAvg(Syu0,nBands);
mSyy0=freqAvg(Syy0,nBands);
mfv0=freqAvg(fv0,nBands);

H0=Syu0./Suu0;
mH0=mSyu0./mSuu0;
mCoh0=abs(mSyu0).^2./(mSuu0.*mSyy0);

% Plot FRF
figure(2)
set(gcf,'Name','FRFs')
subplot(321)
    loglog(fv0,abs(H0),'linewidth',1,'linestyle','--'); hold on;
    loglog(mfv0,abs(mH0),'r','linewidth',2);
    legend('H0','mH0');
    ylabel('Gain', 'fontsize', 18); xlim([1/T0 fs/2]);
subplot(323)
    semilogx(fv0,unwrap(angle(H0))*180/pi,...
        'linewidth',1,'linestyle','--'); hold on;
    semilogx(mfv0,unwrap(angle(mH0))*180/pi,'r','linewidth',2);
    ylabel('Phase [deg]', 'fontsize', 18);
    xlim([1/T0 fs/2]); ylim([-360 360]);
subplot(325)
    semilogx(mfv0,mCoh0,'linewidth',2); hold on;
    xlabel('Frequency [Hz]', 'fontsize', 18);
    ylabel('Coherence', 'fontsize', 18);
    xlim([1/T0 fs/2]); ylim([0 1.1]);

% FRF per segment
segments=[2 22 ;27 47;52 72;77 97];

for ii=1:4
    tidx=segments(ii,1)*fs:(segments(ii,2)*fs);
    t=matrix(1,tidx)';
    u=u0(tidx);
    y=y0(tidx);
    T=segments(ii,2)-segments(ii,1);
    fv=(0:1/T:fs)';

    U=fft(detrend(u));
    Y=fft(detrend(y));

    Suu=U.*conj(U);
    Syy=Y.*conj(Y);
    Syu=Y.*conj(U);

    mSuu=freqAvg(Suu,nBands);
    mSyu=freqAvg(Syu,nBands);
    mSyy=freqAvg(Syy,nBands);
    mfv=freqAvg(fv,nBands);

    H=Syu./Suu;
    mH=mSyu./mSuu;
    mCoh=abs(mSyu).^2./(mSuu.*mSyy);

    % Plot FRF
    figure(2)
    set(gcf,'Name','FRFs')
    subplot(322)
        loglog(fv,abs(H),'color',colset(ii,:),'linewidth',1,'linestyle','--');
        hold on;
        loglog(mfv,abs(mH),'color',colset(ii,:),'linewidth',2);
        ylabel('Gain', 'fontsize', 18); xlim([1/T0 fs/2]);
    subplot(324)
        semilogx(fv,unwrap(angle(H))*180/pi,'color',colset(ii,:),...
            'linewidth',1,'linestyle','--');
        hold on;
        semilogx(mfv,unwrap(angle(mH))*180/pi,'color',colset(ii,:),...
            'linewidth',2);
        ylabel('Phase [deg]', 'fontsize', 18); xlim([1/T0 fs/2]);
    subplot(326)
        semilogx(mfv,mCoh,'linewidth',2,'color',colset(ii,:)); hold on;
        xlabel('Frequency [Hz]', 'fontsize', 18);
        ylabel('Coherence', 'fontsize', 18);
        xlim([1/T0 fs/2]); ylim([0 1.1]);
    figure(1)
    subplot(212)
    for jj=1:2
        line([segments(ii,jj) segments(ii,jj)],[3 -3],'color',colset(ii,:))
    end
end
figure(1); eps_save('wristdata')
figure(2); eps_save('wristfrf')

%% PART B : Estimate MBK model
% for each segment
figure(3); clf
colset = cool(12);
err = zeros(4, 2);

for ii=1:4
    tidx=segments(ii,1)*fs:(segments(ii,2)*fs);
    t=matrix(1,tidx)';
    u=u0(tidx);
    y=y0(tidx);
    T=segments(ii,2)-segments(ii,1);
    fv=(0:1/T:fs)';

    U=fft(detrend(u));
    Y=fft(detrend(y));

    Suu=U.*conj(U);
    Syy=Y.*conj(Y);
    Syu=Y.*conj(U);

    mSuu=freqAvg(Suu,nBands);
    mSyu=freqAvg(Syu,nBands);
    mSyy=freqAvg(Syy,nBands);
    mfv=freqAvg(fv,nBands);

    H=Syu./Suu;
    mH=mSyu./mSuu;
    mCoh=abs(mSyu).^2./(mSuu.*mSyy);

    % LSQ-routine
    p0 = [0.002, 0.1, 5];
    lb0 = [0.001, 0.01, 1];
    ub0 = [0.005, 1, 20];

    %Estimate on these frequencies:
    estfreq=[1/T fs];
    estidx=find(mfv<estfreq(2)&mfv>estfreq(1));
    est2idx = find(mfv < 20  & mfv > 2);

    options = optimset('lsqnonlin');
    options = optimset(options, 'display', 'iter', 'TolX', 1e-15,...
        'diffminchange', 1e-4, 'maxfuneval', 500, 'TolFun', 1e-12);
    [P,~,R,~,~,~,J] = lsqnonlin(@errfunMBK, p0, lb0, ub0, options,...
        mH(estidx), mfv(estidx), mCoh(estidx));
    plsq(ii,:) = P;
    SEM(ii,:) = full(sqrt(R'*R*(diag(inv(J'*J)))/length(R)))';

    [E,Hid] = errfunMBK(plsq(ii,:), mH, mfv, mCoh);

    % compute improved parameter estimates
    [P2,~,R,~,~,~,J] = lsqnonlin(@errfunMBK, p0, lb0, ub0, options,...
        mH(est2idx), mfv(est2idx), mCoh(est2idx), true);
    [~, Hid2] = errfunMBK(P2, mH, mfv, mCoh);
    SEM2(ii,:) = full(sqrt(R'*R*(diag(inv(J'*J)))/length(R)))';

    e1 = errfunMBK(P, mH(est2idx), mfv(est2idx), mCoh(est2idx), true);
    e2 = errfunMBK(P2, mH(est2idx), mfv(est2idx), mCoh(est2idx), true);
    err(ii, :) = [e1'*e1, e2'*e2];

    % Plot FRF
    figure(3)
    set(gcf,'Name','MBK Estimates')
    subplot(311)
        loglog(mfv(est2idx),abs(mH(est2idx)),'x',...
            'color',colset((ii - 1)*3 + 1,:),...
            'linewidth', 4, 'linestyle', ':');
        hold on;
        loglog(mfv(est2idx),abs(Hid(est2idx)),...
            'color',colset((ii - 1)*3 + 2,:),...
            'linewidth',2,'linestyle','--');
        loglog(mfv(est2idx),abs(Hid2(est2idx)),...
            'color',colset((ii - 1)*3 + 3,:),...
            'linewidth',2,'linestyle','-');
        ylabel('Gain', 'fontsize', 18); xlim([1.9 21]);
    subplot(312)
        semilogx(mfv(est2idx),(angle(mH(est2idx)))*180/pi,'x',...
            'color',colset((ii - 1)*3 + 1,:),...
            'linewidth',4, 'linestyle', ':'), hold on;
        semilogx(mfv(est2idx),(angle(Hid(est2idx)))*180/pi,...
            'color',colset((ii - 1)*3 + 2,:),...
            'linewidth',2,'linestyle','--');
        semilogx(mfv(est2idx),(angle(Hid2(est2idx)))*180/pi,...
            'color',colset((ii - 1)*3 + 3,:),...
            'linewidth',2,'linestyle','-');
        ylabel('Phase [deg]', 'fontsize', 18); xlim([1.9 21]);
    subplot(313)
        semilogx(mfv(est2idx),mCoh(est2idx),'x','linewidth',4,...
            'linestyle', ':',...
            'color',colset((ii - 1)*3 + 1,:)); hold on;
        xlabel('Frequency [Hz]', 'fontsize', 18);
        ylabel('Coherence', 'fontsize', 18);
        ylim([0 1.1]); xlim([1.9 21]);
end
subplot(312)
l = legend({...
    'segment 1, measured', 'segment 1, estimate',...
    'segment 1, improved estimate',...
    'segment 2, measured', 'segment 2, estimate',...
    'segment 2, improved estimate',...
    'segment 3, measured', 'segment 3, estimate',...
    'segment 3, improved estimate',...
    'segment 4, measured', 'segment 4, estimate',...
    'segment 4, improved estimate'})
set(l, 'fontsize', 14);
figure(3); eps_save('mckfrf');

%% PART C:Initial values for Kv model
%p =[M B K Kv td w]
p02 = [0.003, 0.01, 5, -0.1, 0.04, 30];
lb02 = [0.001, 0.001, 1, -1e2, 0.02, 10];
ub02 = [0.005, 1, 100, -1e-2, 0.08, 50];

figure(5); clf
colset = cool(12);
err = zeros(4, 2);
plsq = zeros(4, 6)

for ii=1:4
    tidx=segments(ii,1)*fs:(segments(ii,2)*fs);
    t=matrix(1,tidx)';
    u=u0(tidx);
    y=y0(tidx);
    T=segments(ii,2)-segments(ii,1);
    fv=(0:1/T:fs)';

    U=fft(detrend(u));
    Y=fft(detrend(y));

    Suu=U.*conj(U);
    Syy=Y.*conj(Y);
    Syu=Y.*conj(U);

    mSuu=freqAvg(Suu,nBands);
    mSyu=freqAvg(Syu,nBands);
    mSyy=freqAvg(Syy,nBands);
    mfv=freqAvg(fv,nBands);

    H=Syu./Suu;
    mH=mSyu./mSuu;
    mCoh=abs(mSyu).^2./(mSuu.*mSyy);

    P = lsqnonlin(@errfunMBK, p0, lb0, ub0, options,...
        mH(est2idx), mfv(est2idx), mCoh(est2idx), true);

    % compute improved parameter estimates
    [P2,~,R,~,~,~,J] = lsqnonlin(@errfunMBKKv, p02, lb02, ub02, options,...
        mH(est2idx), mfv(est2idx), mCoh(est2idx));
    plsq(ii,:) = P2;
    SEM3(ii,:) = full(sqrt(R'*R*(diag(inv(J'*J)))/length(R)))';

    [e1, Hid] = errfunMBK(P, mH(est2idx), mfv(est2idx), mCoh(est2idx), true);
    [e2, Hid2] = errfunMBKKv(P2, mH(est2idx), mfv(est2idx), mCoh(est2idx));
    err(ii, :) = [e1'*e1, e2'*e2];

    % Plot FRF
    figure(5)
    set(gcf,'Name','MBKKv Estimates')
    subplot(211)
        loglog(mfv(est2idx),abs(mH(est2idx)),'x',...
            'color',colset((ii - 1)*3 + 1,:),...
            'linewidth', 4, 'linestyle', ':');
        hold on;
        loglog(mfv(est2idx),abs(Hid),...
            'color',colset((ii - 1)*3 + 2,:),...
            'linewidth',2,'linestyle','--');
        loglog(mfv(est2idx),abs(Hid2),...
            'color',colset((ii - 1)*3 + 3,:),...
            'linewidth',3,'linestyle','-');
        ylabel('Gain', 'fontsize', 18); xlim([1.9 21]);
    subplot(212)
        semilogx(mfv(est2idx),(angle(mH(est2idx)))*180/pi,'x',...
            'color',colset((ii - 1)*3 + 1,:),...
            'linewidth',4, 'linestyle', ':'), hold on;
        semilogx(mfv(est2idx),(angle(Hid))*180/pi,...
            'color',colset((ii - 1)*3 + 2,:),...
            'linewidth',3,'linestyle','--');
        semilogx(mfv(est2idx),(angle(Hid2))*180/pi,...
            'color',colset((ii - 1)*3 + 3,:),...
            'linewidth',2,'linestyle','-');
        ylabel('Phase [deg]', 'fontsize', 18); xlim([1.9 21]);
        xlabel('Frequency [Hz]', 'fontsize', 18);
end
subplot(211)
l = legend({...
    'segment 1, measured',...
    'segment 1, 3 param estimate',...
    'segment 1, 6 param estimate',...
    'segment 2, measured',...
    'segment 2, 3 param estimate',...
    'segment 2, 6 param estimate',...
    'segment 3, measured',...
    'segment 3, 3 param estimate',...
    'segment 3, 6 param estimate',...
    'segment 4, measured',...
    'segment 4, 3 param estimate',...
    'segment 4, 6 param estimate',...
})
set(l, 'fontsize', 14);
figure(5); eps_save('6parfrf');

%% PART D Determining the fit in the time domain.

figure(4); clf
colset = cool(4);
for ii=1:4
    tidx=segments(ii,1)*fs:(segments(ii,2)*fs);
    t=matrix(1,tidx)';
    ut=detrend(u0(tidx));
    yt=detrend(y0(tidx));
    M=plsq(ii,1);
    B=plsq(ii,2);
    K=plsq(ii,3);
    Kv=plsq(ii,4);
    td=plsq(ii,5);
    w=plsq(ii,6);
    set_param('modmbkkv/invM', 'gain', num2str(1/M));
    set_param('modmbkkv/B', 'gain', num2str(B));
    set_param('modmbkkv/K', 'gain', num2str(K));
    set_param('modmbkkv/Kv', 'gain', num2str(Kv));
    set_param('modmbkkv/delay', 'delay', num2str(td));
    set_param('modmbkkv/Hact', ...
        'numerator', ['[ ' num2str(w.^2) ' ]']);
    set_param('modmbkkv/Hact',...
        'denominator', ['[ 1 ' num2str(2*0.7*w) ' ' num2str(w.^2) ']']);

    [~,~,ysim] = sim('modmbkkv', t, [], [t ut]);
    figure(4)
    %subplot(211);
    %plot(t,ut)
    %subplot(212);
    %plot(t,yt,t,ysim)
    subplot(220 + ii);
    plot(t, yt, 'color', colset(2, :), 'linewidth', 2); hold on
    plot(t, ysim, 'color', colset(4, :), 'linewidth', 2, 'linestyle', '--');
    hold off
    ylabel('Angle [rad]', 'fontsize', 18);
    xlabel('Time [s]', 'fontsize', 18);
    l = legend({'measured', 'estimated'});
    set(l, 'fontsize', 14);
    title(sprintf('Segment %d', ii), 'fontsize', 18);
    VAF(yt,ysim)
end
eps_save('segrot')
for ii = 1:4
    tidx=segments(ii,1)*fs:(segments(ii,2)*fs);
    t=matrix(1,tidx)';
    subplot(220 + ii);
    xlim([t(1500), t(2500)]);
end
eps_save('segrotzoom')


%% PART E: Estimating in the time-domain

% filter to get rid of the low frequent drift
[bfil,afil]=butter(2,2/(.5*fs),'high');
u0f=filtfilt(bfil,afil,u0);
y0f=filtfilt(bfil,afil,y0);

for ii=1:4
    tidx=segments(ii,1)*fs:(segments(ii,2)*fs);
    t=matrix(1,tidx)';
    ut=detrend(u0f(tidx));
    yt=detrend(y0f(tidx));

    % LSQ-routine
    p0= [0.003 0.01    5   -0.1   0.04  30]; %
    lb0=[0.001 0.001   1  -1e2    0.02  10];
    ub0=[0.005 1     100  -1e-2   0.08  50];
    options = optimset('lsqnonlin');
    options = optimset(options, 'display', 'iter', 'TolX', 1e-3,...
        'diffminchange', 1e-2, 'maxfuneval', 200, 'TolFun', 1e-6);
    [P,~,R,~,~,~,J]= lsqnonlin(@errfunMBKKvTD, p0, lb0,ub0, options,ut,yt,t);
    plsqTD(ii,:)=P;
    SEMTD(ii,:)=full(sqrt(R'*R*(diag(inv(J'*J)))/length(R)))';

    M=plsqTD(ii,1);
    B=plsqTD(ii,2);
    K=plsqTD(ii,3);
    Kv=plsqTD(ii,4);
    td=plsqTD(ii,5);
    w=plsqTD(ii,6);

    set_param('modmbkkv/invM', 'gain', num2str(1/M));
    set_param('modmbkkv/B', 'gain', num2str(B));
    set_param('modmbkkv/K', 'gain', num2str(K));
    set_param('modmbkkv/Kv', 'gain', num2str(Kv));
    set_param('modmbkkv/delay', 'delay', num2str(td));
    set_param('modmbkkv/Hact',...
        'numerator', ['[ ' num2str(w.^2) ' ]']);
    set_param('modmbkkv/Hact',...
        'denominator', ['[ 1 ' num2str(2*0.7*w) ' ' num2str(w.^2) ']']);

    [~,~,ysim] = sim('modmbkkv', t, [], [t ut]);
    figure(4)
    subplot(211)
    plot(t,ut)
    subplot(212)
    plot(t,yt,t,ysim)
    VAF(yt,ysim)
end
