%% SIPE (2014-2015): Assignment 4

clear all
close all
clc

%% Load Dataset

global t y

% Data1
load('Data1');   % output y; timevector t

% Plot data
Timefit = figure; clf;
plot(t,y,'r','linewidth',4); hold on

%% Grid search (1)

% Start Time
tic

% Upper & Lower Bounds [a,b]
lb = [0 0];
ub = [10 10];

% Resolution
ares = lb(1):0.1:ub(1);
bres = lb(2):0.1:ub(2);

% Error & Parameter Calculation
e1 = zeros(length(ares),length(bres));
for qq = 1:length(ares)
    for ss = 1:length(bres)
        par = [ares(qq) bres(ss)];
        e1(qq,ss) = errfun(par);    !! Fill in your errfun (sum of squares)
    end
end

% Error plot (watch out use e transposed)
errsurf = figure; clf;
surfc(ares,bres,e1.'); hold on;
xlabel('a'); ylabel('b'); zlabel('e')

% Find Optimum (parest & error residual)
amin = []; bmin = [];               !! Define the best solution (lowest error)
parest1 = [amin bmin];
err1 = [];                          !! Define lowest error
clear tmp

% Estimate & Plot yest
[~,yest1] = errfun(parest1);

figure(Timefit);
plot(t,yest1,'b--','linewidth',2); hold on;
legend('Measurements','Grid Search')

% Stop Time
time1 = toc;

%% Gradient Search (2)

% Start Time
tic

% Gradient Search Options (Use 'doc optimset' for more options)
options = optimset('lsqnonlin');
options = optimset(options, 'Display','iter');
options = optimset(options,'outputfcn',@storepathway);

% Set pathway variables
global pathway
pathway.x = []; pathway.y = []; pathway.z = [];

% Initial Parameter Value (First step)
par0 = [5 5];

% LSQnonlin
[parest2,~,e2,~,~,~,J]=lsqnonlin(@errfun2,par0,lb,ub,options); !! Build your own errfun2
err2 = e2*e2';  % sum of squares (=sum(e2.^2))

% Plot Pathway
figure(errsurf)
plot3(pathway.x, pathway.y, pathway.z, 'g.-', 'linewidth', 2)
plot3(pathway.x(1), pathway.y(1), pathway.z(1), 'g.', 'markersize', 25)
plot3(pathway.x(end), pathway.y(end), pathway.z(end), 'g.', 'markersize', 50)

% Estimate & Plot yest
[~,yest2] = errfun2(parest2);

figure(Timefit);
plot(t,yest2,'g--','linewidth',2); hold on;
legend('Measurements','Grid Search','Gradient Search')

% Stop Time
time2 = toc;

%% Generic Search (3)

% Start Time
tic

% Generic Algorithm Search (Check what they do with 'help gaoptimset')
options = gaoptimset;
options = gaoptimset(options,'populationsize', 50);
options = gaoptimset(options,'elitecount', 15);
options = gaoptimset(options,'crossoverfraction', 0.8);
options = gaoptimset(options,'popinitrange', [0; 8]);
options = gaoptimset(options,'plotfcn', @customplotfcn);
options = gaoptimset(options,'plotinterval', 1);

% Error Contour plot (watch out use e transposed)
figure(99); clf;
contour(ares,bres,e1.'); hold on;
daspect([1 1 1])

% Generic Algorithm
[parest3, err3] = ga(@errfun, 2,[],[],[],[],lb,ub,[],options);

% Estimate & Plot yest
[~,yest3] = errfun(parest3);

figure(Timefit);
plot(t,yest3,'k--','linewidth',2); hold on;
legend('Measurements','Grid Search','Gradient Search','Generic Algorithm')

% Stop Time
time3 = toc;

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
subplot(212)
plot(t0,u0)
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
    ylabel('Gain'); xlim([1/T0 fs/2]);
subplot(323)
    semilogx(fv0,unwrap(angle(H0))*180/pi,'linewidth',1,'linestyle','--'); hold on;
    semilogx(mfv0,unwrap(angle(mH0))*180/pi,'r','linewidth',2);
    ylabel('Phase [deg]'); xlim([1/T0 fs/2]); ylim([-360 360]);
subplot(325)
    semilogx(mfv0,mCoh0,'linewidth',2); hold on;
    xlabel('Frequency [Hz]'); ylabel('Coherence');
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
        loglog(fv,abs(H),'color',colset(ii,:),'linewidth',1,'linestyle','--'); hold on;
        loglog(mfv,abs(mH),'color',colset(ii,:),'linewidth',2);
        ylabel('Gain'); xlim([1/T0 fs/2]);
    subplot(324)
        semilogx(fv,unwrap(angle(H))*180/pi,'color',colset(ii,:),'linewidth',1,'linestyle','--'); hold on;
        semilogx(mfv,unwrap(angle(mH))*180/pi,'color',colset(ii,:),'linewidth',2);
        ylabel('Phase [deg]'); xlim([1/T0 fs/2]);
    subplot(326)
        semilogx(mfv,mCoh,'linewidth',2,'color',colset(ii,:)); hold on;
        xlabel('Frequency [Hz]'); ylabel('Coherence');
        xlim([1/T0 fs/2]); ylim([0 1.1]);
    figure(1)
    subplot(212)
    for jj=1:2
        line([segments(ii,jj) segments(ii,jj)],[3 -3],'color',colset(ii,:))
    end
end

%% PART B : Estimate MBK model
clear all
close all
clc

% for each segment
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
    p0=[0.002   0.1  5 ];
    lb0=[0.001  0.01 1 ];
    ub0=[0.005  1    20];

    %Estimate on these frequencies:
    estfreq=[1/T fs];
    estidx=find(mfv<estfreq(2)&mfv>estfreq(1));

    options = optimset('lsqnonlin');
    options = optimset(options, 'display', 'iter','TolX',1e-15,'diffminchange',1e-4,'maxfuneval',500,'TolFun',1e-12);
    [P,~,R,~,~,~,J]= lsqnonlin(@errfunMBK, p0, lb0,ub0, options,mH(estidx),mfv(estidx),mCoh(estidx));
    plsq(ii,:)=P;
    SEM(ii,:)=full(sqrt(R'*R*(diag(inv(J'*J)))/length(R)))';

    [E,Hid]=errfunMBK(plsq(ii,:),mH,mfv,mCoh);

    % Plot FRF
    figure(3)
    set(gcf,'Name','MBK Estimates')
    subplot(311)
        loglog(mfv(estidx),abs(mH(estidx)),'color',colset(ii,:),'linewidth',2), hold on;
        loglog(mfv(estidx),abs(Hid(estidx)),'color',colset(ii,:),'linewidth',2,'linestyle','--');
        ylabel('Gain'); %xlim([1/T0 fs/2]);
    subplot(312)
        semilogx(mfv(estidx),(angle(mH(estidx)))*180/pi,'color',colset(ii,:),'linewidth',2), hold on;
        semilogx(mfv(estidx),(angle(Hid(estidx)))*180/pi,'color',colset(ii,:),'linewidth',2,'linestyle','--');
        ylabel('Phase [deg]'); %xlim([1/T0 fs/2]);
    subplot(313)
        semilogx(mfv(estidx),mCoh(estidx),'linewidth',2,'color',colset(ii,:),'linestyle','--'); hold on;
        xlabel('Frequency [Hz]'); ylabel('Coherence');
        ylim([0 1.1]); %         xlim([1/T0 fs/2]);
end
%% PART C:Initial values for Kv model
%     p=[M B K Kv td w]
%     p0=  [0.003 0.01    5   -0.1   0.04  30];
%     lb0=[0.001 0.001   1  -1e2    0.02  10];
%     ub0=[0.005 1     100  -1e-2   0.08  50];

%% PART D Determining the fit in the time domain.

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
    set_param('modmbkkv/Hact', 'numerator', ['[ ' num2str(w.^2) ' ]']);
    set_param('modmbkkv/Hact', 'denominator', ['[ 1 ' num2str(2*0.7*w) ' ' num2str(w.^2) ']']);

    [~,~,ysim] = sim('modmbkkv', t, [], [t ut]);
    figure(4)
    subplot(211)
    plot(t,ut)
    subplot(212)
    plot(t,yt,t,ysim)
    VAF(yt,ysim)
end

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
    options = optimset(options, 'display', 'iter','TolX',1e-3,'diffminchange',1e-2,'maxfuneval',200,'TolFun',1e-6);
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
    set_param('modmbkkv/Hact', 'numerator', ['[ ' num2str(w.^2) ' ]']);
    set_param('modmbkkv/Hact', 'denominator', ['[ 1 ' num2str(2*0.7*w) ' ' num2str(w.^2) ']']);

    [~,~,ysim] = sim('modmbkkv', t, [], [t ut]);
    figure(4)
    subplot(211)
    plot(t,ut)
    subplot(212)
    plot(t,yt,t,ysim)
    VAF(yt,ysim)
end
