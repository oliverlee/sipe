close all
clear all
clc

axfs = 18;
lgndfs = 16;

%% A. create system & data
Nt=1000;                            % No. Samples
dt=0.01;                            % Sample time
t=(0:Nt-1).'*dt;                    % Time vector
plotcol='b';                        % Plot color

M=0.5;B=1;K=150;                    % System parameters
% M=0.1;B=0.5;K=450;
sys0=c2d(ss(tf(1,[M B K])),dt);     % discrete system representation
A0=sys0.a;B0=sys0.B;C0=sys0.c;D0=sys0.d;

uk=[randn(Nt,1)];                   % white noise input
y0=lsim(sys0,uk,t);                 % clean output
yk=y0;% + 1e-3*randn(size(y0));       % output noise added

l=size(yk,2);                       % number of outputs
r=size(uk,2);                       % number of inputs

figure(1)
plot(t,y0,t,yk), hold on;
xlabel('t')
ylabel('y')
legend('y0','yk')

%% B. Data matrices
% Hankel matrix starting at yk(ii) size s x N
s=20;ii=0;N=300;
YsN=zeros(l*s,N);
UsN=zeros(r*s,N);

for ii=1:s
    YsN((ii-1)*l+1:ii*l,:)=yk(ii:N+ii-1,:)';
    UsN((ii-1)*r+1:ii*r,:)=uk(ii:N+ii-1,:)';
end
%% C. Column space Os + SVD
[U,S,V]=svd(YsN,'econ');                % Singular value decomposition

figure(2)
colset = cool(4);
set(gcf,'name','Data Space')
plot3(YsN(1,:),YsN(2,:),YsN(3,:),'color', colset(2, :),'linewidth',2), hold on

PsN=eye(N)-UsN'*inv(UsN*UsN')*UsN;     % Orthogonal projection on UsN
YsNP=YsN*PsN;                          
plot3(YsNP(1,:),YsNP(2,:),YsNP(3,:), 'color', colset(4, :),'linewidth',2), hold on
[U,S,V]=svd(YsNP,'econ');              % Singular value decomposition
xlabel('y(k) [m]', 'fontsize', axfs);
ylabel('y(k + 1) [m]', 'fontsize', axfs);
zlabel('y(k + 2) [m]', 'fontsize', axfs);
l2 = legend({'YsN', 'YsNP'});
set(l2, 'fontsize', lgndfs);
eps_save('dataspace2d', figure(2));

%plot colums of U
%Us=U*S;
%figure(2)
%for jj=1:s
%    line([0 Us(1,jj)],[0 Us(2,jj)],[0 Us(3,jj)],'color',plotcol,'linestyle','--','linewidth',2)
%end

figure(3)
semilogy(S,plotcol,'marker','*'), hold on

%% D. Retrieve system from SVD
n = 2;
Un=U(:,1:n);                        % Reduced output singular vectors
Vn=V(:,1:n);                        % Reduced input singular vectors
Sn=S(1:n,1:n);                      % Reduced singular value matrix
SnVtn=Sn*Vn';

Aid=Un(1:(s-1)*l,:)\Un(l+1:s*l,:);
Cid=Un(1:l,:);

Bid=SnVtn(:,1:r);
Did=yk(1);

%% E. Least squares routine to estimate B and D out of A,C,u,y
Phik=zeros(N*l,n+(n+l)*r);
for kk=1:N-1
    clear dummykron
    dummykron=zeros((kk)*l,n*r);
    dummysum=reshape([ones(1,kk);zeros(l-1,kk)],l*kk,1);
    dummyprod=toeplitz(dummysum(1:l),dummysum(1:l*kk));
    for tau=0:kk-1
        dummykron(tau*l+1:(tau+1)*l,:)=kron(uk(tau+1,:),Cid*Aid^(kk-tau-1));
    end
    Phik(kk*l+1:(kk+1)*l,:)=[Cid*Aid^kk dummyprod*dummykron kron(uk(kk+1,:),eye(l))];
end
% Least squares solution min(yk-Phik*theta)
theta=((1/N)*(Phik'*Phik))\(1/N*Phik')*reshape(yk(1:N,:),l*N,1);

% parameters out of parameter vector
x0=theta(1:n);
Bid=reshape(theta(n+1:n+n*r),n,r);
Did=reshape(theta(n+n*r+1:n+n*r+r*l),l,r);
%% F. plot end result
sysid=ss(Aid,Bid,Cid,Did,dt);
yi=lsim(sysid,uk,t);

figure(1)
plot(t,yi,'r--')
legend('y0','yk','yi')

%%
ykid = yi;
yid = lsim(sysid_noiseless, uk, t);
figure(1); clf
colset = cool(4);
plot(t, y0, 'color', colset(1, :), 'linewidth', 2); hold on;
plot(t, yk, 'color', colset(2, :), 'linewidth', 2);
plot(t, yid, '--', 'color', colset(3, :), 'linewidth', 2);
plot(t, ykid, '--', 'color', colset(4, :), 'linewidth', 2); hold off;
l1 = legend('y0', 'yk', 'y_{id}', 'yk_{id}');
set(l1, 'fontsize', lgndfs);
xlabel('t [s]', 'fontsize', axfs);
ylabel('y [m]', 'fontsize', axfs);
eps_save('output2e', figure(1));

