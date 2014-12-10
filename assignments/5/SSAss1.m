clear all
clc

%% A. create system & data
Nt=1000;                                % No. Samples of simulation
dt=0.01;                                % Sample time
t=(0:Nt-1).'*dt;                        % Time vector
plotcol='g';                            % Plot color

% M=0.5;B=1;K=150;                      % System parameters
M=0.1;B=0.5;K=450;
sys0=c2d(ss(tf(1,[M B K])),dt);         % discrete system representation
A0=sys0.a;B0=sys0.b;C0=sys0.c;D0=sys0.d;

uk=zeros(Nt,1);                         % zero input
x0=[1 0]';                              % initial condition
y0=lsim(sys0,zeros(size(t)),t,x0);      % clean output
yk=y0;%+1e-2*randn(size(y0));           % output noise added

l=size(yk,2);                           % number of outputs
r=size(uk,2);                           % number of inputs

figure(1)
set(gcf,'name','Output time domain')
plot(t,y0,t,yk), hold on;
xlabel('t')
ylabel('y')
legend('y0','yk')

%% B. Data matrices
s=3;ii=0;N=200;                                     % Hankel matrix starting at yk(ii) size s x N
YsN=hankel(yk(ii+1:ii+1+s-1),yk(ii+s:ii+N+s-1));    % Output Hankel matrix
UsN=hankel(uk(ii+1:ii+1+s-1),uk(ii+s:ii+N+s-1));    % Input Hankel matrix

figure(2)
set(gcf,'name','Data Space')
plot3(YsN(1,:),YsN(2,:),YsN(3,:),plotcol,'linewidth',2), hold on

%% C. Column space Os + SVD
[U,S,V]=svd(YsN,'econ');           % Singular value decomposition
 
Us=U*S;

figure(2)
for jj=1:s
    line([0 Us(1,jj)],[0 Us(2,jj)],[0 Us(3,jj)],'color',plotcol,'linestyle','--','linewidth',2)
end
 
figure(3)
semilogy(S,plotcol,'marker','*'), hold on
xlim([0 s+1])

%% D. Retrieve system from SVD
n=[];                               % System order
Un=U(:,1:n);                        % Reduced output singular vectors
Vn=V(:,1:n);                        % Reduced input singular vectors
Sn=S(1:n,1:n);                      % Reduced singular value matrix
SnVtn=Sn*Vn';

Aid=Un(1:(s-1)*l,:)\Un(l+1:s*l,:);
Bid=SnVtn(:,1:r);
Cid=Un(1:l,:);
Did=yk(1);
sysid=ss(Aid,Bid,Cid,Did,dt);       % Identified system out of SVD


%% E. plot end result
yi=lsim(sysid,zeros(size(t)),t,x0);

figure(1)
plot(t,yi,'r--')
legend('y0','yk','yi')
