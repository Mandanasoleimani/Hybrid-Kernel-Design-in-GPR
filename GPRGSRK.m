clc
close all
clear all
% the GPR based on the GSRK method
% Calls on: xcdist, halton, testfunction, pdnlGSRKMLE, DistanceMatrix, jit_chol 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long
global DM noise y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_n=10^(-14);  % variance of noise
%sigma_n=0;
nsim=500;  % number of realizations based on the simulation design
% var1=.0001;  % variance of training data
%x01=mvnrnd(mean1,var1);% Randomly choose an initial hyperparameter from
x01=10;
% normal prior distribution for numerical optimizationalgorithm in MATLAB 
% that minimize negative log normal prior-MLE
nd = 2;  % dimension of problem
mn =30;  % 30 The number of points in one dimension for training data
mne =10;  % 10 The number of points in one dimension for test data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A star shaped domain has a boundary that can be expressed as r(theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an irregular domain
% rb = @(theta) 1+0.8*(cos(4*theta));
% R = 1.8; % The maximum radius
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an irregular domain
% rb = @(theta) 1+0.1*(sin(6*theta)+sin(3*theta));
% R = 1.2; %Actually less, this is an upper bound
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an irregular domain
% rb = @(theta) (0.1e1 + cos((2 * theta)).^2);
% R = 2.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an irregular domain
% rb = @(theta) 1*(cos(3*theta)+sqrt(4-(sin(3*theta)).^2)).^(1/3);
% R = 1.7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an irregular domain
rb = @(theta) sin(theta).^2+cos(theta).^2;
R = 1.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------SET UP THE INFORMATION NEEDED FOR NODE GENERATION-------------
%
% We use Halton nodes and we need to know what number to choose
%
N0 = 200; % Used just for testing the geometry
x0 = R*(-1 + 2*halton(N0,nd));
%
% Estimate the fraction of the square covered by the domain.
%
[th,r]=cart2pol(x0(:,1),x0(:,2));
in = find(r<rb(th));
Aomega = length(in)/N0;
r0 = xcdist(x0,x0);
%
% Average distance between nearest nodes
%
h0 = sum(min(r0+eye(N0)))/N0;
%
% We will use estimated arclength to distribute boundary nodes
%
theta0=0:2*pi/N0:2*pi;
r = rb(theta0);

x0 = r.*cos(theta0);
y0 = r.*sin(theta0);
%
% Measure the Euclidian length of each theta interval
%
l = sqrt(diff(x0).^2 + diff(y0).^2);
L = [0 cumsum(l)];
%
% L(theta) is a monotone function that can be inverted to give theta(L).
% That is, we can find theta values that give a uniform distribution of
% arclength between boundary nodes through
%
uniarcl = @(lb) interp1(L,theta0,lb);
%
%------------ CREATE AND PLOT THE EVAL GRID -----------------------------
%
% Errors will be measured over a uniform polar grid.
%
Ne_t = 60;
Ne_r = 20;

th = 0:2*pi/Ne_t:2*pi;
r = 0:1/(Ne_r-1):1;
[rr,tt]=meshgrid(r,th);
rr = diag(rb(th))*rr;
[xx,yy] = pol2cart(tt,rr);
xe = [xx(:) yy(:)];
Ne = size(xe,1);

figure(1)
% mesh(rr.*cos(tt),rr.*sin(tt),0*rr)
axis equal
view(2)
hold on
plot(x0,y0,'k-')
box on
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate training nodes

% Generate interior node
N=(mn)^2;
%%%%%%%%%%%%%%%%%%%%%
% Halton point
% x = R*(-1 + 2*halton(ceil(N/Aomega),nd));
%%%%%%%%%%%%%%%%
% Uniform point
x1 = linspace(0,1,ceil(sqrt(N/Aomega)))';
[xm1,ym1] = meshgrid(x1,x1); x = R*(-1 + 2*[xm1(:) ym1(:)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[th,r]=cart2pol(x(:,1),x(:,2));
in = find(r<rb(th));
x = x(in,:);
Nreal = size(x,1);
N = Nreal;
%
% Generate nodes on the boundary. Estimate the distance between nodes
% through scaling of h0.
%
h=h0*sqrt(N0/N);
%
% Divide the full arclength into pieces of size h
%
Nb = ceil(L(end)/h) + 1;
lg = linspace(0,L(end),Nb);
lg = lg(1:end-1)'; % Last point same as first
Nb = Nb-1;
%
thb = uniarcl(lg);
xb = rb(thb).*cos(thb);
yb = rb(thb).*sin(thb);
xb = [xb,yb];
rbc = xcdist(xb,x) + 2*R*eye(size(xb,1),N);
bdist = min(rbc);
pos = find(bdist>0.5*h);
x = x(pos,:);
data = [xb;x];
X=sortrows(data);  % training data
inb=length(xb);
ind=length(x);
n=inb+ind  % The number of training data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate test nodes

% Generate interior node
Nne=(mne)^2;
%%%%%%%%%%%%%%%%%%%%%
% Halton point
% x = R*(-1 + 2*halton(ceil(N/Aomega),nd));
%%%%%%%%%%%%%%%%
% Uniform point
x11 = linspace(0,1,ceil(sqrt(Nne/Aomega)))';
[xm11,ym11] = meshgrid(x11,x11); x111 = R*(-1 + 2*[xm11(:) ym11(:)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[th1,r]=cart2pol(x111(:,1),x111(:,2));
in1 = find(r<rb(th1));
x111 = x111(in1,:);
Nreal1 = size(x111,1);
Nne = Nreal1;
%
% Generate nodes on the boundary. Estimate the distance between nodes
% through scaling of h0.
%
h1=h0*sqrt(N0/Nne);
%
% Divide the full arclength into pieces of size h
%
Nb1 = ceil(L(end)/h1) + 1;
lg1 = linspace(0,L(end),Nb1);
lg1 = lg1(1:end-1)';  % Last point same as first
Nb1 = Nb1-1;
  %
thb1 = uniarcl(lg1);
xb1 = rb(thb1).*cos(thb1);
yb1 = rb(thb1).*sin(thb1);
xb1 = [xb1,yb1];
rbc1 = xcdist(xb1,x111) + 2*R*eye(size(xb1,1),Nne);
bdist1 = min(rbc1);
pos1 = find(bdist1>0.5*h);
x111 = x111(pos1,:);
data1 =[xb1;x111];
x_star=sortrows(data1);  % test data
inb1=length(xb1);
ind1=length(x111);
neval=inb1+ind1  % The number of test data
axis equal
plot(X(:,1),X(:,2),'bo',x_star(:,1),x_star(:,2),'+r');
%  xlabel('$\sqrt{N}$','interpreter','latex')
%     ylabel('$\|e\|_\infty/\|u\|_\infty$','interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DM= DistanceMatrix(X,X);  % Compute distance matrix of training data
I=eye(n);
noise=sigma_n*I;
rnd1=randn(n,1);
epsil=sqrt(sigma_n)*rnd1;
f= testfunctionsD(X);  % Evaluates the testfunction at training data
y=f+epsil;  % noisy training data
%y=dataTrain(:,3);  % noisy training data
%exact_f_star= dataTest(:,3);  % Evaluates the testfunction at test
% data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical optimization algorithm in MATLAB that minimize negative log
% MLE using initial hyperparameter as the starting value and
% obtain an estimate of the hyperparameter, you can choose fminunc or minimize
% options = optimoptions('fminunc','GradObj','on','Display','iter',...
%    'Algorithm','trust-region','Diagnostics','on','DerivativeCheck','on',...
%    'FinDiffType','central');
% ep = fminunc(@pdnlGSRKMLE,x01,options)  % optimal hyperparameter
options = optimoptions('fmincon','Display','iter','Algorithm','trust-region-reflective','DerivativeCheck','off','SpecifyObjectiveGradient',true,'Diagnostics','on');
lb = 0;
ub = 20;
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
%    ep = fmincon(@pdnlGSRKMLE,x01,A,b,Aeq,beq,lb,ub,nonlcon,options)  % optimal hyperparameter
% [ep,~,~] = minimize(5, @pdnlCSRKMLE, 200)
% [NLML,ff]=pdnlGSRKMLE((ep))
% ep=abs(ep)
 ep=1.5;
 %epp=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%K1= exp(-ep*DM).*(ep^2*DM.^2+3*ep*DM+3);
%%%%%%%%%%%%%%%Cubic
%K1=DM.^3;
%%%%%%%%%%%%%%hybrid kernel
%p=0;
%K1=(exp(-ep^2*DM.^2))+p.*(exp(-epp*DM).*(1+(epp.*DM)+(1/3).*((epp.*DM).^2)));
%K1=(1-p).*(exp(-ep^2*DM.^2))+p.*(exp(-ep*DM).*(1+(ep.*DM)+(1/3).*((ep.*DM).^2)));
%Matern Kernel(first order polynomial)
%K1=exp(-(ep.*DM)).*(1+(ep.*DM));
%Matern Kernel(cubic order polynomial)
  %K1=exp(-(ep.*DM)).*(1+(ep.*DM)+(2/5).*(ep.*DM).^2+(1/15).*(ep.*DM).^3);
%Matern32
%K1=exp(-ep*DM).*(1+(ep.*DM)+(1/3).*((ep.*DM).^2));
%K1= exp(-ep*DM).*(ep^2*DM.^2+3*ep*DM+3);
%Gaussian
K1=exp(-ep^2*DM.^2);
Ky1=K1+noise;  % Compute covariance matrix of training data
Ky2=K1;
condnoise=condest(Ky1)
cond=condest(Ky2)
% Compute condition number for covariance matrix of training
% data
%figure(2)
%spy(Ky1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DM1=DistanceMatrix(x_star,X);  % Compute distance matrix between test data
% and training data
%k_star =  exp(-ep*DM1).*(ep^2*DM1.^2+3*ep*DM1+3);
%gaussian
k_star =  exp(-ep^2*DM1.^2);
%hybrid kernel
%k_star=(exp(-ep^2*DM1.^2))+p.*(exp(-epp*DM1).*(1+(epp.*DM1)+(1/3).*((epp.*DM1).^2)));
%k_star=(1-p).*(exp(-ep^2*DM1.^2))+p.*(exp(-ep*DM1).*(1+(ep.*DM1)+(1/3).*((ep.*DM1).^2)));
%Matern Kernel(first order polynomial)
%k_star=exp(-(ep.*DM1)).*(1+(ep.*DM1));
%Matern Kernel(cubic order polynomial)
  %k_star=exp(-(ep.*DM1)).*(1+(ep.*DM1)+(2/5).*(ep.*DM1).^2+(1/15).*(ep.*DM1).^3);
%k_star=(1+(ep.*DM1)+(2/5).*(ep.*DM1).^2+(1/15).*(ep.*DM1).^3).*exp(-(ep.*DM1));
%matern32
%k_star=exp(-ep*DM1).*(1+(ep.*DM1)+(1/3).*((ep.*DM1).^2));
%cubic
%k_star=DM1.^3;
% Compute 
% covariance matrix between test data and training data
L1 = chol(Ky1 ,'lower');
% pL=pinv(L1);
m_f_s= L1'\(L1\y);
% m_f_s= pL'*pL*y;
mean_f_star= k_star*m_f_s;  % mean of predictive distribution
DM2=DistanceMatrix(x_star,x_star);  % Compute distance matrix of test data
%gaussian
k_star_star=exp(-ep^2*DM2.^2);
%cubic
%k_star_star=DM2.^3;
%matern32
%k_star_star= exp(-ep*DM2).*(ep^2*DM2.^2+3*ep*DM2+3);
%k_star_star=exp(-ep*DM2).*(1+(ep.*DM2)+(1/3).*((ep.*DM2).^2));
%Matern Kernel(first order polynomial)
%k_star_star=exp(-(ep.*DM2)).*(1+(ep.*DM2));
%Matern Kernel(cubic order polynomial)
   % k_star_star=exp(-(ep.*DM2)).*(1+(ep.*DM2)+(2/5).*(ep.*DM2).^2+(1/15).*(ep.*DM2).^3);
%k_star_star=(1+(ep.*DM2)+(2/5).*(ep.*DM2).^2+(1/15).*(ep.*DM2).^3).*exp(-(ep.*DM2));
%hybrid kernel
%k_star_star=(exp(-ep^2*DM2.^2))+p.*(exp(-epp*DM2).*(1+(epp.*DM2)+(1/3).*((epp.*DM2).^2)));
%k_star_star=(1-p).*(exp(-ep^2*DM2.^2))+p.*(exp(-ep*DM2).*(1+(ep.*DM2)+(1/3).*((ep.*DM2).^2)));
% Compute 
% covariance matrix of test data
v=L1\k_star';
% v=pL*k_star';
V_f_star=k_star_star-v'*v;  % variance of predictive distribution
p_mean_f_star=zeros(neval,nsim);
exact_f_star= testfunctionsD(x_star);  % Evaluates the testfunction at test
% data
error2 = exact_f_star - mean_f_star;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RMSerror2=rms(error2)  % root mean squared error of predictive
for j1=1:nsim
p_mean_f_star(:,j1)=mvnrnd(mean_f_star,V_f_star);  % predictive of test data
end
 
pre_mean_f_star=mean(p_mean_f_star,2);  % mean of predictive
error1 = exact_f_star - pre_mean_f_star; 
error2 = exact_f_star - mean_f_star;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RMSerror=rms(error1)  % root mean squared error of predictive
abserror=abs(error1); 
maxerror=max(abserror)  % maximum absolute error of predictive
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pre_var_f_star=var(p_mean_f_star,0,2);  % variance of predictive
lower=pre_mean_f_star-1.96*sqrt(pre_var_f_star);  % lower bound of predictive
upper=pre_mean_f_star+1.96*sqrt(pre_var_f_star);  % upper bounds of predictive
pre_std_f_star=std(p_mean_f_star,0,2);  % standard deviation of predictive
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XX=reshape(x_star(:,1),1,[]);
Y=reshape(x_star(:,2),1,[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
plot3(XX,Y,exact_f_star,'*');
zlabel('Noise-free solution','FontSize',15)
grid on; set(gca,'Fontsize',15);
axis([-1 1 -1 1 -2 2])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4);
plot3(XX,Y,(pre_mean_f_star),'*')
zlabel('Mean solution','FontSize',15)
grid on; set(gca,'Fontsize',15);
axis([-1 1 -1 1 -2 2])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5);
plot3(XX,Y,(abserror),'*')
zlabel('Absolute error','FontSize',15)
grid on; set(gca,'Fontsize',15);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6);
plot3(XX,Y,(pre_std_f_star),'*')
zlabel('Standard deviation','FontSize',15)
grid on; set(gca,'Fontsize',15);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7);
plot3(XX,Y,(lower),'*')
zlabel('Lower bound','FontSize',15)
grid on; set(gca,'Fontsize',15);
axis([-1 1 -1 1 -2 2])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(8);
plot3(XX,Y,(upper),'*')
zlabel('Upper bound','FontSize',15)
grid on; set(gca,'Fontsize',15);
axis([-1 1 -1 1 -2 2])
%gaussian cond noise
%x=[7.8479e+04,7.5354e+06,9.0918e+08,9.6373e+10,9.4351e+12,1.0531e+15,1.8343e+17];
%y= [1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
% z=-log10(y);
% xx=log10(x);
%gaussian RMSE & noise
%x=[0.0424,0.0038,3.9232e-04,4.2894e-05,5.5119e-06,4.5723e-07,6.1164e-08];
%xx=log10(x);
%y= [1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
% z=-log10(y);
%plot(z,xx,'-o');
%matern cond & noise
%x=[1.8725e+05,1.7274e+07,1.1109e+09,1.1590e+10,1.5304e+10,1.5322e+10,1.5322e+10];
%xx=-log10(x);
%y= [1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
% z=-log10(y);
%plot(z,xx,'-o');
%matern RMSE & noise
%x=[0.1432,0.0096,7.0842e-04,1.6438e-04,1.3977e-04,1.3840e-04,1.3701e-04];
%xx=-log10(x);
%y= [1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
% z=-log10(y);
%plot(z,xx,'-o');
%%C=[7.8479e+04,7.5558e+06,8.9378e+08,9.3658e+10,8.4856e+12,3.8247e+14,3.5002e+15,7.8481e+04,8.0032e+06,8.8113e+08,9.3033e+10,4.9886e+12,1.5019e+14,3.6023e+14,7.8497e+04,7.4693e+06,9.2219e+08,7.2963e+10,4.3556e+12,3.0557e+13,3.9850e+13,7.8654e+04,7.8378e+06,7.7781e+08,4.8387e+10,1.8755e+12,3.9567e+12,4.0072e+12,8.0221e+04,8.1890e+06,6.0862e+08,3.6240e+10,3.6152e+11,4.2364e+11,4.0931e+11];
% cond=log10(C);
%r=[0.03607,0.0036,4.4834e-04,4.9393e-05,4.4509e-06,6.9282e-07,8.8039e-08,0.03760,0.0040,4.4995e-04,5.5343e-05,5.2854e-06,8.6043e-07,1.6160e-07,0.0368,0.0036,4.3085e-04,4.7537e-05,7.3136e-06,9.6875e-07,3.5777e-07,0.0363,0.0040,4.4254e-04,4.9074e-05,7.3209e-06,1.2956e-06,1.0211e-06,0.0329,0.0039,4.8295e-04,6.6012e-05,9.7327e-06,2.7807e-06,2.3487e-06];
%rms=log10(r);
%N=[1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14,1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14, 1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14,1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14,1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
% Noise=log10(N);
%p=-6.*ones(1,7);
%s=-4.*ones(1,7);
%u=-2.*ones(1,7);
%w=-5.*ones(1,7);
%d=-3.*ones(1,7);
%Wei=[p,w,s,d,u];
%Weigh=meshgrid(Wei);
%[X,Y]=meshgrid(Noise,rms);
%figure(3);
%surf(X,Y,Weigh)
%plot3(Noise,cond,Wei,'-o')
%plot(Noise(1:7),cond(1:7),'-o')
%plot(Noise(1:7),cond(8:14),'-o')
%plot(Noise(1:7),cond(15:21),'-o')
%plot(Noise(1:7),cond(22:28),'-o')
%plot(Noise(1:7),cond(29:35),'-o')