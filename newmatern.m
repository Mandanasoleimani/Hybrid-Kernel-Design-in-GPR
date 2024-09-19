clc
%clear all
% the GPR based on the GSRK method
% Calls on: xcdist, halton, testfunction, pdnlGSRKMLE, DistanceMatrix, jit_chol 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long
global DM noise y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_n=0.01;  % variance of noise
nsim=500;  % number of realizations based on the simulation design
% var1=.0001;  % variance of training data
%x01=mvnrnd(mean1,var1);% Randomly choose an initial hyperparameter from
x01=10;
% normal prior distribution for numerical optimizationalgorithm in MATLAB 
% that minimize negative log normal prior-MLE
nd = 2;  % dimension of problem
mn =30;  % The number of points in one dimension for training data
mne =10;  % The number of points in one dimension for test data
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
 ep=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K1= exp(-ep*DM).*(ep^2*DM.^2+3*ep*DM+3);
%K1=exp(-ep^2*DM.^2);
Ky1=K1+noise;  % Compute covariance matrix of training data
Ky2=K1;
condest(Ky1)
condest(Ky2)
% Compute condition number for covariance matrix of training
% data
%figure(2)
%spy(Ky1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DM1=DistanceMatrix(x_star,X);  % Compute distance matrix between test data
% and training data
k_star =  exp(-ep*DM1).*(ep^2*DM1.^2+3*ep*DM1+3); 
%k_star =  exp(-ep^2*DM1.^2);
% Compute 
% covariance matrix between test data and training data
L1 = chol(Ky1 ,'lower');
% pL=pinv(L1);
m_f_s= L1'\(L1\y);
% m_f_s= pL'*pL*y;
mean_f_star= k_star*m_f_s;  % mean of predictive distribution
DM2=DistanceMatrix(x_star,x_star);  % Compute distance matrix of test data
%k_star_star=exp(-ep^2*DM2.^2);
k_star_star= exp(-ep*DM2).*(ep^2*DM2.^2+3*ep*DM2+3);  % Compute 
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