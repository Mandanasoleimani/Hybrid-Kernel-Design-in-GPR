clc
%clear all
% the GPR based on the GSRK method
% Calls on: xcdist, halton, testfunction, pdnlGSRKMLE, DistanceMatrix, jit_chol 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long
global DM noise y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_n=10^(-2);  % variance of noise
nsim=500;  % number of realizations based on the simulation design
x01=1;  % choose an initial hyperparameter for numerical 
% optimization algorithm in MATLAB that minimize negative log MLE
nd = 7;  % dimension of problem
% data=co2mmmlo{:,[3 4]};
% % [C,ia] = unique(data(:,1:8));
% % data = data{ia,:};
% cv=cvpartition(size(data,1),'HoldOut',.2);
% idx=cv.test;
% dataTrain=data(~idx,:);
% dataTrain=sortrows(dataTrain);
% dataTest=data(idx,:);
% dataTest=sortrows(dataTest);
% X=dataTrain(:,1);  % training data
% n=length(X);  % The number of training data
% x_star=dataTest(:,1);  % test data
% neval=length(x_star);  % The number of test data
data=CO2EmissionsCanada(:,[4 5 8 9 10 11 12]);
[C,ia] = unique(data(:,1:6));
data = data{ia,:};
data = zscore(data);
cv=cvpartition(size(data,1),'HoldOut',.3);
idx=cv.test;
dataTrain=data(~idx,:);
dataTrain=sortrows(dataTrain);
dataTest=data(idx,:);
dataTest=sortrows(dataTest);
X=dataTrain(:,1:6);  % training data
n=length(X);  % The number of training data
x_star=dataTest(:,1:6);  % test data
neval=length(x_star);  % The number of test data
%plot(X(:,1),X(:,2),'bo',x_star(:,1),x_star(:,2),'+r');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DM= DistanceMatrix(X,X);  % Compute distance matrix of training data
I=eye(n);
noise=sigma_n*I;
y=dataTrain(:,7);  % noisy training data
exact_f_star= dataTest(:,7);  % Evaluates the testfunction at test
% data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical optimization algorithm in MATLAB that minimize negative log
% MLE using initial hyperparameter as the starting value and
% obtain an estimate of the hyperparameter, you can choose fminunc or minimize
% options = optimoptions('fminunc','GradObj','on','Display','iter',...
%    'Algorithm','trust-region','Diagnostics','on','DerivativeCheck','on',...
%    'FinDiffType','central');
% ep = fminunc(@pdnlGSRKMLE,x01,options)  % optimal hyperparameter
options = optimoptions('fmincon','Display','iter','Algorithm','trust-region-reflective','SpecifyObjectiveGradient',true,'Diagnostics','on');
lb = 0;
ub = 2;
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
ep=1.5;
  % ep = fmincon(@pdnlGSRKMLE,x01,A,b,Aeq,beq,lb,ub,nonlcon,options)  % optimal hyperparameter
% [ep,~,~] = minimize(5, @pdnlCSRKMLE, 200)
%[NLML,ff]=pdnlGSRKMLE((ep))
% ep=abs(ep)
%hybrid kernel
%p=10^(-2);
%K1=(exp(-ep^2*DM.^2))+p.*(exp(-ep*DM).*(1+(ep.*DM)+(1/3).*((ep.*DM).^2)));
%K1=(1-p).*(exp(-ep^2*DM.^2))+p.*(exp(-ep*DM).*(1+(ep.*DM)+(1/3).*((ep.*DM).^2)));
%Matern Kernel(cubic order polynomial)
 % K1=exp(-(ep.*DM)).*(1+(ep.*DM)+(2/5).*(ep.*DM).^2+(1/15).*(ep.*DM).^3);
%Matern32
%K1=(1+((sqrt(3).*DM)/ep)).*exp((-sqrt(3).*DM)./(ep));
%K1=exp(-ep*DM).*(1+(ep.*DM)+(1/3).*((ep.*DM).^2));
%K1= exp(-ep*DM).*(ep^2*DM.^2+3*ep*DM+3);
%Gaussian
K1=exp(-ep^2*DM.^2);
%Matern Kernel(first order polynomial)
%K1=exp(-(ep.*DM)).*(1+(ep.*DM));
%squared exponential
%K1=exp((-1/2).*(1/(ep))^2*DM.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ky1=K1+noise;  % Compute covariance matrix of training data
Ky2=K1;
condnoise=condest(Ky1)
cond=condest(Ky2)
% data
%figure(2)
%spy(Ky1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DM1=DistanceMatrix(x_star,X);  % Compute distance matrix between test data
% and training data
%matern32
%k_star=(1+((sqrt(3).*DM1)/ep)).*exp((-sqrt(3).*DM1)./(ep));
%k_star=exp(-ep*DM1).*(1+(ep.*DM1)+(1/3).*((ep.*DM1).^2));
%Matern Kernel(cubic order polynomial)
  %k_star = exp(-(ep.*DM1)).*(1+(ep.*DM1)+(2/5).*(ep.*DM1).^2+(1/15).*(ep.*DM1).^3);
%hybrid kernel
%k_star=(exp(-ep^2*DM1.^2))+p.*(exp(-ep*DM1).*(1+(ep.*DM1)+(1/3).*((ep.*DM1).^2)));
  %gaussian
k_star =  exp(-ep^2*DM1.^2);
%matern first order
%k_star=exp(-(ep.*DM1)).*(1+(ep.*DM1));
%squared exponential
%k_star=exp((-1/2).*(1/(ep))^2*DM1.^2);
%k_star =  exp(-ep*DM1).*(ep^2*DM1.^2+3*ep*DM1+3);  % Compute 
% covariance matrix between test data and training data
L1 = chol(Ky1 ,'lower');
% pL=pinv(L1);
m_f_s= L1'\(L1\y);
% m_f_s= pL'*pL*y;
mean_f_star= k_star*m_f_s;  % mean of predictive distribution
DM2=DistanceMatrix(x_star,x_star);  % Compute distance matrix of test data 
%k_star_star= exp(-ep*DM2).*(ep^2*DM2.^2+3*ep*DM2+3);
%matern32
%k_star_star=exp(-ep*DM2).*(1+(ep.*DM2)+(1/3).*((ep.*DM2).^2));
%k_star_star= exp(-ep*DM2).*(ep^2*DM2.^2+3*ep*DM2+3);
%k_star_star=exp(-ep*DM2).*(1+(ep.*DM2)+(1/3).*((ep.*DM2).^2));
%Matern Kernel(cubic order polynomial)
 %k_star_star = exp(-(ep.*DM2)).*(1+(ep.*DM2)+(2/5).*(ep.*DM2).^2+(1/15).*(ep.*DM2).^3);
%hybrid kernel
%k_star_star=(exp(-ep^2*DM2.^2))+p.*(exp(-ep*DM2).*(1+(ep.*DM2)+(1/3).*((ep.*DM2).^2)));
%squared exponential
%k_star_star=exp((-1/2).*(1/(ep))^2*DM2.^2);
%gaussian
k_star_star=exp(-ep^2*DM2.^2);
%matern first order
%k_star_star=exp(-(ep.*DM2)).*(1+(ep.*DM2));
% Compute 
% covariance matrix of test data
v=L1\k_star';
% v=pL*k_star';
V_f_star=k_star_star-v'*v;  % variance of predictive distribution
p_mean_f_star=zeros(neval,nsim);
error2 = exact_f_star - mean_f_star;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RMSerror2=rms(error2)
MAPerror2=mape(mean_f_star,exact_f_star)
%%%DM1=DistanceMatrix(X,X);  % Compute distance matrix between test data
% and training data
%matern32
%k_star=(1+((sqrt(3).*DM1)/ep)).*exp((-sqrt(3).*DM1)./(ep));
%k_star=exp(-ep*DM1).*(1+(ep.*DM1)+(1/3).*((ep.*DM1).^2));
%gaussian
%%%k_star =  exp(-ep^2*DM1.^2);
%squared exponential
%k_star=exp((-1/2).*(1/(ep))^2*DM1.^2);
%k_star =  exp(-ep*DM1).*(ep^2*DM1.^2+3*ep*DM1+3);  % Compute 
% covariance matrix between test data and training data
%%%L1 = chol(Ky1 ,'lower');
% pL=pinv(L1);
%%%m_f_s= L1'\(L1\y);
% m_f_s= pL'*pL*y;
%%%mean_f_star= k_star*m_f_s;  % mean of predictive distribution
%%%DM2=DistanceMatrix(X,X);  % Compute distance matrix of test data 
%k_star_star= exp(-ep*DM2).*(ep^2*DM2.^2+3*ep*DM2+3);
%matern32
%k_star_star=exp(-ep*DM2).*(1+(ep.*DM2)+(1/3).*((ep.*DM2).^2));
%k_star_star= exp(-ep*DM2).*(ep^2*DM2.^2+3*ep*DM2+3);
%k_star_star=exp(-ep*DM2).*(1+(ep.*DM2)+(1/3).*((ep.*DM2).^2));
%squared exponential
%k_star_star=exp((-1/2).*(1/(ep))^2*DM2.^2);
%gaussian
%%%k_star_star=exp(-ep^2*DM2.^2);
% Compute 
% covariance matrix of test data
%%%v=L1\k_star';
% v=pL*k_star';
%%%V_f_star=k_star_star-v'*v;  % variance of predictive distribution
%%%p_mean_f_star=zeros(neval,nsim);
%%%error2 = y - mean_f_star;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%RMSerror2_train=rms(error2)
for j1=1:nsim
p_mean_f_star(:,j1)=mvnrnd(mean_f_star,V_f_star);  % predictive of test data
end
% exact_f_star= testfunctionsD(x_star);  % Evaluates the testfunction at test
% data
pre_mean_f_star=mean(p_mean_f_star,2);  % mean of predictive
error = exact_f_star - pre_mean_f_star; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RMSerror=rms(error)  % root mean squared error of predictive
Maperror=mape(pre_mean_f_star,exact_f_star)
abserror=abs(error); 
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
% eps .5 gaussian
% x=[5.8258e+05,8.1445e+07,9.0376e+09,1.0093e+12,9.0210e+13,9.3466e+15,2.6991e+18];
% y=[1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
% z=log10(y);
% xx=log10(x);
% plot(z,xx,'-o');
% eps .5 matern3/2
%t=[ 1.3982e+06,8.3091e+07,1.0293e+10,5.2225e+11,3.4556e+13,9.6282e+13,8.9414e+13];
% y=[1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
% z=log10(y);
% tt=log10(t);
% plot(z,tt,'-o');
%hybrid  eps .5 gaussian matern32
%s=[6.4613e+05,7.7695e+07,7.8082e+09,1.0837e+12,7.8987e+13,8.5620e+15,1.6030e+18];
%y=[1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
% z=log10(y);
% ss=log10(s);
% plot(z,ss,'-o');
% u=[5.5048e+05,7.6454e+07,9.1049e+09,9.8855e+11,7.2234e+13,5.5825e+15,2.1605e+17];
% y=[1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
% z=log10(y);
% uu=log10(u);
% plot(z,uu,'-o');
% w=[5.7327e+05,9.4159e+07,8.2569e+09,6.0867e+11,5.5543e+13,2.6649e+15,4.2654e+15];
% y=[1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
% z=log10(y);
% ww=log10(w);
% plot(z,ww,'-o');
% eps 1
% x=[2.5434e+05,3.5333e+07,4.1053e+09,3.4422e+11,3.1466e+13,2.3562e+15,2.3866e+17];
% y=[1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
% z=log10(y);
% xx=log10(x);
% plot(z,xx,'-o');
% hold on
% t=[5.3477e+05,5.5386e+07,5.0308e+09,2.8195e+11,2.1230e+12,1.9985e+12,1.6773e+12];
% y=[1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
% z=log10(y);
% tt=log10(t);
% plot(z,tt,'-o');
% s=[2.5795e+05,3.5509e+07,3.1834e+09,3.7201e+11,2.6249e+13,1.9019e+15,1.0356e+17];
% y=[1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
% z=log10(y);
% ss=log10(s);
% plot(z,ss,'-o');
% u=[2.4918e+05,3.5965e+07,2.9878e+09,3.1617e+11,2.4386e+13,1.3141e+15,9.5740e+15];
% z=log10(y);
% uu=log10(u);
% plot(z,uu,'-o');
% w=[2.4446e+05,3.6722e+07,3.3988e+09,2.4802e+11,1.5747e+13,8.4364e+13,1.2078e+14];
% y=[1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
% z=log10(y);
% ww=log10(w);
% plot(z,ww,'-o');
%eps 1.5
% x=[1.3080e+05,1.6652e+07,1.9483e+09,1.5152e+11,1.3358e+13,1.1338e+15,7.9710e+16];
% y=[1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
% z=log10(y);
% xx=log10(x);
% plot(z,xx,'-o');
% hold on
% t=[3.3077e+05,5.2625e+07,2.5963e+09,1.2903e+11,2.2952e+11,2.2025e+11,2.2194e+11];
% y=[1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
% z=log10(y);
% tt=log10(t);
% plot(z,tt,'-o');
% s=[1.3318e+05,1.7681e+07,1.6961e+09,1.6930e+11,1.3426e+13,9.0632e+14,1.1623e+16];
% y=[1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
% z=log10(y);
% ss=log10(s);
% plot(z,ss,'-o');
% w=[1.7240e+05,1.7217e+07,1.6271e+09,1.2089e+11,4.7811e+12,7.5473e+12,7.1817e+12];
% y=[1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
% z=log10(y);
% ww=log10(w);
% plot(z,ww,'-o');
%rms ep 0.5 Ga  mat  hyb6 hyb4 hyb2
% eps 0.5
% x=[0.2058,0.2261,0.6000,3.3967,22.0193,50.9889];
%y=[1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12];
% z=log10(y);
% xx=log10(x);
% plot(z,xx,'-o');
% hold on
% t=[0.2276,0.1989,0.2256,0.4462,0.3723,0.5857,0.3219];
% tt=log10(t);
% yy=[1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
%zz=log10(yy);
% plot(zz,tt,'-o');
% s=[0.2202,0.3379,0.5062,1.2068,31.8561,17.0874];
% ss=log10(s);
% plot(z,ss,'-o');
% u=[0.1970,0.5184,0.7978,11.1515,2.3886,19.5345,28.6932];
% uu=log10(u);
% plot(zz,uu,'-o');
% w=[0.2282,0.4743,0.5800,1.0252,7.7686,3.9060,3.3517];
% ww=log10(w);
% plot(zz,ww,'-o');
% eps 1
% x=[0.1947,0.4093,0.8051,2.8614,17.5473,38.9516,1.2031e+02];
% y=[1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
% z=log10(y);
% xx=log10(x);
% plot(z,xx,'-o');
% hold on
% t=[0.2292,0.2356,0.3395,0.5038,0.4906,0.9188,0.4542];
% tt=log10(t);
%yy=[1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
%zz=log10(yy);
% plot(zz,tt,'-o');
% s=[0.2609,0.3172,0.3662,1.6449,5.9747,61.0585,1.180104120743491e+02];
% ss=log10(s);
% plot(zz,ss,'-o');
% u=[0.2324,0.2784,0.5549,4.1119,7.4376,17.1330,22.9713];
% uu=log10(u);
% plot(zz,uu,'-o');
% w=[0.2022,0.2559,0.7184,4.7440,2.3681,4.6816,6.8983];
% ww=log10(w);
% plot(zz,ww,'-o');
% eps 1.5
% x=[0.2217,0.3186,1.5166,8.4933,11.0073,52.1122,73.6934];
% y=[1.0e-02,1.0e-04,1.0e-06,1.0e-08,1.0e-10,1.0e-12,1.0e-14];
% z=log10(y);
% xx=log10(x);
% plot(z,xx,'-o');
% hold on
% t=[0.2326,0.2759,0.2901,0.3019,0.3731,0.3265,0.3724];
% tt=log10(t);
% plot(zz,tt,'-o');
% s=[0.2828,0.2451,1.3962,2.2688,35.6450,60.8844,67.7139];
% ss=log10(s);
% plot(zz,ss,'-o');
% u=[0.1967,0.2580,0.7326,2.2725,12.6707,14.0919,26.0446];
% uu=log10(u);
% plot(zz,uu,'-o');
% w=[0.2750,0.3114,0.7265,2.5488,3.8967,6.2826,13.2128];
% ww=log10(w);
% plot(zz,ww,'-o');