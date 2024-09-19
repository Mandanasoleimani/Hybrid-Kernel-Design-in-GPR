%% Pre-processing
% clear all
clc
close all

addpath ../../../Utilities
addpath ../../../Utilities/modeling

global ModelInfo

%% Setup

D = 2;
lb = zeros(1,D);
ub = ones(1,D);
%  ModelInfo.jitter = eps;
noise_u = 0;
noise_f = 0;
%epvec=[1:50];
% for n_u=1:length(epvec)
%     for n_f=1:length(epvec)
%         i=n_u;j=n_f;
n_f =80;
n_u = 80;
%% Generate Data
rng(1111)
% Data on u(x)
ModelInfo.x_u = bsxfun(@plus,lb,bsxfun(@times,   lhsdesign(n_u,D)    ,(ub-lb)));
% ModelInfo.y_u = u(ModelInfo.x_u) + noise_u*randn(n_u,1);
ModelInfo.y_u = u(ModelInfo.x_u);

rng(2222)
% Data on f(x)
ModelInfo.x_f = bsxfun(@plus,lb,bsxfun(@times,   lhsdesign(n_f,D)    ,(ub-lb)));
% ModelInfo.y_f = f(ModelInfo.x_f) + noise_f*randn(n_f,1);
ModelInfo.y_f = f(ModelInfo.x_f);

%% Optimize model
% hyp = [logsigma logtheta alpha beta logsigma_n_u logsigma_n_f]
N=50;
ModelInfo.zeta=1e-3;
%1e-0;
zeta=ModelInfo.zeta;
%epvec=[logspace(-4,0,N/2),logspace(0,1,N/2)];
 epvec=linspace(.1,1,N);
 epvecc=linspace(.5,1.5,N);
for i=1:length(epvec)
    for j=1:length(epvec)        
hyp = [ log([.1 epvec(i) epvec(i)]) epvecc(j) noise_u noise_f];
[Ksigmafree(i,j)  KWithnoise(i,j) ]=likelihood(hyp);
    end
end
%epvec=1./epvec;
[X,Y] = meshgrid(epvec,epvecc);
% figure
% subplot(1,2,1),surf(X,Y,NLML_sigmafree)
% xlabel('\epsilon_x'),ylabel('\epsilon_t'),zlabel('NLML_ Withoutnoise')
% subplot(1,2,2),surf(X,Y,abs(log10(NLML_sigmafree)))
% xlabel('\epsilon_x'),ylabel('\epsilon_t'),zlabel('log(NLML)_ Withoutnoise')
% figure
% subplot(1,2,1),surf(X,Y,NLML_Withnoise)
% xlabel('\epsilon_x'),ylabel('\epsilon_t'),zlabel('NLML_ Withnoise')
% subplot(1,2,2),surf(X,Y,abs(log10(NLML_Withnoise)))
% xlabel('\epsilon_x'),ylabel('\epsilon_t'),zlabel('log(NLML)_ Withnoise')
figure
subplot(1,2,1),surf(X,Y,log10(Ksigmafree))
xlabel('\epsilon'),ylabel('\alpha'),zlabel('log_{10}(Cond(K))_ using_ hybrid_ kernel')
subplot(1,2,2),surf(X,Y,log10(KWithnoise))
xlabel('\epsilon'),ylabel('\alpha'),zlabel('log_{10}(Cond(K))_ using_ Gaussian_ kernel')
% figure
% subplot(1,2,1),surf(X,Y,Ksigmafree)
% xlabel('\epsilon'),ylabel('\alpha'),zlabel('Cond(K)_ noise_ free')
% subplot(1,2,2),surf(X,Y,KWithnoise)
% xlabel('\epsilon'),ylabel('\alpha'),zlabel('Cond(K)_ With_ noise')
colormap jet
view([-.2 -1 1.3])
%%
% figure
% subplot(1,2,1),[M,c]=contour(X,Y,NLML_sigmafree);c.LineWidth = 2;xlabel('\epsilon_x');ylabel('\epsilon_t');l = colorbar;set(get(l,'Ylabel'),'String','NLML_ Withoutnoise');
% subplot(1,2,2),[M,c]=contour(X,Y,abs(log(NLML_sigmafree)));c.LineWidth = 2;xlabel('\epsilon_x');ylabel('\epsilon_t');l = colorbar;set(get(l,'Ylabel'),'String','log(NLML)_ Withoutnoise');
% figure
% subplot(1,2,1),[M,c]=contour(X,Y,NLML_Withnoise);c.LineWidth = 2;xlabel('\epsilon_x');ylabel('\epsilon_t');l = colorbar;set(get(l,'Ylabel'),'String','NLML_ Withnoise');
% subplot(1,2,2),[M,c]=contour(X,Y,abs(log10(NLML_Withnoise)));c.LineWidth = 2;xlabel('\epsilon_x');ylabel('\epsilon_t');l = colorbar;set(get(l,'Ylabel'),'String','log(NLML)_ Withnoise');
figure
subplot(1,2,2),[M,c]=contour(X,Y,log10(KWithnoise));c.LineWidth = 2;xlabel('\epsilon');ylabel('\alpha');l = colorbar;set(get(l,'Ylabel'),'String','log_{10}(Cond(K))_ using_ Gaussian_ kernel');
% subplot(1,2,1),[M,c]=contour(X,Y,Ksigmafree);c.LineWidth = 2;xlabel('\epsilon_x');ylabel('\epsilon_t');l = colorbar;set(get(l,'Ylabel'),'String','Cond(K)_ Withoutnoise');
subplot(1,2,1),[M,c]=contour(X,Y,log10(Ksigmafree));c.LineWidth = 2;xlabel('\epsilon');ylabel('\alpha');l = colorbar;set(get(l,'Ylabel'),'String','log_{10}(Cond(K))_ using_ hybrid_ kernel');
% figure
% subplot(1,2,1),[M,c]=contour(X,Y,KWithnoise);c.LineWidth = 2;xlabel('\epsilon_x');ylabel('\epsilon_t');l = colorbar;set(get(l,'Ylabel'),'String','Cond(K)_ Withnoise');
% subplot(1,2,2),[M,c]=contour(X,Y,log10(KWithnoise));c.LineWidth = 2;xlabel('\epsilon_x');ylabel('\epsilon_t');l = colorbar;set(get(l,'Ylabel'),'String','log(Cond(K))_ Withnoise');

