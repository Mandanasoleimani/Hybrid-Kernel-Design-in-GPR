clear all
clc
close all
format short
addpath ../../../Utilities
addpath ../../../Utilities/modeling

global ModelInfo

%% Setup
n_f = 10;
n_u = 10;
D = 2;
lb = zeros(1,D);
ub = ones(1,D);

noise_u = 0.0;
noise_f = 0.0;

%% Generate Data
rng(1111)
% Data on u(x)
ModelInfo.x_u = bsxfun(@plus,lb,bsxfun(@times,   lhsdesign(n_u,D)    ,(ub-lb)));
ModelInfo.y_u = u(ModelInfo.x_u) + noise_u*randn(n_u,1);

rng(2222)
% Data on f(x)
ModelInfo.x_f = bsxfun(@plus,lb,bsxfun(@times,   lhsdesign(n_f,D)    ,(ub-lb)));
ModelInfo.y_f = f(ModelInfo.x_f) + noise_f*randn(n_f,1);
ModelInfo.zeta=1e-1;
%1e-1;

ModelInfo.jitter = 1e-16;
 %% Optimize model
% hyp = [logsigma logtheta alpha beta logsigma_n_u logsigma_n_f]
hyp = log([1 .01 .01 exp(2)]);

[ModelInfo.hyp,~,~] = minimize(hyp, @likelihoodoptimal, -10000);

fprintf(1,'alpha = %f\n\n', ModelInfo.hyp(4));

%% Make Predictions & Plot results
save_plots = 0;
time_dep = 1;
PlotResults

%% Post-processing
rmpath ../../../Utilities
rmpath ../../../Utilities/modeling