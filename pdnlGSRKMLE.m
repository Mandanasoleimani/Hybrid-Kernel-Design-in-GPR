function [NLML,ff] = pdnlGSRKMLE(ep)
% Compute the negative log MLE and its gradient w.r.t. the hyperparameter
% Inputs:
% ep: hyperparameter of the kernel
% Calls on: jit_chol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global DM noise y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,D] = size(y);
% K1=exp(-ep*DM);
% K2=(ep^2*DM.^2+3*ep*DM+3);
% K11=K1.*K2;
K11=exp(-ep^2*DM.^2);
cond(K11)
Ky1=(K11+noise);
cond(Ky1)
L = chol(Ky1);
p1=L'\(L\y);
NLML = 0.5*y'*p1 + sum(log(diag(L)));  % Compute the negative log MLE
% NLML = 0.5*y'*p1 + sum(log(diag(L))) + log(2*pi)*n/2;
pK=L'\(L\eye(n));
% p1=pK*y;
p2=-(p1*p1'-pK);
% n1=-DM.*K11;
% n2=K1.*(2*ep*DM.^2+3*DM);
n=(-2*ep*DM.^2).*K11;
% ff=sum(sum(p2.*(n1+n2)))/2;
ff=sum(sum(p2.*(n)))/2;  % Compute the gradient of the negative log
% MLE w.r.t. the hyperparamete


