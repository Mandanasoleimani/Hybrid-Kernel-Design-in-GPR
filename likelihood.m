function [condknoisefree  condkwnoise  ]=likelihood(hyp)

global ModelInfo
x_u = ModelInfo.x_u;
x_f = ModelInfo.x_f;
zeta=ModelInfo.zeta;
% y_u = ModelInfo.y_u;
% y_f = ModelInfo.y_f;

%y=[y_u;y_f];

% jitter = ModelInfo.jitter;

 sigma_n_f = (hyp(end));
 sigma_n_u = (hyp(end-1));

n_u = size(x_u,1);
[n_f,D] = size(x_f);
n = n_u+n_f;

K_uu = k_uu(x_u, x_u, hyp(1:D+1),0);
 K_uf = k_uf(x_u, x_f, hyp(1:end-2),0);
%K_uf = k_uf(x_u, x_f, hyp(1:end),0);
K_fu = K_uf';
 K_ff = k_ff(x_f, x_f, hyp(1:end-2),0);
%K_ff = k_ff(x_f, x_f, hyp(1:end),0);
K = [K_uu K_uf;
    K_fu K_ff];
% K = K + eye(n).*jitter;
% Cholesky factorisation
% [L,p]=chol(K,'lower');
% 
% ModelInfo.L = L;
% 
% if p > 0
%     fprintf(1,'Covariance is ill-conditioned\n');
% end
% 
% %alpha = L'\(L\y);
%  pinvk=inv(K);alpha=pinvk*y;
% NLMLnoisefree = 0.5*y'*alpha + sum(log(diag(L))) + log(2*pi)*n/2;
condknoisefree=cond(K)
 K_uu = K_uu;
 K_ff = K_ff;
 K_uu = K_uu + eye(n_u).*sigma_n_u;
 K_ff = K_ff + eye(n_f).*sigma_n_f;
 K = [K_uu K_uf;
     K_fu K_ff];

%K = K + eye(n).*jitter;
% Cholesky factorisation
% [L,p]=chol(K,'lower');
% 
% ModelInfo.L = L;
% 
% if p > 0
%     fprintf(1,'Covariance is ill-conditioned\n');
% end
% 
% %alpha = L'\(L\y);
% pinvk=inv(K);alpha=pinvk*y;
% NLMLwnoise = 0.5*y'*alpha + sum(log(diag(L))) + log(2*pi)*n/2;
 condkwnoise=(cond(K))
end