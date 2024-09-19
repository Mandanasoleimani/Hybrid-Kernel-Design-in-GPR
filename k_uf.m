function K_uf = k_uf(xx, yy, hyp,i)

logsigma = hyp(1);
logthetat = hyp(2);
logthetax = hyp(3);
alpha = hyp(4);

t = xx(:,1);
x = xx(:,2);
s = yy(:,1);
y = yy(:,2);
global ModelInfo
zeta=ModelInfo.zeta;
if i == 0 || i == 1
    
    K_uf = exp(logsigma) .* (0.10e1 .* bsxfun(@minus,t,s') ./ exp(logthetat) .* exp(-0.5e0 .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) - 0.5e0 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) + zeta .* (0.1e1 + bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) + bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) .^ (-0.3e1 ./ 0.2e1) .* bsxfun(@minus,t,s') ./ exp(logthetat)) - alpha .* exp(logsigma) .* (-0.10e1 ./ exp(logthetax) .* exp(-0.5e0 .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) - 0.5e0 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) + 0.100e1 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax) .^ 2 .* exp(-0.5e0 .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) - 0.5e0 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) + 0.3e1 .* zeta .* (0.1e1 + bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) + bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) .^ (-0.5e1 ./ 0.2e1) .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax) .^ 2 - zeta .* (0.1e1 + bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) + bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) .^ (-0.3e1 ./ 0.2e1) ./ exp(logthetax));
elseif i == 2
    
    K_uf = exp(logsigma) .* (-0.10e1 .* bsxfun(@minus,t,s') ./ exp(logthetat) .* exp(-0.5e0 .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) - 0.5e0 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) + 0.50e0 .* bsxfun(@minus,t,s') .^ 3 ./ exp(logthetat) .^ 2 .* exp(-0.5e0 .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) - 0.5e0 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) + 0.3e1 ./ 0.2e1 .* zeta .* (0.1e1 + bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) + bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) .^ (-0.5e1 ./ 0.2e1) .* bsxfun(@minus,t,s') .^ 3 ./ exp(logthetat) .^ 2 - zeta .* (0.1e1 + bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) + bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) .^ (-0.3e1 ./ 0.2e1) .* bsxfun(@minus,t,s') ./ exp(logthetat)) - alpha .* exp(logsigma) .* (-0.50e0 ./ exp(logthetax) .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) .* exp(-0.5e0 .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) - 0.5e0 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) + 0.500e0 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax) .^ 2 .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) .* exp(-0.5e0 .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) - 0.5e0 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) + 0.15e2 ./ 0.2e1 .* zeta .* (0.1e1 + bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) + bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) .^ (-0.7e1 ./ 0.2e1) .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax) .^ 2 .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) - 0.3e1 ./ 0.2e1 .* zeta .* (0.1e1 + bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) + bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) .^ (-0.5e1 ./ 0.2e1) ./ exp(logthetax) .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat));
elseif i == 3
    
    K_uf = exp(logsigma) .* (0.50e0 .* bsxfun(@minus,t,s') ./ exp(logthetat) .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax) .* exp(-0.5e0 .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) - 0.5e0 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) + 0.3e1 ./ 0.2e1 .* zeta .* (0.1e1 + bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) + bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) .^ (-0.5e1 ./ 0.2e1) .* bsxfun(@minus,t,s') ./ exp(logthetat) .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) - alpha .* exp(logsigma) .* (0.10e1 ./ exp(logthetax) .* exp(-0.5e0 .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) - 0.5e0 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) - 0.250e1 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax) .^ 2 .* exp(-0.5e0 .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) - 0.5e0 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) + 0.500e0 .* bsxfun(@minus,x,y') .^ 4 ./ exp(logthetax) .^ 3 .* exp(-0.5e0 .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) - 0.5e0 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) + 0.15e2 ./ 0.2e1 .* zeta .* (0.1e1 + bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) + bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) .^ (-0.7e1 ./ 0.2e1) .* bsxfun(@minus,x,y') .^ 4 ./ exp(logthetax) .^ 3 - 0.15e2 ./ 0.2e1 .* zeta .* (0.1e1 + bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) + bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) .^ (-0.5e1 ./ 0.2e1) .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax) .^ 2 + zeta .* (0.1e1 + bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) + bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) .^ (-0.3e1 ./ 0.2e1) ./ exp(logthetax));
elseif i == 4
    
    K_uf = -exp(logsigma) .* (-0.10e1 ./ exp(logthetax) .* exp(-0.5e0 .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) - 0.5e0 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) + 0.100e1 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax) .^ 2 .* exp(-0.5e0 .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) - 0.5e0 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) + 0.3e1 .* zeta .* (0.1e1 + bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) + bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) .^ (-0.5e1 ./ 0.2e1) .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax) .^ 2 - zeta .* (0.1e1 + bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) + bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) .^ (-0.3e1 ./ 0.2e1) ./ exp(logthetax));
end

