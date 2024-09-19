function K_uu = k_uu( xx, yy, hyp, i )
global ModelInfo
zeta=ModelInfo.zeta;
logsigma = hyp(1);
logthetat = hyp(2);
logthetax = hyp(3);
t = xx(:,1);
x = xx(:,2);
s = yy(:,1);
y = yy(:,2);
if i == 0 || i == 1

K_uu=exp(logsigma) .* (exp(-0.5e0 .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) - 0.5e0 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) + zeta .* (0.1e1 + bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) + bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) .^ (-0.1e1 ./ 0.2e1));
elseif i == 2
  
K_uu=exp(logsigma) .* (0.5e0 .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) .* exp(-0.5e0 .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) - 0.5e0 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) + zeta .* (0.1e1 + bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) + bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) .^ (-0.3e1 ./ 0.2e1) .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) ./ 0.2e1);
elseif i == 3
                
   
K_uu=exp(logsigma) .* (0.5e0 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax) .* exp(-0.5e0 .* bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) - 0.5e0 .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) + zeta .* (0.1e1 + bsxfun(@minus,t,s') .^ 2 ./ exp(logthetat) + bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax)) .^ (-0.3e1 ./ 0.2e1) .* bsxfun(@minus,x,y') .^ 2 ./ exp(logthetax) ./ 0.2e1);
end
end

