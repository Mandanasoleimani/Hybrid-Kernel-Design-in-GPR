% tf = testfunctionsD(x)
% Evaluates one of the testfunctions below at s-dimensional points

function testfunction = testfunctionsD(x)

[N,s] = size(x);
% Franke type function
%a1 = repmat([7 sqrt(10)],N,floor((s+1)/2));
%a2 = repmat([7 3],N,floor((s+1)/2));
%a3 = repmat([4 7],N,floor((s+1)/2));
%testfunction = 0.75*exp(-sum((9*x-2).^2,2)/4) + 0.75*exp(-sum(((9*x+1)./a1(:,1:s)).^2,2)) ...
%+ 0.5*exp(-sum((9*x-a2(:,1:s)).^2,2)/4) - 0.2*exp(-sum((9*x-a3(:,1:s)).^2,2));

% tensor-product multi-linear function 
%testfunction = 4^s*prod(x.*(1-x),2);
    
% s-D sinc function
%testfunction = prod(sinc(x),2);

% tensor-product weighted
%alpha = .1;
%testfunction = prod(1-alpha*x.*(1-x),2);

% borehole function
%rw=0.05 + 0.1*x(:,1);
%if s<2; r=(50000+100)/2; else r=100 + (50000-100)*x(:,2); end
%if s<3; Tu=(115600-63070)/2; else Tu=63070 + (115600-63070)*x(:,3); end
%if s<4; Tl=(116+63.1)/2; else Tl=63.1 + (116-63.1)*x(:,4); end
%if s<5; Hu=(1110+990)/2; else Hu=990 + (1110-990)*x(:,5); end
%if s<6; Hl=(820+700)/2; else Hl=700 + (820-700)*x(:,6); end
%if s<7; L=(1680+1120)/2; else L=1120 + (1680-1120)*x(:,7); end
%if s<8; Kw=(12045+9855)/2; else Kw=9855 + (12045-9855)*x(:,8); end
%logrrw=log(r./rw);
%testfunction = 2*pi*Tu.*(Hu-Hl)./(logrrw.*(1+2*L.*Tu./(logrrw.*rw.*rw.*Kw) + Tu./Tl));

% robotic arm function
%z = floor(s/2);
%theta = (2*pi)*x(:,1:z);
%len = x(:,(z+1):(2*z));
%sumtheta = cumsum(theta,2);
%u = sum(len.*cos(sumtheta),2);
%v = sum(len.*sin(sumtheta),2);
%testfunction = sqrt(u.*u+v.*v);

% product function
%alpha = .1;
%testfunction = prod(1+alpha*x,2);

% only 2D
%testfunction = (tanh(9*(x(:,2)-x(:,1)))+1)/(tanh(9)+1);
%f1 = 0.75*exp(-((9*x(:,1)-2).^2+(9*x(:,2)-2).^2)/4);
%f2 = 0.75*exp(-((9*x(:,1)+1).^2/49+(9*x(:,2)+1).^2/10));
%f3 = 0.5*exp(-((9*x(:,1)-7).^2+(9*x(:,2)-3).^2)/4);
%f4 = 0.2*exp(-((9*x(:,1)-4).^2+(9*x(:,2)-7).^2));
%testfunction = f1+f2+f3-f4;
% testfunction = 5*(x(:,2)).^4+2*(x(:,1)).^4;
%testfunction = atan(2*(x(:,2).^2+2*x(:,1)));
%testfunction = cos(x(:,2))+sin(x(:,1));
%testfunction = 1-alpha*( x(:,1).* x(:,2) - x(:,1) .* x(:,2) .^ 2 - x(:,1) .^ 2 .* x(:,2) + x(:,1) .^ 2 .* x(:,2) .^ 2);
%testfunction = x(:,1).* x(:,2) - x(:,1) .* x(:,2) .^ 2 - x(:,1) .^ 2 .* x(:,2) + x(:,1) .^ 2 .* x(:,2) .^ 2;
%testfunction = exp(( x(:,1) + x(:,2)));
%testfunction = exp(-x(:,2)).*sin(x(:,1));
%testfunction = sin(x(:,1)).*cos(x(:,2));
%testfunction =0.25e2 ./ (0.25e2 + (x(:,1) - 0.2e0) .^ 2 + (2 * x(:,2) .^ 2));
% testfunction =cos(pi * x(:,2)/ 0.2e1) .* sin(pi *x(:,1));
%testfunction=x(:,2) .* sin(pi * x(:,1)) + x(:,1) .* cos(pi * x(:,2));
%%%%%%%%%%%%%%%%%%%%%%
% testfunction = sin(pi*x(:,1)).*cosh(x(:,2))-cos(pi*x(:,1)).*sinh(x(:,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%
% testfunction=16*x(:,1).*(1-x(:,2)).*x(:,2).*(1-x(:,1));
testfunction = cos(2*pi*(x(:,2)+x(:,1)));
%testfunction = sin(pi*x(:,1)).*cosh(x(:,2))-cos(pi*x(:,1)).*sinh(x(:,2));
%testfunction =exp(-0.1e2 * (x(:,1)- 0.5e0) .^ 2 - 0.1e2 * (x(:,2) - 0.5e0) .^ 2) / 0.5e1 ;
%testfunction =0.25e1 ./ (0.25e1 + (x(:,1) - 0.8e0) .^ 2 + ( (x(:,2)-.2) .^ 2));