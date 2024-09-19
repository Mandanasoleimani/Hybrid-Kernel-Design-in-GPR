epp=[0.5 1 1.5];
sigma_nn=[10^(-2) 10^(-4) 10^(-6) 10^(-8) 10^(-10) 10^(-12) 10^(-14)];
[X,Y] = meshgrid(epp,sigma_nn);
%xxx=X(1,:);
yyy=Y(:,1);
%Cond gaussi ep 0.1 0.5 1 1.5
Cond_gaussi_01=[2.6498e+05,3.0168e+07,3.6005e+09,3.6790e+11,3.7623e+13,4.0689e+15,2.7200e+19];
Cond_gaussi_05=[2.8386e+05,2.6752e+07,3.4029e+09,3.0139e+11,3.5522e+13,3.2708e+15,2.0384e+20];
Cond_gaussi_1=[2.0584e+05,2.2318e+07,2.3085e+09,2.3047e+11,2.6734e+13,2.6192e+15,1.6637e+18];
Cond_gaussi_15=[1.2187e+05,1.30327e+07,1.4497e+09,1.5728e+11,1.6825e+13,1.7478e+15,5.2126e+17];
Cond_gaussi=[Cond_gaussi_01; Cond_gaussi_05; Cond_gaussi_1; Cond_gaussi_15]';
Cond_matern_01=[2.5631e+05,2.9727e+07,2.8579e+09,2.6140e+11,2.3417e+13,1.9459e+15,7.6377e+16];
Cond_matern_05=[2.3377e+05,2.6232e+07,2.3572e+09,1.9670e+11,7.3393e+12,2.0268e+13,1.5547e+13];
Cond_matern_1=[2.4428e+05,2.2225e+07,2.0081e+09,1.05108e+11,5.8056e+11,5.1691e+11,6.0182e+11];
Cond_matern_15=[2.0705e+05,2.0273e+07,1.6045e+09,4.8413e+10,7.1948e+10,7.2295e+10,7.2299e+10];
%cond gaussi 1.5 matern 1
figure
hold on
plot(log10(yyy),log10(Cond_gaussi_15),'b')
plot(log10(yyy),log10(Cond_matern_1),'r')
xlabel('\sigma')
ylabel('log_{10}(Condition number with noise)')
legend('Condition number for Gaussian (\epsilon=1.5) and Matern 3/2 (\epsilon=1) kernels with different \sigma','Location','southwest')

%gaussi with ep 0.5 1 1.5
RMS_gaussi_01=[0.7021,0.7075,0.7050,0.6821,0.6778,0.6451];
RMS_gaussi_05=[0.6900,0.6357,0.4721,0.3376,0.1487,0.0908];
RMS_gaussi_1=[0.4862,0.2404,0.0648,0.0111,0.0023,2.2331e-04];
RMS_gaussi_15=[0.1524,0.0140,0.0013,1.0390e-04,1.0e-05,1.0e-06];
%matern
RMS_matern_01=[0.6990,0.7060,0.6738,0.2348,0.0161,7.8345e-04,1.0065e-04];
RMS_matern_05=[0.6974,0.3757,0.0340,0.0018,9.4471e-05,1.1539e-04,1.1592e-04];
RMS_matern_1=[0.5379,0.0694,0.0038,1.6559e-04,1.2197e-04,1.2329e-04,1.2327e-04];
RMS_matern_15=[0.2873,0.0201,0.0011,1.4763e-04,1.3031e-04,1.2971e-04,1.3075e-04];
figure
hold on
plot(log10(yyy(1:6)),log10(RMS_gaussi_05),'b')
plot(log10(yyy(1:7)),log10(RMS_matern_05),'r')
xlabel('\sigma')
ylabel('log_{10}(RMS error with noise for \epsilon= 0.5)')
legend('RMS error for Gaussian and Matern 3/2 kernels with different \sigma','Location','southwest')
%%
figure
hold on
plot(log10(yyy(1:6)),log10(RMS_gaussi_1),'b')
plot(log10(yyy(1:7)),log10(RMS_matern_1),'r')
xlabel('\sigma')
ylabel('log_{10}(RMS error with noise for \epsilon= 1)')
legend('RMS error for Gaussian and Matern 3/2 kernels with different \sigma','Location','southwest')
%%
figure
hold on
axis([-14 -2,-8 0])
plot(log10(yyy(1:6)),log10(RMS_gaussi_15),'b')
plot(log10(yyy(1:7)),log10(RMS_matern_15),'r')
xlabel('\sigma')
ylabel('log_{10}(RMS error with noise for \epsilon= 1.5)')
legend('RMS error for Gaussian and Matern 3/2 kernels with different \sigma','Location','southwest')
 plot(log10(yyy(6)),log10(RMS_gaussi_15(6)),'r*');
%gaussi 1.5 matern 1
figure
hold on
plot(log10(yyy(1:6)),log10(RMS_gaussi_15),'b')
plot(log10(yyy(1:7)),log10(RMS_matern_1),'r')
xlabel('\sigma')
ylabel('log_{10}(RMS error with noise for \epsilon= 1.5)')
legend('RMS error for Gaussian (\epsilon=1.5) and Matern 3/2 (\epsilon=1) kernels with different \sigma','Location','southwest')
%3d plot gaussi
epp=[0.1 0.5 1 1.5];
RMS_gaussi=[RMS_gaussi_01; RMS_gaussi_05; RMS_gaussi_1; RMS_gaussi_15]';
[X,Y] = meshgrid(epp,sigma_nn(1:6));
figure
surf(X,log10(Y),log10(RMS_gaussi))
xlabel('\epsilon'),ylabel('log_{10}(\sigma)'),zlabel('log_{10}(RMS_gaussi with noise)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%hybrid
legend('RMS error for Gaussian kernel with different \epsilon and noise','Location','southwest')
figure
[X,Y] = meshgrid(epp,sigma_nn);
%surf(X,log10(Y),log10(Cond_gaussi))
xlabel('\epsilon'),ylabel('log_{10}(\sigma)'),zlabel('log_{10}(Cond_gaussi) with noise')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%hybrid
legend('Condition number for Gaussian kernel with different \epsilon and noise','Location','southwest')
rho=[10^(-6) 10^(-5) 10^(-4) 10^(-3) 10^(-2)];
[X,Y] = meshgrid(rho,sigma_nn);
rms_hybrid_ep05_p6=[0.6900,0.6357,0.4717,0.3182,0.0550,0.0027,3.4156e-04];
rms_hybrid_ep1_p6=[0.4862,0.2407,0.0648,0.0111,0.0020,1.1364e-04,4.5108e-06];
rms_hybrid_ep15_p6=[0.1546,0.0142,0.0013,9.5025e-05,9.2208e-06,7.4047e-07,8.1469e-08];
rms_hybrid_ep05_p5=[0.6907,0.6355,0.4675,0.2145,0.0132,6.9095e-04,8.3782e-05];
rms_hybrid_ep1_p5=[0.4875,0.2404,0.0647,0.0107,0.0011,6.6347e-05,1.60018e-05];
rms_hybrid_ep15_p5=[0.1528,0.0142,0.0013,9.9362e-05,8.4259e-06,7.9070e-07,1.1715e-07];
rms_hybrid_ep05_p4=[0.6915,0.6335,0.4304,0.0636,0.0028,1.2570e-04,7.1741e-05];
rms_hybrid_ep1_p4=[0.4892,0.2403,0.0641,0.0087,4.2035e-04,4.0367e-05,2.1707e-05];
rms_hybrid_ep15_p4=[0.1556,0.0146,0.0013,9.9517e-05,7.2868e-06,8.3090e-07,1.5806e-07];
rms_hybrid_ep05_p3=[0.6899,0.6147,0.2535,0.0138,6.5995e-04,6.9181e-05,8.0536e-05];
rms_hybrid_ep1_p3=[0.4870,0.2394,0.0584,0.0056,2.6199e-04,3.4297e-05,3.8300e-05];
rms_hybrid_ep15_p3=[0.1520,0.0145,0.0014,9.1643e-05,7.1800e-06,9.2749e-07,3.5382e-07];
rms_hybrid_ep05_p2=[0.6881,0.4895,0.0708,0.0034,1.1197e-04,9.9340e-05,1.0097e-04];
rms_hybrid_ep1_p2=[0.4869,0.2319,0.0336,0.0021,1.2075e-04,1.06519e-04,1.0742e-04];
rms_hybrid_ep15_p2=[0.1532,0.0140,0.0013,8.5469e-05,9.5575e-06,2.4193e-06,2.2375e-06];
rms_hybrid_05=[rms_hybrid_ep05_p6;rms_hybrid_ep05_p5;rms_hybrid_ep05_p4;rms_hybrid_ep05_p3;rms_hybrid_ep05_p2]';
rms_hybrid_1=[rms_hybrid_ep1_p6;rms_hybrid_ep1_p5;rms_hybrid_ep1_p4;rms_hybrid_ep1_p3;rms_hybrid_ep1_p2]';
rms_hybrid_15=[rms_hybrid_ep15_p6;rms_hybrid_ep15_p5;rms_hybrid_ep15_p4;rms_hybrid_ep15_p3;rms_hybrid_ep15_p2]';
xxx=X(1,:);
yyy=Y(:,1);
%%%2D
figure
plot(log10(yyy(1:6)),log10(RMS_gaussi_05))
hold on
plot(log10(yyy),log10(rms_hybrid_ep05_p6))
plot(log10(yyy),log10(rms_hybrid_ep05_p5))
plot(log10(yyy),log10(rms_hybrid_ep05_p4))
plot(log10(yyy),log10(rms_hybrid_ep05_p3))
plot(log10(yyy),log10(rms_hybrid_ep05_p2))
legend('condition number for Gaussian and hybrid kernel with different \rho, \sigma and \epsilon=0.5.','Location','southwest')
hold off
%%
figure
plot(log10(yyy(1:6)),log10(RMS_gaussi_1))
hold on
plot(log10(yyy),log10(rms_hybrid_ep1_p6))
plot(log10(yyy),log10(rms_hybrid_ep1_p5))
plot(log10(yyy),log10(rms_hybrid_ep1_p4))
plot(log10(yyy),log10(rms_hybrid_ep1_p3))
plot(log10(yyy),log10(rms_hybrid_ep1_p2))
legend('condition number for Gaussian and hybrid kernel with different \rho, \sigma and \epsilon=1.','Location','southwest')
hold off
%%
figure
plot(log10(yyy(1:6)),log10(RMS_gaussi_15))
hold on
plot(log10(yyy),log10(rms_hybrid_ep15_p6))
plot(log10(yyy),log10(rms_hybrid_ep15_p5))
plot(log10(yyy),log10(rms_hybrid_ep15_p4))
plot(log10(yyy),log10(rms_hybrid_ep15_p3))
plot(log10(yyy),log10(rms_hybrid_ep15_p2))
legend('condition number for Gaussian and hybrid kernel with different \rho, \sigma and \epsilon=1.5.','Location','southwest')
hold off
%%3D
figure
surf(log10(X),log10(Y),log10(rms_hybrid_05))
xlabel('log_{10}(\rho)'),ylabel('log_{10}(\sigma)'),zlabel('log_{10}(RMSE_ using_ hybrid_ kernel_ with_ noise')
legend('RMS error for hybrid kernel with different \sigma, \rho and \epsilon= 0.5.','Location','best')
figure
surf(log10(X),log10(Y),log10(rms_hybrid_1))
xlabel('\rho'),ylabel('log_{10}(\sigma)'),zlabel('log_{10}(RMSE_ using_ hybrid_ kernel_ with_ noise')
legend('RMS error for hybrid kernel with different \sigma, \rho and \epsilon= 1.','Location','best')
figure
surf(log10(X),log10(Y),log10(rms_hybrid_15))
xlabel('\rho'),ylabel('log_{10}(\sigma)'),zlabel('log_{10}(RMSE_ using_ hybrid_ kernel_ with_ noise')
legend('RMS error for hybrid kernel with different \sigma, \rho and \epsilon= 1.5.','Location','best')
