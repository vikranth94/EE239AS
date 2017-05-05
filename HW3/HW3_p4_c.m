clc;
clear all;
close all;
s = [30 70 110 150 190 230 310 350];
f_rate_mean = [1.7031 2.0156 2.4219 2.9453 5.0547 6.1641 4.6328 2.1641];
plot(s,f_rate_mean,'rx')
%f = k0+k1*sind(x)+k2*cosd(x);
F = @(x,xdata)x(1) + x(2)*cosd(xdata - x(3));
x0 = [1 1 100];
[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,s,f_rate_mean);
hold on
t = 0:1:360;
plot(t,F(x,t))
legend('Data Points','Tuning Curve')
xlabel('Angle (Degrees)')
ylabel('Firing Rate (Hz)')
title('Tuning Curve (Question - 2g)')
hold off
c0 = x(1)
c1 = x(2)
theta0 = x(3)