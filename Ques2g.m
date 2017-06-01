clc;
clear all;
close all;
theta = [0 60 120 180 240 300];
y = [25 40 70 30 10 15];
plot(theta,y,'rx')
%f = k0+k1*sind(x)+k2*cosd(x);
F = @(x,xdata)x(1) + x(2)*cosd(xdata - x(3));
x0 = [1 1 100];
[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,theta,y);
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