function [c0, c1, theta0] = tcFit(theta, y, plotFlag)

% take out the center 0, 0 target 

%f = k0+k1*sind(x)+k2*cosd(x);
F = @(x,xdata)x(1) + x(2)*cosd(xdata - x(3));
x0 = [1 1 100];
[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,theta,y);
hold on

c0 = x(1);
c1 = x(2);
theta0 = x(3);

if plotFlag == 1
    figure
    plot(theta,y,'rx')
    hold on
    t = 0:360;
    plot(t,F(x,t))
    legend('Data Points','Tuning Curve')
    xlabel('Angle (Degrees)')
    ylabel('Firing Rate (Hz)')
    title('Part A: Tuning Curve')
    hold off
end
end

