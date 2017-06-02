theta_0 = [0; 60; 120; 180; 240; 300];      % angles recorded
y_g = [25; 40; 70; 30; 10; 15];         % firing rates at angles recorded
coeff_g = [coeff(1,:); 1 sqrt(3)/2 0.5; coeff(2,:); 1 0 -1; coeff(3,:); 1 -sqrt(3)/2 0.5];
% coefficients of y from set of linear equations

k_g = coeff_g\y_g;                      % solve for coefficients of k
display(k_g)

f_g = k_g(1) + k_g(2)*sind(x) + k_g(3)*cosd(x);         % tuning curve

figure(2)                               % plot tuning curve from 0 to 360
                                        % as well as data points
plot(theta_0,y_g, 'o', x,f_g)
legend('y(0), y(60), y(120), y(180), y(240), and y(360)', 'Tuning Curve')
xlabel('Angle')
ylabel('Firing Rate')
title('Tuning Curve (2g)')

theta_0_g = atan2d(k_g(2),k_g(3));            % calculate theta_0, c_0, and c_1
c_0_g = k_g(1);
c_1_g = k_g(2)/sind(theta_0_g);

display(theta_0_g)
display(c_0_g)
display(c_1_g)
