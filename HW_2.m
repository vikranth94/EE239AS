% EE239AS.2 HW #2, Problem 2

%% Part E

coeff = [1 0 1; 1 sqrt(3)/2 -0.5; 1 -sqrt(3)/2 -0.5];          
% linear equations for coefficient of k
y_coeff = inv(coeff);
% solve for the coefficients of y
display(y_coeff)

%% Part F

% fit f(t) to the given data

theta = [0; 120; 240];          % angles recorded
y = [25; 70; 10];               % firing rates at angles recorded

k = coeff\y;                    % k coefficients
display(k)

x = linspace(0, 360, 100);      % for plotting
f = k(1) + k(2)*sind(x) + k(3)*cosd(x);         % tuning curve

theta_0 = atan2d(k(2),k(3));    % calculate theta_0, c_0, and c_1
c_0 = k(1);
c_1 = k(2)/sind(theta_0);

display(theta_0)
display(c_0)
display(c_1)

figure(1)                       % plot the tuning curve from 0 to 360 degrees
                                % as well as data points
plot(theta,y, 'o', x,f)
legend('y(0), y(120), and y(240)', 'Tuning Curve')
xlabel('Angle')
ylabel('Firing Rate')
title('Tuning Curve (2f)')

%% Part G

% fit f(t) to the given additional data

t_g = [0; 60; 120; 180; 240; 300];      % angles recorded
y_g = [25; 40; 70; 30; 10; 15];         % firing rates at angles recorded
coeff_g = [coeff(1,:); 1 sqrt(3)/2 0.5; coeff(2,:); 1 0 -1; coeff(3,:); 1 -sqrt(3)/2 0.5];
% coefficients of y from set of linear equations

k_g = coeff_g\y_g;                      % solve for coefficients of k
display(k_g)

f_g = k_g(1) + k_g(2)*sind(x) + k_g(3)*cosd(x);         % tuning curve

figure(2)                               % plot tuning curve from 0 to 360
                                        % as well as data points
plot(t_g,y_g, 'o', x,f_g)
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
