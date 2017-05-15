% EE239AS Homework 4

clc
clear
close all
%% Problem 3: Simulated Neural Data

ps3_data = importdata('ps4_simdata.mat');

% 20x3 struct
% rows = data point
% columns = class

%% Part A: 2D Plot

D_trial = 2;
n_class = size(ps3_data,2);
n_trial = size(ps3_data,1);
data = cell(1, n_class);

for i = 1:n_trial
    for j = 1:n_class
        data{1,j} = [data{1,j}, ps3_data(i,j).x];
        % organize data into cell for easier access
    end
end

figure(1)
plotData(data)
% plot neurons in each class in different colors
title('2D Firing Rate Data Representation')
xlabel('Number of Spikes (Neuron 1)')
ylabel('Number of Spikes (Neuron 2)')
legend('Class 1','Class 2','Class 3')

%% Part B: ML Parameters

% Model (i) Gaussian, Shared Covariance
N_k = n_trial;
N = N_k*n_class;
P_Ck = N_k/(n_class*N_k);
% calculate the prior probabilities of each class (equal for all classes)

mu_i = zeros(D_trial, n_class);
S_k_i = cell(1, n_class);
sigma_i = zeros(D_trial, D_trial);

for i = 1:n_class
    mu_i(:,i) = 1/(N_k)*sum(data{1,i},2);
    cov_trial_i = zeros(D_trial, D_trial);
    for j = 1:n_trial
        cov_trial_i = cov_trial_i + (data{1,i}(:,j)-mu_i(:,i))*(data{1,i}(:,j)-mu_i(:,i))';
        % sum the (x-mu)*(x-mu)' matrices for each trial 
    end
    S_k_i{i} = 1/N_k * cov_trial_i;
    % calculate the S_k for each class and store into cell
    sigma_i = sigma_i + N_k/N * S_k_i{i};
    % calculate sigma (weighted sum of S_k)
end

fprintf('Model (i) Gaussian, Shared Covariance\n------------------------\n\n')
disp('Probability of Each Class:')
disp(P_Ck)
disp('Means:')
disp(mu_i)
disp('Covariance Matrix:')
disp(sigma_i)

% Model (ii) Gaussian, Class Specific Covariance

% Class probabilities and mean are the same as Model (i).
% The covariance matrices are specific to each class, as opposed to the
% weighted sum in Model (i).
fprintf('Model (ii) Gaussian, Class Specific Covariance\n------------------------\n\n')
disp('Probability of Each Class:')
disp(P_Ck)
disp('Means:')
disp(mu_i)
disp('Covariance Matrix (Class 1):')
disp(S_k_i{1})
disp('Covariance Matrix (Class 2):')
disp(S_k_i{2})
disp('Covariance Matrix (Class 3):')
disp(S_k_i{3})

% Model (iii) Poisson

% Class probabilities and mean firing rate are the same as Model (i).
fprintf('Model (iii) Poisson, Class Specific Covariance\n------------------------\n\n')
disp('Probability of Each Class:')
disp(P_Ck)
disp('Mean Firing Rates:')
disp(mu_i)

%% Part C: Firing Rate Means

figure(2)
plotData(data)
title('2D Firing Rate: Means')
xlabel('Number of Spikes (Neuron 1)')
ylabel('Number of Spikes (Neuron 2)')
legend('Class 1','Class 2','Class 3')

hold on
% plot means of each class (same for each model)
plotMeans(mu_i)
hold off

%% Part D, E: Means, Covariances, and Decision Boundaries

% Model (i) Gaussian, Shared Covariance

figure(3)
plotData(data)
title('Model (i) Firing Rate: Means, Covariances, and Decision Boundaries')
xlabel('Number of Spikes (Neuron 1)')
ylabel('Number of Spikes (Neuron 2)')
legend('Class 1','Class 2','Class 3')

hold on
% plot means of each class (same for each model)
plotMeans(mu_i)

% plot covariance ellipses for each class (shared covariance)
plotContour(mu_i(:,1)',sigma_i,'r');
plotContour(mu_i(:,2)',sigma_i,'g');
plotContour(mu_i(:,3)',sigma_i,'b');

hold off

% Model (ii) Gaussian, Class Specific Covariance

figure(4)
plotData(data)
title('Model (ii) Firing Rate: Means, Covariances, and Decision Boundaries')
xlabel('Number of Spikes (Neuron 1)')
ylabel('Number of Spikes (Neuron 2)')
legend('Class 1','Class 2','Class 3')

hold on
% plot means of each class (same for each model)
plotMeans(mu_i)
plotContour(mu_i(:,1)',S_k_i{1},'r');
plotContour(mu_i(:,2)',S_k_i{2},'g');
plotContour(mu_i(:,3)',S_k_i{3},'b');

hold off

%%

for i = 1:n_class
    
    w_k(:,i) = inv(sigma_i)*mu_i(:,i);
    w_k0(i) = ((-1/2) * mu_i(:,i)'*inv(sigma_i)*mu_i(:,i) + log(P_Ck));

end

w_1 = [w_k(:,2) - w_k(:,1), w_k(:,3) - w_k(:,2), w_k(:,3) - w_k(:,1)];
w_0 = [w_k0(2) - w_k0(1), w_k0(3) - w_k0(2), w_k0(3) - w_k0(1)];

% coefficients for between class comparisons
% 1st column: Class 1 (red) and Class 2 (green)
% 2nd column: Class 3 (blue) and Class 2 (green)
% 3rd column: Class 3 (blue) and Class 1 (red)

x = 0:0.1:20;
y_1_2 = (-w_1(1,1)*x - w_0(1))/(w_1(2,1));
y_3_2 = (-w_1(1,2)*x - w_0(2))/(w_1(2,2));
y_3_1 = (-w_1(1,3)*x - w_0(3))/(w_1(2,3));

figure(6)

plotData(data)
hold on
plotMeans(mu_i)
plotContour(mu_i(:,1)',S_k_i{1},'r');
plotContour(mu_i(:,2)',S_k_i{2},'g');
plotContour(mu_i(:,3)',S_k_i{3},'b');

plot(x,y_1_2)
plot(x,y_3_2)
plot(x,y_3_1)
hold off
