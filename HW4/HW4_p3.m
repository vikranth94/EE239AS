% EE239AS Homework 4

%% Problem 3: Simulated Neural Data

ps3_data = importdata('ps4_simdata.mat');

% 20x3 struct
% rows = data point
% columns = class

%% Part A: 2D Plot

D = 2;          % number of neurons
n_class = size(ps3_data,2);
n_trial = size(ps3_data,1);
data = cell(1, n_class);

for i = 1:n_trial
    for j = 1:n_class
        data{1,j} = [data{1,j}, ps3_data(i,j).x];
    end
end

figure(1)
plot(data{1,1}(1,:),data{1,1}(2,:),'xr', data{1,2}(1,:),data{1,2}(2,:), '+g', ...
     data{1,3}(1,:),data{1,3}(2,:), 'ob')

title('2D Firing Rate Data Representation')
xlabel('Number of Spikes (Neuron 1)')
ylabel('Number of Spikes (Neuron 2)')
legend('k=1','k=2','k=3')
axis([0 20 0 20])

%% Part B: ML Parameters

% Model (i) Gaussian, Shared Covariance
% Model (ii) Gaussian, Class Specific Covariance
% Model (iii) Poisson


