% EE239AS Homework 4

clc
clear
close all

%% Problem 3: Simulated Neural Data

ps4_data = importdata('ps4_realdata.mat');

% 20x3 struct
% rows = data point
% columns = class

%% Part A: ML Parameters

train_data = ps4_data.train_trial;
test_data = ps4_data.test_trial;
% since sizes of test and train data sets are the same, extract number of
% spikes in the same loop

D_trial = 97;           % number of neurons per trial
n_class = size(train_data,2);
n_trial = size(train_data,1);
n_spikes_train = cell(1, n_class);
n_spikes_test = n_spikes_train;

for i = 1:n_trial
    for j = 1:n_class
        n_spikes_train{1,j} = [n_spikes_train{1,j}, sum(train_data(i,j).spikes(:,351:550),2)];
        n_spikes_test{1,j} = [n_spikes_test{1,j}, sum(test_data(i,j).spikes(:,351:550),2)]; 
        % only include the spikes in the middle of the trial due to
        % experimental procedure
    end
end

%% Part A: Model(i) Classification

% Model (i) Gaussian, Shared Covariance

N_k = n_trial;
N = N_k*n_class;
P_Ck = N_k/(n_class*N_k);
% calculate the prior probabilities of each class (equal for all classes)

mu_i = zeros(D_trial, n_class);
S_k_i = cell(1, n_class);
sigma_i = zeros(D_trial, D_trial);

for i = 1:n_class
    mu_i(:,i) = 1/(N_k)*sum(n_spikes_train{1,i},2);
    cov_trial_i = zeros(D_trial, D_trial);
    for j = 1:n_trial
        cov_trial_i = cov_trial_i + (n_spikes_train{1,i}(:,j)-mu_i(:,i))*(n_spikes_train{1,i}(:,j)-mu_i(:,i))';
        % sum the (x-mu)*(x-mu)' matrices for each trial 
    end
    S_k_i{i} = 1/N_k * cov_trial_i;
    % calculate the S_k for each class and store into cell
    sigma_i = sigma_i + N_k/N * S_k_i{i};
    % calculate sigma (weighted sum of S_k)
end

k = zeros(n_trial, n_class);
errors = zeros(1,n_class);

for n = 1:n_class
    xy = n_spikes_test{n}';
    % matrix of neuron firing rates across all trials (97x91) transposed
    for j = 1:n_trial
        for i = 1:n_class
            k(j,i) = log(P_Ck)+ mu_i(:,i)'*inv(sigma_i)*xy(j,:)'-0.5*mu_i(:,i)'*inv(sigma_i)*mu_i(:,i);
        end
    end
    [m,idx] = max(k, [], 2);
    errors(n)= length(idx(idx~=n));     
    % number of incorrect classifications (not matching the specified reach
    % angle)
end
accuracy = 1-sum(errors)/(D_trial*n_class);
% calculate accuracy of classification
%% Part B: Model (ii) Classification

% Model (ii) Gaussian, Class Specific Covariance

% Some neurons do not spike enough across all trials, and therefore make
% the covariance matrix non positive definite due to their small values. We 
% can remove these in order to ensure a positive definite sigma matrix. 

n_spikes_total = zeros(D_trial, n_class);
for i = 1:n_class
    % sum the number of spikes for each neuron across the entire train
    n_spikes_class = sum(n_spikes_train{1,i},2);
    % sum the number of spikes for each neuron across all trials
    n_spikes_total(:,i) = n_spikes_class;
end

% set the minimum spike threshold number to be 5
% find the neurons which do not spike at least 5 times in each class, and
% find the unique rows
thres = 5;
[row,col] = find(n_spikes_total<=thres);
row_del = unique(row);

% create a second spike train matrix to store the neurons above the
% threshold
n_spikes_train_ii = n_spikes_train;
n_spikes_test_ii = n_spikes_test;

D_trial_ii = D_trial-length(row_del);
n_spikes_total_ii = zeros(D_trial_ii, n_class);

for i = 1:n_class
    % remove the below-threshold neurons across all classes
    n_spikes_train_ii{1,i}(row_del,:) = [];
    n_spikes_test_ii{1,i}(row_del,:) = [];

    % check that removing these neurons worked
    n_spikes_class_ii = sum(n_spikes_train_ii{1,i},2);
    n_spikes_total_ii(:,i) = n_spikes_class_ii;
    %n_spikes_trial_ii{1,i} = sum(n_spikes_train_ii{1,i},2);
    %n_spikes_total_ii = [n_spikes_total_ii, n_spikes_trial_ii{1,i}];
end

mu_ii = zeros(D_trial_ii, n_class);
S_k_ii = cell(1, n_class);
sigma_ii = zeros(D_trial_ii, D_trial_ii);

for i = 1:n_class
    mu_ii(:,i) = 1/(N_k)*sum(n_spikes_train_ii{1,i},2);
    cov_trial_ii = zeros(D_trial_ii, D_trial_ii);
    for j = 1:n_trial
        cov_trial_ii = cov_trial_ii + (n_spikes_train_ii{1,i}(:,j)-mu_ii(:,i))...
            *(n_spikes_train_ii{1,i}(:,j)-mu_ii(:,i))';
        % sum the (x-mu)*(x-mu)' matrices for each trial 
    end
    S_k_ii{i} = 1/N_k * cov_trial_ii;
    % calculate the S_k for each class and store into cell
    sigma_ii = sigma_ii + N_k/N * S_k_ii{i};
    % calculate sigma (weighted sum of S_k)
end

k_ii = zeros(n_trial, n_class);
errors_ii = zeros(1,n_class);

for n = 1:n_class
    xy = n_spikes_test_ii{n}';
    for j = 1:n_trial
        
        % Warning before removing offending neurons: 
        % Matrix is close to singular or badly scaled. Results may be inaccurate.
        % This happens because there are neurons that do not fire very much across
        % all the trials, causing the covariance matrix to be close to singular.
        
        for i = 1:n_class
            k_ii(j,i) = log(P_Ck)+ mu_ii(:,i)'*inv(S_k_ii{i})*xy(j,:)'-0.5*mu_ii(:,i)'*....
                inv(S_k_ii{i})*mu_ii(:,i) - 0.5*xy(j,:)*inv(S_k_ii{i})*xy(j,:)';
        end
    end
    [m_ii,idx_ii] = max(k_ii, [], 2);
    errors_ii(n)= length(idx_ii(idx_ii~=n));
end
accuracy_ii = 1-sum(errors_ii)/(D_trial_ii*n_class);
