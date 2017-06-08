%% EE239AS HW 7
clc
clear all
close all

load('hw7_data.mat')

%% Problem 2

%% Part A: PCA Visualization

Y = Xsim';

Y0 = bsxfun(@minus, Y, mean(Y,2));
S = (1/size(Y0,2))*(Y0*Y0');

[u1,v1] = eigs(S,1);

Xsim_hat = Y0'*u1;

figure(1)
dim_reduce_plot(Y0',Xsim_hat,u1)
title('PCA Visualization')
xlabel('Dimension 1')
ylabel('Dimension 2')

%% Part B: PPCA EM Algorithm

D = 1;      % low dimensional space
[LL, W, s2] = ppca_nsp(Y, D);
figure(2);
plot(LL)
title('Part B: Log Data Likehood (PPCA)')
xlabel('Iteration')
ylabel('Log Data Likelihood')

%% Part C: PPCA Covariance

[N, K] = size(Xsim');

cov_sample  = cov(Xsim, 1);
fprintf('\nSample Covariance:\n')
disp(cov_sample)

cov_PPCA = (W*W' + s2*eye(N));
fprintf('\nPPCA Covariance:\n')
disp(cov_PPCA)

% The sample covariance and PPCA covariance are very similar.

%% Part D: PPCA Visualization

mu = mean(Y, 2);
mu_mat = repmat(mu,1,K);

Es = W'*inv(W*W' + s2*eye(N))*(Y-mu_mat);

figure(3)
dim_reduce_plot(Xsim,Es',W)
title('PPCA Visualization')
xlabel('Dimension 1')
ylabel('Dimension 2')
        
%% Part E: FA EM Algorithm

D = 1;      % low dimensional space
[LL_FA, W_FA, psi] = fa_nsp(Y, D);
figure(4);
plot(LL_FA)
title('Part D: Log Data Likehood (FA)')
xlabel('Iteration')
ylabel('Log Data Likelihood')

% The sample covariance and FA covariance are very similar.

%% Part F: FA Covariance

[N, K] = size(Xsim');

cov_sample  = cov(Xsim, 1);
fprintf('\nSample Covariance:\n')
disp(cov_sample)

cov_FA = (W_FA*W_FA' + psi);
fprintf('\nFA Covariance:\n')
disp(cov_FA)

% The sample covariance and FA covariance are very similar.

%% Part G: FA Visualization

mu = mean(Y, 2);
mu_mat = repmat(mu,1,K);

Es_FA = W_FA'*inv(W_FA*W_FA' + psi)*(Y-mu_mat);

figure(5)
dim_reduce_plot(Xsim,Es_FA',W_FA)
title('FA Visualization')
xlabel('Dimension 1')
ylabel('Dimension 2')