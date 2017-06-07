%% EE239AS HW 7
clc
clear all
close all

load('hw7_data.mat')

%% Problem 2

%% Part A: PCA Plots

Y = Xsim';

Y0 = bsxfun(@minus, Y, mean(Y,2));
S = (1/size(Y0,2))*(Y0*Y0');

[u1,v1] = eigs(S,1);

Xsim_hat = Y0'*u1;

figure(1)
dim_reduce_plot(Y0',Xsim_hat,u1)

%% Part B: PPCA EM Algorithm

D = 1;      % low dimension
[LL, W, s2] = ppca_nsp(Y, D);
figure(2);
plot(LL);