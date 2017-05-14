% EE239AS Homework 4

clc
clear
close all

%% Problem 3: Simulated Neural Data

ps4_data = importdata('ps4_realdata.mat');

% 20x3 struct
% rows = data point
% columns = class
train_data = ps4_data.train_trial;


sum(train_data(1,1).spikes(1,:))

