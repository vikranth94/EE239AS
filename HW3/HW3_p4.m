% EE239.2 HW 3

%% Problem 4

clc
clear 
close all
%% Part a)

ps3_data = importdata('ps3_data.mat');
T_cell = {};

spike_loc = {};
max_trains = zeros(8,5);
max_spikes = zeros(8,5);

for i = 1:size(ps3_data,2)
    for k = 1:size(ps3_data,1)
        spike_cell{i,k} = ps3_data(k,i).spikes;
        T_cell{i,k} = (find(spike_cell{i,k} > 0)-1)/1000;
        % subtract 1 because time indexing starts at t = 0 and MATLAB
        % indexing starts at 1
        % divide by 1000 to convert to seconds
        
        spike_sum(i,k) = sum(spike_cell{i,k});
        
%         if (spike_sum(i,k) >= 5)
%             spike_loc{i} = [i, k];
%         end
        
    end
    [B,I] = sort(spike_sum(i,:),'descend');
    max_trains(i,:) = I(1:5);
    max_spikes(i,:) = B(1:5);

    for j = 1:size(max_trains,2)
        T_cell_plot{i,j} = T_cell{i,max_trains(i,j)};
    end

end

% T_mat_plot = T_mat{1,:};

% for i = 1:5

figure(1)
subplotRaster(T_cell_plot)

%% Part b)

bins = 0:0.020:1;

% h=histogram(T_cell{1,1}, bins)
% counts = h.Values
counts_sum = zeros(8,length(bins));
for i=1:100
    for j=1:8
        counts = histc(T_cell{j,i}, bins);
        counts_sum(j,:) = counts_sum(j,:) + counts;
    end
end

figure(2)
counts_sum = counts_sum(:,1:50);
counts_avg = counts_sum/100;
bar(bins(1:end-1)*1000,counts_avg(1,:),'histc')
xlabel('Time in ms');
ylabel('Firing Rate');

%counts = histc(T_cell{1,1}, bins);
% for j = size(T_cell,1)
%     for k = size(T_cell,2)
%         for i = 1:length(bins)-1
%             numSpikes = length(find(T_cell{j,k} > bins(i) & T_cell{j,k} < bins(i+1)));
%         end
%     end
% end

%% Part c

f_rate = zeros(1,800);
rate = zeros(8,100);
for i=1:8
    for j=1:100
       % f_rate((i-1)*100+j) = length(T_cell{i,j});
       rate(i,j) = length(T_cell{i,j});
    end
end

f_rate = reshape(rate',[1,800]);

s_rep = repmat(s,100,1);
s_rep = reshape(s_rep, 1, numel(s_rep));

scatter(s_rep,f_rate,'x');
f_rate_mean = sum(rate,2)/100;
hold on;
scatter(s,f_rate_mean,'x')
hold on
plot(s,lambda)

%% Part d)

counts_dist = [];

for i = 1:size(counts_avg,1)
    counts_dist(i,:) = hist(counts_avg(i,:))/length(counts_avg);
end

figure(3)

subplotHist(counts_dist, 0.02*lambda)       % account for 50 bins

% exp(-lambda)*(lambda.^x)/fact(x)

%% Part e)

counts_mean = mean(rate,2);
counts_var = var(rate,1,2);
figure(4)
scatter(counts_var, counts_mean)

%% Part f)

ISI = cell(8,1);
ISI_dist = cell(8,1);

for i = 1:length(mu)
    for k = 1:numTrials
        ISI{i} = [ISI{i}, diff(T_cell{i,k})];
    end
    ISI_dist{i} = 100*histc(ISI{i},bins)/sum(histc(ISI{i},bins));
    % normalize by dividing by sum (to get pdf) and multiply by 100 (100
    % trials)
end

figure(5)
subplotISI(ISI_dist,bins,mu)

%% Part g)

for i = 1:length(mu)
    ISI_mean(i) = mean(ISI{i});
    ISI_CV(i) = std(ISI{i})/mean(ISI{i});
end

scatter(ISI_mean,ISI_CV)
ylabel('ISI CV')
xlabel('ISI Mean')