% EE239.2 HW 3

clc
clear all
close all
%% Problem 3
%% Part a)

t = linspace(0,1,8);
s = (t.^2)*180;

r_0 = 35;
r_max = 60;
s_max = 90;

lambda = r_0 + (r_max - r_0) * cosd(s - s_max);

% lambda_max = r_max;

n = 0;
T = 0;

time = 1;           % spike trains last 1 second
mu = 1./lambda;
T_vec = [];
numTrials = 100;

T_cell = {};         % matrix of 100 spike trains per lambda

for i = 1:length(mu)
    for k = 1:numTrials
        
        while ( T < time )
            dt = exprnd(mu(i));
            n = n + 1;
            T = T + dt;
            
            u = rand(1);                % sample from uniform distribution
            s_Tn = (T^2)*180;
            lambda_Tn = r_0 + (r_max - r_0) * cosd(s_Tn - s_max);
            % calculate threshold value
            if (lambda_Tn/r_max) >= u
                T_vec(n) = T;
            end
            
        end
        
        T_vec = T_vec(:,1:end-1);
        
        % delete last value because it will go over the set time
        
        T_cell{i,k} = T_vec;
        n = 0;       % reset n and t
        T = 0;
    end
    
end



% T_mat_plot = T_mat{1,:};

% for i = 1:5
%
figure(1)
subplotRaster(T_cell)

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

%% Part c)

counts_dist = [];

for i = 1:size(counts_avg,1)
    counts_dist(i,:) = hist(counts_avg(i,:))/length(counts_avg);
end

figure(3)

subplotHist(counts_dist, 0.02*lambda)       % account for 50 bins

% exp(-lambda)*(lambda.^x)/fact(x)

%% Part d)

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
