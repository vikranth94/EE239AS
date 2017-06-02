function spike_count = binFunc(raster_mat, dt)

n_neurons = size(raster_mat,1);
time = size(raster_mat,2);
n_bins = ceil(time/dt);
n_trials = size(raster_mat,3);

spike_count = zeros(n_neurons,n_bins,n_trials);
% size of spike_count
% number of trials x time/dt
i = 1;          % index for time interval
j = 1;          % index for spike_count
while j <= floor(time/dt)
    spike_count(:,j,:) = sum(raster_mat(:, i:i+dt-1,:),2);
    i = i+dt;
    j = j+1;
end
% all multiples of dt get saved
if i-1 < time
    % because matlab starts indexing at 1
    spike_count(:,j,:) = sum(raster_mat(:, i-dt:end,:),2);
end
% save last (not multiple of dt)

end