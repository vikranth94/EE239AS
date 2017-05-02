% EE239.2 HW 3

%% Part a)

s = [0 45 90 135 180 225 270 315];

r_0 = 35;
r_max = 60;
s_max = pi/2;

lambda = r_0 + (r_max - r_0) * cosd(s - s_max);

n = 0;
t = 0;

time = 1;           % spike trains last 1 second
mu = 1./lambda;
T_vec = [];
numTrials = 100;

T_mat = {};         % matrix of 100 spike trains per lambda

for k = 1:numTrials
    for i = 1:length(mu)

      while ( t < time )
        dt = exprnd(mu(i));
        n = n + 1;
        t = t + dt;
        T_vec(n) = t;
      end
        T_vec = T_vec(:,1:end-1);           
        % delete last value because it will go over the set time 
        
        T_mat{i,k} = T_vec;
        n = 0;       % reset n and t 
        t = 0;
    end 
end

T_mat_plot = T_mat{1,:};

for i = 1:5
 
T_mat_plot = T_mat
plotRaster(T_mat_plot)

end