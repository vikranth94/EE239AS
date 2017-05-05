r_0 = 35;
r_max = 60;
s_max = 90;

t = linspace(0,1,100);
lambda = r_0 + (r_max - r_0) * cosd(t.^2*180 - s_max);

plot(t, lambda)