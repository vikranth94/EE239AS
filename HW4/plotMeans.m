function plotMeans(mu_i)

% plot means of each class (same for each model)
plot(mu_i(1,1), mu_i(2,1),'.r','markersize',50)
plot(mu_i(1,2), mu_i(2,2),'.g','markersize', 50)
plot(mu_i(1,3), mu_i(2,3),'.b','markersize', 50)

end