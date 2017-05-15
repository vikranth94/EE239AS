function plotData(data)

plot(data{1,1}(1,:),data{1,1}(2,:),'xr', data{1,2}(1,:),data{1,2}(2,:), '+g', ...
     data{1,3}(1,:),data{1,3}(2,:), 'ob')

axis([0 20 0 20])
end