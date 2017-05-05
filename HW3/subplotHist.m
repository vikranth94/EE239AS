function subplotHist(rate, lambda)

x = 0:100;           % 50 bins

subplot(5,3,9)
hold on
histogram(rate(1,:),'Normalization','pdf')
%histogram(bar_func(1,:),'Normalization','pdf')
plot(x, poisspdf(x,lambda(1)))

subplot(5,3,6)
hold on
histogram(rate(2,:),'Normalization','pdf')
plot(x, poisspdf(x,lambda(2)))

subplot(5,3,2)
hold on
histogram(rate(3,:),'Normalization','pdf')
plot(x, poisspdf(x,lambda(3)))

subplot(5,3,4)
hold on
histogram(rate(4,:),'Normalization','pdf')
plot(x, poisspdf(x,lambda(4)))

subplot(5,3,7)
hold on
histogram(rate(5,:),'Normalization','pdf')
plot(x, poisspdf(x,lambda(5)))

subplot(5,3,10)
hold on
histogram(rate(6,:),'Normalization','pdf')
plot(x, poisspdf(x,lambda(6)))

subplot(5,3,14)
hold on
histogram(rate(7,:),'Normalization','pdf')
plot(x, poisspdf(x,lambda(7)))

subplot(5,3,12)
hold on
histogram(rate(8,:),'Normalization','pdf')
plot(x, poisspdf(x,lambda(8)))
