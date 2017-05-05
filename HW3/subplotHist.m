function subplotHist(bar_func, lambda)

x = 0:10;           % 50 bins

subplot(5,3,9)
hold on
histogram(bar_func(1,:),'Normalization','pdf')
plot(x, poisspdf(x,lambda(1)))

subplot(5,3,6)
hold on
bar(bar_func(2,:))
plot(x, poisspdf(x,lambda(2)))

subplot(5,3,2)
hold on
bar(bar_func(3,:))
plot(x, poisspdf(x,lambda(3)))

subplot(5,3,4)
hold on
bar(bar_func(4,:))
plot(x, poisspdf(x,lambda(4)))

subplot(5,3,7)
hold on
bar(bar_func(5,:))
plot(x, poisspdf(x,lambda(5)))

subplot(5,3,10)
hold on
bar(bar_func(6,:))
plot(x, poisspdf(x,lambda(6)))

subplot(5,3,14)
hold on
bar(bar_func(7,:))
plot(x, poisspdf(x,lambda(7)))

subplot(5,3,12)
hold on
bar(bar_func(8,:))
plot(x, poisspdf(x,lambda(8)))
