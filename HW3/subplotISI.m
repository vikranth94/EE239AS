function subplotISI(ISI_dist,bins,mu)

subplot(5,3,9)
hold on
bar(bins,ISI_dist{1},'histc')
plot(bins, exppdf(bins,mu(1)),'r')
xlim([0 1])

subplot(5,3,6)
hold on
bar(bins,ISI_dist{2},'histc')
plot(bins, exppdf(bins,mu(2)),'r')
xlim([0 1])

subplot(5,3,2)
hold on
bar(bins,ISI_dist{3},'histc')
plot(bins, exppdf(bins,mu(3)),'r')
xlim([0 1])

subplot(5,3,4)
hold on
bar(bins,ISI_dist{4},'histc')
plot(bins, exppdf(bins,mu(4)),'r')
xlim([0 1])

subplot(5,3,7)
hold on
bar(bins,ISI_dist{5},'histc')
plot(bins, exppdf(bins,mu(5)),'r')
xlim([0 1])

subplot(5,3,10)
hold on
bar(bins,ISI_dist{6},'histc')
plot(bins, exppdf(bins,mu(6)),'r')
xlim([0 1])

subplot(5,3,14)
hold on
bar(bins,ISI_dist{7},'histc')
plot(bins, exppdf(bins,mu(7)),'r')
xlim([0 1])

subplot(5,3,12)
hold on
bar(bins,ISI_dist{8},'histc')
plot(bins, exppdf(bins,mu(8)),'r')
xlim([0 1])

