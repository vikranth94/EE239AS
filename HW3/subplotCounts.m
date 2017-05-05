function subplotCounts(Counts,bins)

subplot(5,3,9)
bar(bins(1:end-1)*1000,Counts(1,:),'histc')

subplot(5,3,6)
bar(bins(1:end-1)*1000,Counts(2,:),'histc')

subplot(5,3,2)
bar(bins(1:end-1)*1000,Counts(3,:),'histc')

subplot(5,3,4)
bar(bins(1:end-1)*1000,Counts(4,:),'histc')

subplot(5,3,7)
bar(bins(1:end-1)*1000,Counts(5,:),'histc')

subplot(5,3,10)
bar(bins(1:end-1)*1000,Counts(6,:),'histc')

subplot(5,3,14)
bar(bins(1:end-1)*1000,Counts(7,:),'histc')

subplot(5,3,12)
bar(bins(1:end-1)*1000,Counts(8,:),'histc')


