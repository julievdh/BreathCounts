% for each whale, plot IBI after dive
%plot(files(i).resp(2:end),files(i).tdiff,'o')

for i = 19:length(files)
% set a threshold for dive duration - 95% of tdiff 
th = quantile(files(i).tdiff,0.95);
% find above threshold
ii = find(files(i).tdiff > th)+1;
% plot each *interval*
figure(8), clf, hold on
for k = 1:length(ii)-1;
plot(files(i).resp(ii(k):ii(k+1))-files(i).resp(ii(k)),files(i).tdiff(ii(k):ii(k+1)),'o')
files(i).chunks{k}(:,1) = files(i).resp(ii(k):ii(k+1))-files(i).resp(ii(k));
files(i).chunks{k}(:,2) = files(i).tdiff(ii(k):ii(k+1));
end

xlabel('Time (sec)'), ylabel('IBI (sec)')
pause
end
