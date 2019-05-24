load('countdata')
figure(9), clf
% for each file
fnum = [9 7 19 37 46 57 81 86 28]; %1:length(files)   
for i = 1:length(fnum); 
    f = fnum(i); 
    % find breaks
    a = diff(files(f).resp');
    % separate based on breaks - longer than 30 s? 
    b = find([a inf]>mean(a));
    % plot minus break start - get IBI within dives
    c = diff([0 b]); % length of the sequences
    d = cumsum(c); % endpoints of the sequences
    
    % plot something here
    figure(9), hold on
    subplot(3,3,i), hold on, title(regexprep(files(f).tag,'_',' '))
    for j = 1:length(d)-1
    h = scatter(files(f).resp(d(j)+1:d(j+1))-files(f).resp(d(j)+1),1:c(j+1),'o','markeredgecolor',files(f).col);
    %h.MarkerFaceAlpha = 0.4; 
    h.MarkerEdgeAlpha = 0.3; 
    end
end

ylabel('Number of breaths'), xlabel('Time (sec)')
print('PostDiveIBI_9spp','-dpng','-r300')

xlim([0 1500]), ylim ([0 50])
linkaxes
print('PostDiveIBI_9spp_axes','-dpng','-r300')