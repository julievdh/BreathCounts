% gm11_248c figure -- this must exist somewhere 
tag = 'gm11_248c';
loadprh(tag)
R = loadaudit(tag); 

figure(1), clf, hold on 
set(gcf,'position',[68.3333  377.6667  580.0000  240.0000],...
    'paperpositionmode','auto')
t = (1:length(p))/(fs)/3600;   % time in seconds
plot(t,-p,'k','linewidth',1)
xlabel('Time (hours)'), ylabel('Depth (m)')
ylim([-1000 200])

print('gm11_248c_depth','-dpng','-r300')

% find breaths
[~,breath] = findbreathcues(R);
ct = resprate(breath,[0 breath.cue(end,1)]);
tdiff = diff(unique(breath.cue(:,1)));
exth = find(tdiff > 1); % remove anything closer than 1 s 
tdiff = tdiff(exth); resp = breath.cue([exth(1) exth'+1],1);
iRR = 60./tdiff; 
iRR(iRR < 0.5) = NaN;

% plot(resp(1:end-1,1)/3600,iRR*20,'color',[0.5 0.5 0.5])

plot(ct(2,:)/3600,ct(1,:)*20,'color',[0.5 0.5 0.5])

print('gm11_248c_depth_resp','-dpng','-r300') 

figure(3), clf, hold on
histogram(iRR,'facecolor',[0.5 0.5 0.5])
axis ij, axis off 
print('gm11_248c_depth_resp_hist','-dpng','-r300')

%% zoom 
figure(2), clf, hold on 
set(gcf,'position',[68.3333  377.6667  852.0000  240.0000],...
    'paperpositionmode','auto')
t = (1:length(p))/(fs)/3600;   % time in seconds
plot(t,-p,'k','linewidth',1)
xlabel('Time (hours)'), ylabel('Depth (m)')
ylim([-1000 200])

plot(resp(1:end-1,1)/3600,iRR,'.-','color',[0.5 0.5 0.5])
xlim([10.8 12.2]), ylim([-50 20])

print('gm11_248c_depth_resp_zoom','-dpng','-r300')

wthn = find(iswithin(resp/3600,[10.4 12.2]));
figure(3), hold on
histogram(iRR(wthn),'facecolor',[0.15 0.15 0.15])
print('gm11_248c_depth_resp_hist_zoom','-dpng','-r300')
