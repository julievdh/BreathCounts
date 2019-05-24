function tdiff = waittime(R)

% calculate waiting time between breaths
tdiff = R.cue(1:end-1,1) - R.cue(2:end,1);
tdiff = R.cue(2:end,1) - R.cue(1:end-1,1);

% plot
% figure(1), 
% plot(tdiff)

% figure(2),
% cdfplot(tdiff)

