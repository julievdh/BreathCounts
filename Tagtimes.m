% Plot wake and sleep times of different individuals on clock face with
% overlapping patches

% load tagtimes (example from gm data Tarifa from Frants)
% tag on time hour of day, minute of day, duration of recording
load('tagtimes')

for i = 1:length(tagtimes)
    tagtimes(i,4) = datenum([0 0 0 tagtimes(i,1) tagtimes(i,2) 0]);
    
    % add tag duration to tag on time to get end time
    warning off
    tagtimes(i,5) = addtodate(tagtimes(i,4),tagtimes(i,3)*60,'min');
    
    % get end time hour and minute
    tagoff(i,:) = datevec(tagtimes(i,5));
    tagtimes(i,6:7) = tagoff(i,4:5);
    % IF NEXT DAY, FLAG AND ADD 24 H
    if tagoff(i,3) == 1;
        tagtimes(i,6) = tagtimes(i,6)+24;
    end
end

rad = 1.5; % radius of plot

figure(1); clf; % hold on

for i = 1:length(tagtimes)
    timeClockPer(tagtimes(i,1:2),tagtimes(i,6:7),rad,0.05)
end

% add clock labels
text(0,rad*1.2,'24','FontSize',14)
text(0,-rad*1.2,'12','FontSize',14)
text(rad*1.2,0,'6','FontSize',14)
text(-rad*1.2,0,'18','FontSize',14)

axis off
title('Times of Pilot Whale Tag Deployments','FontSize',14)
