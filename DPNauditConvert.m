% tagcue conversion for DPN audits from past
% cue = [];
tag = 'eg05_224a';
% cue is matrix of chip number and seconds in chip, copy and paste from DPN audit .xls

N = makecuetab(tag); % make cuetab
for i = 1:max(cue(:,1))
    ii = find(cue(:,1) == i);
    ctime(ii,1) = N.N(i,2)+cue(ii,2);
end

% ctime then becomes the cue time used for audit file

%% check audit?
R = loadaudit(tag);
loadprh(tag)
figure(18), clf
t = (1:length(p))/fs;   % time in seconds
plot(t,-p)
hold on
for i = 1:length(R.cue)
    plot(R.cue(i,1),0,'k*')
end
[~,NB] = findaudit(R,'NB'); % all Noisy Breaths
for i = 1:length(NB.cue(:,1))
    plot(NB.cue(i,1),0,'g*')
end
[~,Q] = findaudit(R,'B?'); % all Questionable Breaths
for i = 1:length(Q.cue)
    plot(Q.cue(i,1),0,'r*')
end

%% compute jerk and find peaks
j = njerk(Aw,fs);
plot([1 t(end)],[quantile(j,0.9) quantile(j,0.9)])
plot(t,j)
% find peaks in jerk signal
[pks,locs] = findpeaks(j,t(2:end),'minpeakheight',quantile(j,0.9),'minpeakdist',10);
plot(locs,pks,'o')

% find peaks within surfacings
T = finddives(p,fs,5,1,0);
S = nan(length(T),2);
S(1,1) = 1; S(1,2) = T(1,1); % assume tag begins at surface before first dive
S(2:length(T),1) = T(1:end-1,2);
S(2:length(T),2) = T(2:end,1); % surfacings

for i = 1:length(S) % all surfacings
    ii = find(iswithin(locs,S(i,:)));
    if i == 1
    surfjerk = ii;
    else surfjerk = horzcat(ii,surfjerk); 
    end
end

plot(locs(surfjerk),pks(surfjerk)','k.','markersize',15)
    
