% minke surface detector 
clear, % close all
% load data for file 
filename = 'bb190309-52';
load([gettagpath('PRH') '/' filename ' 10Hzprh']) 

%% correct depth and find surfacings 
surfs = detectresp(p,fs,0.2,1); 

figure(1), clf, hold on 
plot((1:length(p))/fs,-p,'k')
plot(mean(surfs(:,1:2)'),zeros(length(surfs),1),'ro')
%% put in audit structure

R.cue = [surfs(:,1) surfs(:,2)-surfs(:,1)]; % time and duration of detected surfacings
R.stype = cell(1,length(surfs)); 
R.stype(:) = {'dbreath'}; 

%% or load audit
% R = loadaudit(filename); 

%% audit checker
R = point_audit(filename,1,p,fs,R);

%% 
saveaudit(tag,R)