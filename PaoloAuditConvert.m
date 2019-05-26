% PaoloAuditConvert
clear 
% convert .csv files from Paolo to tag audits
% load csv 
tag = 'bw170813-44';
fs = 10; 

ftoread = strcat('c:/tag/audit/',tag,'prh10_breaths.csv');
% ftoread = strcat('c:/tag/audit/',tag,'_tmp_breaths.csv');
fid = fopen(ftoread);
fgetl(fid); %reads line but does nothing with it
fgetl(fid);
fgetl(fid);
M = textscan(fid, '%f', 'Delimiter',','); % you will need to change the number   of values to match your file %f for numbers and %s for strings.
fclose (fid);

%% remove bad sections or CEE 
bad = 92650+(36*60*60*fs); % end of CEE + 24h 
ii = find(M{1,1}(:,1) > bad); 
M{1,1} = M{1,1}(ii,1); 

%% put in audit structure
R.cue = M{1,1}/fs; % convert, as audit is in samples, not seconds
(R.cue(end,1)-R.cue(1,1))/3600 % check duration 
R.stype = cell(1,length(M{1,1})); 
R.stype(:) = {'breath'}; 

saveaudit(tag,R)