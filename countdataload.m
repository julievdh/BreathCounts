% files for BreathCounts

th = 1; % second threshold between breaths  

% tt14_126a
files(1).tag = 'tt14_126a';
files(1).path = 'C:\tag\audit';
files(1).lnth = 274; % length in cm
files(1).wt = 256; % weight in kg

% tt14_126b
files(2).tag = 'tt14_126b';
files(2).path = 'C:\tag\audit';
files(2).lnth = 285;
files(2).wt = 291; % weight in kg

% tt14_127a
files(75).tag = 'tt14_127a';
files(75).path = 'C:\tag\audit';
files(75).lnth = 234;
files(75).wt = 169; % weight in kg

% % tt14_127b
% files(4).tag = 'tt14_127b';
% files(4).path = strcat(res,':\tt14\tt14_127b\');

% tt14_128a
files(3).tag = 'tt14_128a';
files(3).path = 'C:\tag\audit';
files(3).lnth = 242;
files(3).wt = 167;

% tt14_129a
files(4).tag = 'tt14_129a';
files(4).path = 'C:\tag\audit';
files(4).lnth = 252;
files(4).wt = 172;

% tt14_129c
files(5).tag = 'tt14_129c';
files(5).path = 'C:\tag\audit';
files(5).lnth = 262;
files(5).wt = 291;

% tt14_129d
files(6).tag = 'tt14_129d';
files(6).path = 'C:\tag\audit';
files(6).lnth = 281;
files(6).wt = 276;

% tt14_125a
files(7).tag = 'tt14_125a';
files(7).path = 'C:\tag\audit';
files(7).lnth = 269;
files(7).wt = 248;

% tt14_125b
files(8).tag = 'tt14_125b';
files(8).path = 'C:\tag\audit';
files(8).lnth = 273;
files(8).wt = 227;

% tt15_134a
fl = 99; 
files(fl).tag = 'tt14_125b';
files(fl).path = 'C:\tag\audit';
files(fl).lnth = 236; % cm
files(fl).wt = 142; % kg

%% pilot whales from Spain
files(17).tag = 'gm13_220b';
files(17).path = 'C:\tag\audit';
files(17).lnth = NaN; % cm estimate
files(17).wt = 1500; % kg estimate, Jefferson et al. 1993
files(17).sex = 'M';

files(56).tag = 'gm15_266a';
files(56).path = 'C:\tag\audit';
files(56).lnth = NaN; % large male
files(56).wt = 1500; % kg estimate, Jefferson et al. 1993
files(56).sex = 'M';

files(57).tag = 'gm15_265a';
files(57).path = 'C:\tag\audit';
files(57).lnth = NaN; % large male, 
files(57).wt = 1500; % kg estimate, Jefferson et al. 1993
files(57).sex = 'M';

fl = 58; 
files(fl).tag = 'gm10_231b';
files(fl).path = 'C:\tag\audit';
files(fl).lnth = NaN; % large male
files(fl).wt = 1500; % kg estimate, Jefferson et al. 1993
files(fl).sex = 'M';

fl = 59; 
files(fl).tag = 'gm11_248c';
files(fl).path = 'C:\tag\audit';
files(fl).lnth = NaN; % adult female
files(fl).wt = 1200; % kg estimate, Jefferson et al. 1993 (2000 for males)
files(fl).sex = 'F';

fl = 92; 
files(fl).tag = 'gm10_238b';
files(fl).path = 'C:\tag\audit';
files(fl).lnth = NaN; % adult female
files(fl).wt = 1500; % kg estimate, Jefferson et al. 1993 (2000 for males)
files(fl).sex = 'M';

fl = 93; 
files(fl).tag = 'gm10_248c';
files(fl).path = 'C:\tag\audit';
files(fl).lnth = NaN; % adult female
files(fl).wt = 1500; % kg estimate, Jefferson et al. 1993 (2000 for males)
files(fl).sex = 'M';

fl = 94; 
files(fl).tag = 'gm11_235d';
files(fl).path = 'C:\tag\audit';
files(fl).lnth = NaN; % adult female
files(fl).wt = 1300; % kg estimate, Jefferson et al. 1993 (2000 for males)
files(fl).sex = NaN;

fl = 95; 
files(fl).tag = 'gm11_241c';
files(fl).path = 'C:\tag\audit';
files(fl).lnth = NaN; % adult female
files(fl).wt = 1500; % kg estimate, Jefferson et al. 1993 (2000 for males)
files(fl).sex = 'M';

fl = 96; 
files(fl).tag = 'gm12_224d';
files(fl).path = 'C:\tag\audit';
files(fl).lnth = NaN; % adult female
files(fl).wt = 1500; % kg estimate, Jefferson et al. 1993 (2000 for males)
files(fl).sex = 'M';

fl = 97; 
files(fl).tag = 'gm11_250b';
files(fl).path = 'C:\tag\audit';
files(fl).lnth = NaN; % adult female
files(fl).wt = 1500; % kg estimate, Jefferson et al. 1993 (2000 for males)
files(fl).sex = 'M';

fl = 98; 
files(fl).tag = 'gm12_207a';
files(fl).path = 'C:\tag\audit';
files(fl).lnth = NaN; % adult female
files(fl).wt = 1500; % kg estimate, Jefferson et al. 1993 (2000 for males)
files(fl).sex = 'M';

%% harbour porpoises
files(9).tag = 'hp12_272a';
files(9).path = 'C:\tag\audit';
files(9).lnth = 122;
files(9).sex = 'F'; 
files(9).wt = 0.000083*files(9).lnth^2.67; % kg estimate

files(10).tag = 'hp12_293a';
files(10).path = 'C:\tag\audit';
files(10).lnth = 163;
files(10).sex = 'F'; 
files(10).wt = 0.000083*files(10).lnth^2.67; % kg estimate

files(11).tag = 'hp13_102a';
files(11).path = 'C:\tag\audit';
files(11).lnth = 114;
files(11).sex = 'M'; 
files(11).wt = 0.00063*files(11).lnth^2.25; % kg estimate Rojano Donate 2018

files(12).tag = 'hp13_145a';
files(12).path = 'C:\tag\audit';
files(12).lnth = 135;
files(12).sex = 'F'; 
files(12).wt = 0.000083*files(12).lnth^2.67; % kg estimate

files(13).tag = 'hp13_170a';
files(13).path = 'C:\tag\audit';
files(13).lnth = 123;
files(13).sex = 'M'; 
files(13).wt = 0.00063*files(13).lnth^2.25; % kg estimate

files(14).tag = 'hp14_226b';
files(14).path = 'C:\tag\audit';
files(14).lnth = 126;
files(14).sex = 'M';
files(14).wt = 0.00063*files(14).lnth^2.25; % kg estimate

files(15).tag = 'hp15_117a';
files(15).path = 'C:\tag\audit';
files(15).lnth = 170;
files(15).sex = 'F'; 
files(15).wt = 0.000083*files(15).lnth^2.67; % kg estimate

files(16).tag = 'hp16_264a';
files(16).path = 'C:\tag\audit';
files(16).lnth = 163;
files(16).sex = 'F'; 
files(16).wt = 0.000083*files(16).lnth^2.67; % kg estimate

fl = 61; 
files(fl).tag = 'hp14_305a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = 'M';
files(fl).lnth = 133;
files(fl).wt = 0.00063*files(fl).lnth^2.25; % kg estimate

fl = 62; 
files(fl).tag = 'hp15_096a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = 'M';
files(fl).lnth = 128;
files(fl).wt = 0.00063*files(fl).lnth^2.25; % kg estimate

fl = 63; 
files(fl).tag = 'hp15_160a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = 'F';
files(fl).lnth = 145;
files(fl).wt = 0.000083*files(fl).lnth^2.67; % kg estimate

fl = 64; 
files(fl).tag = 'hp15_218a';
files(fl).path = 'C:\tag\audit';
files(fl).lnth = 156;
files(fl).sex = 'F';
files(fl).wt = 0.000083*files(fl).lnth^2.67; % kg estimate

fl = 65; 
files(fl).tag = 'hp15_243a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = 'F';
files(fl).lnth = 94;
files(fl).wt = 0.000083*files(fl).lnth^2.67; % kg estimate

fl = 66; 
files(fl).tag = 'hp15_267a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = 'M';
files(fl).lnth = 121;
files(fl).wt = 0.00063*files(fl).lnth^2.25; % kg estimate

fl = 67; 
files(fl).tag = 'hp16_316a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = 'M';
files(fl).lnth = 113;
files(fl).wt = 0.00063*files(fl).lnth^2.25; % kg estimate

fl = 68; 
files(fl).tag = 'hp17_115a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = 'M';
files(fl).lnth = 134;
files(fl).wt = 0.00063*files(fl).lnth^2.25; % kg estimate

fl = 69; 
files(fl).tag = 'hp17_267a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = 'M';
files(fl).lnth = 121;
files(fl).wt = 0.00063*files(fl).lnth^2.25; % kg estimate

%% Short-finned pilot whales
fl = 18; 
files(fl).tag = 'pw03_308c';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\pwrespdata\';
addpath(files(fl).path)
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
eval(strcat([files(fl).tag(1:2) files(fl).tag(6:9)],'resp;')) % run the resp code
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP RESPAURAL 

fl = 19;
files(fl).tag = 'pw04_299c';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\pwrespdata\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
eval(strcat([files(fl).tag(1:2) files(fl).tag(6:9)],'resp;')) % run the resp code
files(fl).resp = sort(unique(resp)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th));
clear RESP resp


fl = 20;
files(fl).tag = 'pw03_306b';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\pwrespdata\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
eval(strcat([files(fl).tag(1:2) files(fl).tag(6:9)],'resp;')) % run the resp code
files(fl).resp = vertcat(RESP,blowsafterlastdive); 
files(fl).tdiff = diff(sort(vertcat(RESP,blowsafterlastdive))); clear RESP RESPAURAL blowsafterlastdive


fl = 21;
files(fl).tag = 'pw03_306d';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\pwrespdata\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
eval(strcat([files(fl).tag(1:2) files(fl).tag(6:9)],'resp;')) % run the resp code
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP 

fl = 22;
files(fl).tag = 'pw03_307a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\pwrespdata\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
eval(strcat([files(fl).tag(1:2) files(fl).tag(6:9)],'resp;')) % run the resp code
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP RESPAURAL resp

fl = 23;
files(fl).tag = 'pw03_307b';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\pwrespdata\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
eval(strcat([files(fl).tag(1:2) files(fl).tag(6:9)],'resp;')) % run the resp code
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP RESPAURAL resp

fl = 24; % lots of 'unclear' breaths 
files(fl).tag = 'pw03_307c';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\pwrespdata\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
eval(strcat([files(fl).tag(1:2) files(fl).tag(6:9)],'resp;')) % run the resp code
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP RESPAURAL resp

fl = 25;
files(fl).tag = 'pw03_308d';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\pwrespdata\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
eval(strcat([files(fl).tag(1:2) files(fl).tag(6:9)],'resp')); % run the resp code
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP RESPAURAL resp

fl = 26;
files(fl).tag = 'pw03_309b';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\pwrespdata\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
eval(strcat([files(fl).tag(1:2) files(fl).tag(6:9)],'resp;')) % run the resp code
files(fl).resp = vertcat(RESP,Respafterlastdive); 
files(fl).tdiff = diff(sort(vertcat(RESP,Respafterlastdive))); clear RESP Respafterlastdive

fl = 27;
files(fl).tag = 'pw03_309c';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\pwrespdata\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
eval(strcat([files(fl).tag(1:2) files(fl).tag(6:9)],'resp;')) % run the resp code
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP RESPAURAL resp

fl = 69; 
files(fl).tag = 'pw04_295b';
files(fl).path = 'C:\tag\audit';
files(fl).sex = 'F';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993

fl = 70; 
files(fl).tag = 'pw04_297f';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993

fl = 71; 
files(fl).tag = 'pw04_297h';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993

fl = 72; 
files(fl).tag = 'pw08_108a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = 'F';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993

fl = 73; 
files(fl).tag = 'pw08_110a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993

fl = 74; 
files(fl).tag = 'pw08_113c';
files(fl).path = 'C:\tag\audit';
files(fl).sex = 'F';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993


%% sperm whales

fl = 28;
files(fl).tag = 'sw05_196a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\sperm whale\sw05196a\';
files(fl).lnth = 16.3; % from IPI, Teloni et al. 2008
files(fl).wt = 51000; % kg estimate, from IPI, Teloni et al. 2008
load([files(fl).path 'sw05196a_resp'])
files(fl).resp = t; 
files(fl).tdiff = diff(t);

fl = 29;
files(fl).tag = 'sw05_199a';
files(fl).path = 'C:\tag\audit';
files(fl).lnth = 14.3; % from IPI, Teloni et al. 2008
files(fl).wt = 36000; % kg estimate, Teloni et al. 2008
sw05_199cRESP
files(fl).resp = cues(:,1); 
files(fl).tdiff = diff(cues(:,1));

fl = 60;
files(fl).tag = 'sw17_204a';
files(fl).path = 'C:\tag\audit';
files(fl).lnth = IPI2BL(677,192012); % from IPI CHECK
files(fl).wt = (1.25*0.0196*(files(fl).lnth)^2.74)*10^3; % kg estimate, Lockyer 1976, Rice 1989 relationship cited in Teloni
R = loadaudit(files(fl).tag);
[cues,breaths] = findbreathcues(R); 
files(fl).resp = cues(:,1); 
files(fl).tdiff = diff(cues(:,1));

% fl = 126;
% files(fl).tag = 'sw17_196a';
% files(fl).path = 'C:\tag\audit';
% files(fl).lnth = 9.3; % from IPI CHECK
% files(fl).wt = 30000; % kg estimate, GUESS RIGHT NOW
% R = loadaudit(files(fl).tag);
% [cues,breaths] = findbreathcues(R); 
% files(fl).resp = cues(:,1); 
% files(fl).tdiff = diff(cues(:,1));



%% ziphius
fl = 30;
files(fl).tag = 'zc03_263a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\resp_ziphius\Zc_OK\';
files(fl).lnth = NaN;
files(fl).wt = 2000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RespDepthOK')])
files(fl).resp = sort(RESP); 
files(fl).tdiff = diff(sort(RESP)); clear RESP

fl = 31;
files(fl).tag = 'zc04_160a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\resp_ziphius\Zc_OK\';
files(fl).lnth = NaN;
files(fl).wt = 2000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RespDepthOK')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP

fl = 32;
files(fl).tag = 'zc04_161a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\resp_ziphius\Zc_OK\';
files(fl).lnth = NaN;
files(fl).wt = 2000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RespDepthOK')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP

fl = 33;
files(fl).tag = 'zc04_161b';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\resp_ziphius\Zc_OK\';
files(fl).lnth = NaN;
files(fl).wt = 2000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RespDepthOK')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP

fl = 34;
files(fl).tag = 'zc04_175a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\resp_ziphius\Zc_OK\';
files(fl).lnth = NaN;
files(fl).wt = 2000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RespDepthOK')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP

fl = 35;
files(fl).tag = 'zc04_179a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\resp_ziphius\Zc_OK\';
files(fl).lnth = NaN;
files(fl).wt = 2000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RespDepthOK')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP

fl = 36;
files(fl).tag = 'zc05_167a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\resp_ziphius\Zc_OK\';
files(fl).lnth = NaN;
files(fl).wt = 2000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RespDepthOK')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP

fl = 37;
files(fl).tag = 'zc05_170a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\resp_ziphius\Zc_OK\';
files(fl).lnth = NaN;
files(fl).wt = 2000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RespDepthOK')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP

fl = 38;
files(fl).tag = 'zc06_204a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\resp_ziphius\Zc_OK\';
files(fl).lnth = NaN;
files(fl).wt = 2000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RespDepthOK')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP

fl = 39;
files(fl).tag = 'zc06_205a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\resp_ziphius\Zc_OK\';
files(fl).lnth = NaN;
files(fl).wt = 2000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RespDepthOK')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP

fl = 40;
files(fl).tag = 'zc15_212a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\resp_ziphius\Zc_OK\';
files(fl).lnth = NaN;
files(fl).wt = 2000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RespDepthOK')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP

fl = 41;
files(fl).tag = 'zc15_212b';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\resp_ziphius\Zc_OK\';
files(fl).lnth = NaN;
files(fl).wt = 2000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RespDepthOK')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th));clear RESP

%% mesoplodon
fl = 42;
files(fl).tag = 'md03_284a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\mesoplodon\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RESP')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP 

fl = 43;
files(fl).tag = 'md03_298a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\mesoplodon\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RESP')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP 

fl = 44;
files(fl).tag = 'md04_287a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\mesoplodon\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RESP')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP 

fl = 45;
files(fl).tag = 'md05_277a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\mesoplodon\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RESP')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP 

fl = 46;
files(fl).tag = 'md05_285a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\mesoplodon\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RESP')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP 

fl = 47;
files(fl).tag = 'md05_294a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\mesoplodon\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RESP')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP 

fl = 48;
files(fl).tag = 'md05_294b';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\mesoplodon\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RESP')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP 

fl = 49;
files(fl).tag = 'md08_137a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\mesoplodon\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RESP')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP 

fl = 50;
files(fl).tag = 'md08_142a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\mesoplodon\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RESP')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP 

fl = 51;
files(fl).tag = 'md08_148a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\mesoplodon\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RESP')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP 

fl = 52;
files(fl).tag = 'md08_289a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\mesoplodon\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RESP')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP 

fl = 53;
files(fl).tag = 'md10_146a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\mesoplodon\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RESP')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP 

fl = 54;
files(fl).tag = 'md10_163a';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\mesoplodon\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
load([files(fl).path strcat(files(fl).tag,'RESP')])
files(fl).resp = sort(unique(RESP)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th)); clear RESP 

fl = 55;
files(fl).tag = 'md296';
files(fl).path = 'O:\ST_Bio-Kaskelot\van der Hoop Other BreathCounts\resps for julie\mesoplodon\md298a\';
files(fl).lnth = NaN;
files(fl).wt = 1000; % kg estimate, Jefferson et al. 1993
load([files(fl).path files(fl).tag])
files(fl).resp = sort(unique(resp)); 
files(fl).tdiff = diff(files(fl).resp); 
files(fl).tdiff = files(fl).tdiff(find(files(fl).tdiff > th))'; clear resp 

%% Hyperoodon - Patrick Miller

fl = 76; 
files(fl).tag = 'ha14_165a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 750; % cm 
files(fl).wt = 0.0000131*files(fl).lnth^3.07; % kg estimate, Jefferson et al. 1993

fl = 77; 
files(fl).tag = 'ha14_166a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 750; % cm 
files(fl).wt = 0.0000131*files(fl).lnth^3.07; % kg estimate, Jefferson et al. 1993

fl = 78; 
files(fl).tag = 'ha14_174b';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 750; % cm 
files(fl).wt = 0.0000131*files(fl).lnth^3.07; % kg estimate, Jefferson et al. 1993

fl = 79; 
files(fl).tag = 'ha14_175a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 750; % cm 
files(fl).wt = 0.0000131*files(fl).lnth^3.07; % kg estimate, Jefferson et al. 1993

fl = 80; 
files(fl).tag = 'ha15_173a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 750; % cm 
files(fl).wt = 0.0000131*files(fl).lnth^3.07; % kg estimate, Jefferson et al. 1993

fl = 81; 
files(fl).tag = 'ha15_174a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 750; % cm 
files(fl).wt = 0.0000131*files(fl).lnth^3.07; % kg estimate, Jefferson et al. 1993

fl = 82; 
files(fl).tag = 'ha14_174b';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 750; % cm 
files(fl).wt = 0.0000131*files(fl).lnth^3.07; % kg estimate, Jefferson et al. 1993

fl = 83; 
files(fl).tag = 'ha16_173a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 750; % cm 
files(fl).wt = 0.0000131*files(fl).lnth^3.07; % kg estimate, Jefferson et al. 1993

%% beluga whales - Aran Mooney, Manolo Castellote

fl = 84; 
files(fl).tag = 'dl16_133a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = 'F';
files(fl).age = 6.3;
files(fl).lnth = 259; % cm 
files(fl).wt = 10^(-3.95)*(259)^(1.08)*(165)^1.71; % kg estimate, Doidge et al. 1990

fl = 85; 
files(fl).tag = 'dl16_134a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = 'F';
files(fl).age = 9.2;
files(fl).lnth = 284; % cm 
files(fl).wt = 10^(-3.95)*(284)^(1.08)*(170)^1.71; % kg estimate, Jefferson et al. 1993

fl = 86; 
files(fl).tag = 'dl16_135a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = 'M';
files(fl).age = 5;
files(fl).lnth = 251; % cm 
files(fl).wt = 10^(-4.33)*(251)^(2.46)*(159)^0.36; % kg estimate, Jefferson et al. 1993

fl = 87; 
files(fl).tag = 'dl16_136a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = 'M';
files(fl).age = 19;
files(fl).lnth = 384; % cm 
files(fl).wt = 10^(-4.33)*(384)^(2.46)*(211)^0.36; % kg estimate, Jefferson et al. 1993

fl = 88; 
files(fl).tag = 'dl16_138a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = 'M';
files(fl).age = 15;
files(fl).lnth = 363; % cm 
files(fl).wt = 10^(-4.33)*(363)^(2.46)*(213)^0.36; % kg estimate, Jefferson et al. 1993

fl = 89; 
files(fl).tag = 'dl14_237a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = 'M';
files(fl).lnth = 366; % cm 
files(fl).wt = 10^(-4.33)*(363)^(2.46)*(213)^0.36; % kg estimate, Jefferson et al. 1993

fl = 90; 
files(fl).tag = 'dl14_238a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = 'M';
files(fl).age = 15;
files(fl).lnth = 351; % cm 
files(fl).wt = 10^(-4.33)*(351)^(2.46)*(226)^0.36; % kg estimate, Jefferson et al. 1993

fl = 91; 
files(fl).tag = 'dl14_240a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = 'M';
files(fl).age = 18;
files(fl).lnth = 386; % cm 
files(fl).wt = 10^(-4.33)*(386)^(2.46)*(241)^0.36; % kg estimate, Jefferson et al. 1993

%% humpback whales
fl = 100; 
files(fl).tag = 'mn07_203a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = NaN;
files(fl).wt = 30000; % kg estimate, Jefferson et al. 1993

fl = 101; 
files(fl).tag = 'mn12_178a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = NaN;
files(fl).wt = 30000; % kg estimate, Jefferson et al. 1993

fl = 102; 
files(fl).tag = 'mn12_184a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = NaN;
files(fl).wt = 30000; % kg estimate, Jefferson et al. 1993

fl = 103; 
files(fl).tag = 'mn12_185a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = NaN;
files(fl).wt = 30000; % kg estimate, Jefferson et al. 1993

%% Minke whales ´- Paolo, Jeremy, Ari
fl = 104; 
files(fl).tag = 'ba13_265a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 800;
files(fl).wt = 9000; % kg estimate, Jefferson et al. 1993

fl = 105; 
files(fl).tag = 'ba14_211a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 800;
files(fl).wt = 9000; % kg estimate, Jefferson et al. 1993

fl = 106; 
files(fl).tag = 'bb-130213-B009'; % BB - ANTARCTIC MINKE
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 800;
files(fl).wt = 9000; % kg estimate, Jefferson et al. 1993

fl = 107; 
files(fl).tag = 'bb-130215-B008'; % BB - ANTARCTIC MINKE
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 800;
files(fl).wt = 9000; % kg estimate, Jefferson et al. 1993

%% blue whales - Paolo, Jeremy
fl = 108; 
files(fl).tag = 'bw10_229b';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 2500;
files(fl).wt = 92000; % kg estimate, Jefferson et al. 1993

fl = 109; 
files(fl).tag = 'bw10_241a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 2500;
files(fl).wt = 92000; % kg estimate, Jefferson et al. 1993

fl = 110; 
files(fl).tag = 'bw14_212a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 2500;
files(fl).wt = 92000; % kg estimate, Jefferson et al. 1993

fl = 111; 
files(fl).tag = 'bw150823-6';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 2500;
files(fl).wt = 92000; % kg estimate, Jefferson et al. 1993

fl = 112; 
files(fl).tag = 'bw160727-10';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 2500;
files(fl).wt = 92000; % kg estimate, Jefferson et al. 1993

fl = 113; 
files(fl).tag = 'bw170619-36';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 2500;
files(fl).wt = 92000; % kg estimate, Jefferson et al. 1993

fl = 114; 
files(fl).tag = 'bw170813-44';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 2500;
files(fl).wt = 92000; % kg estimate, Jefferson et al. 1993

fl = 115; 
files(fl).tag = 'bw170814-40';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 2500;
files(fl).wt = 92000; % kg estimate, Jefferson et al. 1993

fl = 116; 
files(fl).tag = 'bw170814-51';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 2500;
files(fl).wt = 92000; % kg estimate, Jefferson et al. 1993

fl = 117; 
files(fl).tag = 'bw20160717';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 2500;
files(fl).wt = 92000; % kg estimate, Jefferson et al. 1993

%% fin whales - Paolo, Jeremy
fl = 118; 
files(fl).tag = 'bp160912';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 2000;
files(fl).wt = 52000; % kg estimate, Jefferson et al. 1993

fl = 119; 
files(fl).tag = 'bp13_258a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 2000;
files(fl).wt = 52000; % kg estimate, Jefferson et al. 1993

fl = 120; 
files(fl).tag = 'bp13_258b';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 2000;
files(fl).wt = 52000; % kg estimate, Jefferson et al. 1993

fl = 121; 
files(fl).tag = 'bp13_258c';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 2000;
files(fl).wt = 52000; % kg estimate, Jefferson et al. 1993

fl = 122; 
files(fl).tag = 'bp14_259a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 2000;
files(fl).wt = 52000; % kg estimate, Jefferson et al. 1993

fl = 123; 
files(fl).tag = 'bp15_079a';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 2000;
files(fl).wt = 52000; % kg estimate, Jefferson et al. 1993

fl = 124; 
files(fl).tag = 'bp160609-36';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 2000;
files(fl).wt = 52000; % kg estimate, Jefferson et al. 1993

fl = 125; 
files(fl).tag = 'bp170907-41b';
files(fl).path = 'C:\tag\audit';
files(fl).sex = '';
files(fl).lnth = 2000;
files(fl).wt = 52000; % kg estimate, Jefferson et al. 1993

%% grey whales - PS, JG, JC

fl = 126; 
files(fl).tag = 'er150417-3';
files(fl).lnth = 12.76; % length estimate from Selina Agbayani 
files(fl).wt = 19924; % kg estimate from Selina Agbayani @ UBC
files(fl).age = 25; % seen since 1990

fl = 127; 
files(fl).tag = 'er160325-22';
files(fl).lnth = 12.53; % length estimate from Selina Agbayani 
files(fl).wt = 19031; % kg estimate from Selina Agbayani @ UBC
files(fl).age = 16; 

fl = 128; 
files(fl).tag = 'er160406-21';
files(fl).lnth = 12.53; % length estimate from Selina Agbayani 
files(fl).wt = 19031; % kg estimate from Selina Agbayani @ UBC
files(fl).age = 16; 

fl = 129; 
files(fl).tag = 'er160406-22';
files(fl).lnth = 12.8; % length estimate from Selina Agbayani 
files(fl).wt = 19970; % kg estimate from Selina Agbayani @ UBC
files(fl).age = 26; 

fl = 130; 
files(fl).tag = 'er160407-25';
files(fl).lnth = 12.57; % length estimate from Selina Agbayani 
files(fl).wt = 19211; % kg estimate from Selina Agbayani @ UBC
files(fl).age = 17; 

fl = 131; 
files(fl).tag = 'er160505-22';
files(fl).lnth = 12.53; % length estimate from Selina Agbayani 
files(fl).wt = 19031; % kg estimate from Selina Agbayani @ UBC
files(fl).age = 16; 

for fl = 126:131
files(fl).path = '/Users/julievanderhoop/Dropbox (Personal)/tag/audit';
files(fl).sex = '';
end

%% right whales - DPN, SP

fl = 132;
files(fl).tag = 'eg02_221d';
files(fl).sex = 'M';
files(fl).age = 11; 

fl = 133; 
files(fl).tag = 'eg05_219a';
files(fl).sex = 'F';
files(fl).age = 4; 

fl = 134; 
files(fl).tag = 'eg05_224a';
files(fl).sex = 'M';
files(fl).age = 3;

addpath '/Users/julievanderhoop/Documents/MATLAB/NARW_Foraging/'
for fl = 132:134
files(fl).path = '/Users/julievanderhoop/Dropbox (Personal)/tag/audit';
files(fl).lnth = MooreAgeLength(files(fl).age); % in cm 
files(fl).wt = 243.08*exp(0.0036*files(fl).lnth); % in kg -- weight function for NARW from Moore et al. 2005
end

%% risso's dolphins - FV

files(135).tag = 'gg15_229a';
files(136).tag = 'gg15_229c';
files(137).tag = 'gg16_171a';
files(138).tag = 'gg16_176a';
files(139).tag = 'gg17_200a';
files(140).tag = 'gg17_203a';
files(141).tag = 'gg17_213a';

for fl = 135:141
files(fl).path = '/Users/julievanderhoop/Dropbox (Personal)/tag/audit';
files(fl).sex = '';
files(fl).lnth = NaN; % looking for age-length curve
files(fl).wt = 300; % kg estimate
end

%% ziphius - FV
files(142).tag = 'zc15_212b';
files(142).sex = 'M';
files(143).tag = 'zc17_205a';
files(143).sex = 'F';
files(144).tag = 'zc17_227a';
files(144).sex = 'M';


for fl = 142:144
files(fl).path = '/Users/julievanderhoop/Dropbox (Personal)/tag/audit';
files(fl).lnth = NaN; % looking for age-length curve
files(fl).wt = 2000; % kg estimate
end

%% antarctic minkes

files(145).tag = 'bb180125-30';
files(146).tag = 'bb180227-45';
files(147).tag = 'bb180228-42';
files(148).tag = 'bb180304-40';
files(149).tag = 'bb180305-42b';
files(150).tag = 'bb190224-52';
files(151).tag = 'bb190225-55';
files(152).tag = 'bb190228-55b';
files(153).tag = 'bb190304-18';
files(154).tag = 'bb190309-52';

for fl = 145:154
files(fl).path = '/Users/julievanderhoop/Dropbox (Personal)/tag/audit';
files(fl).lnth = 800; % looking for age-length curve
files(fl).wt = 9000; % kg estimate
end


%%
% cd \\uni.au.dk\Users\au575532\Documents\MATLAB\BreathSounds\BreathCounts
save('countdata','files')