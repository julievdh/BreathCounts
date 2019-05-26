% Different spp breath spectrograms

settagpath('audio','D:/')
tag = 'hp13_170a'; 

R = loadaudit(tag); 
R = d3tagaudit(tag,2145,R,20000);

%% beluga

tag = 'dl16_135a'; % 'tt17_135b'; %  
R = loadaudit(tag); 
R = d3tagaudit(tag,90026,R,20000);

%% 
settagpath('Audio','O:\ST_Bio-Kaskelot\beaked whales\')
tag = 'zc10_272a'; 

R = loadaudit(tag); 
R = d3tagaudit(tag,2145,R,20000);



