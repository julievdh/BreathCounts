tag = 'sw05_199c';
prefix = strcat(tag(1:2),tag(6:9));
recdir = strcat(gettagpath('AUDIO'),'/',tag(1:4),'/',tag);
R = loadaudit(tag);

[cues,breaths] = findbreathcues(R); 