% applyLungMass

for fl = 1:length(files)
    files(fl).TLC = 0.135*(files(fl).wt).^0.92; % estimate from Kooyman 
    % but if other data exist, use spp-specific:
    if     strfind(files(fl).tag,'tt') == 1 % tursiops
        files(fl).lungmass = 0.037*files(fl).wt.^(0.948);
        files(fl).TLC = 0.081*files(fl).wt; % piscitelli et al. 2013 supplemental 1
        files(fl).lungsc = 0.081;
    end
    if     strfind(files(fl).tag,'hp') == 1 % phocoena
        files(fl).lungmass = 0.037*files(fl).wt.^(0.948);
        files(fl).TLC = 0.091*files(fl).wt; % piscitelli et al. 2013 supplemental 1
        files(fl).lungsc = 0.091;
    end
    if     strfind(files(fl).tag,'pw') == 1 % short finned pilot whale
        files(fl).lungmass = 0.037*files(fl).wt.^(0.948);
        files(fl).TLC = 0.100*files(fl).wt; % piscitelli et al. 2013 supplemental 1
        files(fl).lungsc = 0.100;
    end
    if     strfind(files(fl).tag,'gm') == 1 % long finned pilot whale
        files(fl).lungmass = 0.037*files(fl).wt.^(0.948);
        files(fl).TLC = 0.100*files(fl).wt; % piscitelli et al. 2013 supplemental 1
        files(fl).lungsc = 0.100;
    end
    if     strfind(files(fl).tag,'sw') == 1 % sperm whale
        files(fl).lungmass = 0.019*files(fl).wt.^(0.900);
        files(fl).TLC = 0.019*files(fl).wt; % piscitelli et al. 2013 supplemental 1
        files(fl).lungsc = 0.019;
    end
    if     strfind(files(fl).tag,'zc') == 1 % cuvier's beaked whale
        files(fl).lungmass = 0.019*files(fl).wt.^(0.900);
        files(fl).TLC = 0.026*files(fl).wt; % piscitelli et al. 2013 supplemental 1
        files(fl).lungsc = 0.026;
    end
    if     strfind(files(fl).tag,'md') == 1 % baird's beaked whale
        files(fl).lungmass = 0.019*files(fl).wt.^(0.900);
        files(fl).TLC = 0.026*files(fl).wt; % piscitelli et al. 2013 supplemental 1
        files(fl).lungsc = 0.026;
    end
    if strfind(files(fl).tag,'dl') == 1 & strfind(files(fl).sex,'M') % beluga
        files(fl).TLC = 0.0837*files(fl).wt;
    end
    if strfind(files(fl).tag,'dl') == 1 & strfind(files(fl).sex,'F') % beluga
        files(fl).TLC = 0.117*files(fl).wt;
    end
    if strfind(files(fl).tag,'ha') == 1 % bottlenosed whale
        files(fl).TLC = 0.025*files(fl).wt;
    end
    if strfind(files(fl).tag,'bp') == 1 % fin whale
        files(fl).TLC = 0.048*files(fl).wt;
    end
end
TLC = [files(:).TLC]; lungmass = [files(:).lungmass];
lungsc = [files(:).lungsc];

% MAKE A SCALING VECTOR FOR LUNG MASS

%% get relationship for TLC to BM
for i = unique(spp(~isnan(spp)))
    sp_wt(i) = mean(wt(find([files(:).spp] == i)));
    sp_TLC(i) = mean(TLC(find([files(:).spp] == i)));
end

[p_dp,s_dp] = polyfit(log10(sp_wt([1 2 6 4])),log10(sp_TLC([1 2 6 4])),1); % fit porpoise and dolphin
[p_kzp,s_kpz] = polyfit(log10(sp_wt([7 5 3])),log10(sp_TLC([7 5 3])),1); % fit deep-divers
figure(8), clf, hold on
for i = 1:length(files)
    scatter(files(i).wt, files(i).TLC,'ko','markerfacecolor',files(i).col,'markerfacealpha',0.5)
end
xlabel('Body Mass (kg)'), ylabel('Total Lung Capacity (L)')
set(gca,'xscale','log','yscale','log')
plot(sp_wt,sp_TLC,'ko','markersize',9,'markerfacecolor','k')
% add labels
% for i = unique([files(:).spp])
%    text(sp_wt(i),sp_TLC(i)*(i/3),files(sp_ind(i)).tag(1:2))
% end


%% data from Tenney and Remmers 1963:
% body mass and lung volume (L) 
tenneyremmers = [0.011185056	4.4005254E-4
0.016005902	6.274557E-4
0.041005343	9.135348E-4
0.16814382	0.011177338
0.35212156	0.014970165
1.3806317	0.16501471
2.2098255	0.16501471
3.7828665	0.13675766
2.8912756	0.23042934
2.0662131	0.31512976
4.3269987	0.6144968
12.967029	1.6732659
23.74031	3.9370067
36.33384	5.164056
68.02755	6.3624644
69.56847	12.150444
390.33698	12.935478
139.3057	10.720412
199.34755	22.72456
266.72888	47.175446]; 

plot(tenneyremmers(:,1),tenneyremmers(:,2),'ko'), set(gca,'xscale','log','yscale','log')

cf = fit(log10(tenneyremmers(:,1)),log10(tenneyremmers(:,2)),'poly1');
cf_coeff = coeffvalues(cf);
a = 10.^(cf_coeff(2)); b = cf_coeff(1);
plot([min(tenneyremmers(:,1)) max(wt)], a*([min(tenneyremmers(:,1)) max(wt)]).^b, 'k--')

adjustfigurefont('Helvetica',14)

print -dpng -r300 BreathCounts_TLCscaling 
