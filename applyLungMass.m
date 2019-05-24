% applyLungMass

for fl = 1:length(files)
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
end
TLC = [files(:).TLC]; lungmass = [files(:).lungmass]; 
lungsc = [files(:).lungsc];

% MAKE A SCALING VECTOR FOR LUNG MASS

%% get relationship for TLC to BM
for i = unique([files(:).spp])
    sp_wt(i) = mean(wt(find([files(:).spp] == i)));
    sp_TLC(i) = mean(TLC(find([files(:).spp] == i)));
    sp_lungmass(i) = mean(lungmass(find([files(:).spp] == i)));
end

[p_dp,s_dp] = polyfit(log10(sp_wt([1 2 6 4])),log10(sp_TLC([1 2 6 4])),1);
[p_kzp,s_kpz] = polyfit(log10(sp_wt([7 5 3])),log10(sp_TLC([7 5 3])),1);
figure(8), clf, hold on
for i = 1:length(files)
    plot(files(i).wt, files(i).TLC,'ko','markerfacecolor',files(i).col)
end
xlabel('Body Mass (kg)'), ylabel('Total Lung Capacity (L)')
set(gca,'xscale','log','yscale','log')
plot(sp_wt,sp_TLC,'ko','markersize',8,'markerfacecolor','k')
% add labels
for i = unique([files(:).spp])
    text(sp_wt(i),sp_TLC(i)*(i/3),files(sp_ind(i)).tag(1:2))
end


print -dpng -r300 BreathCounts_TLCscaling
