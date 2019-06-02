clear, close all

% load count data from tags
load('countdata')
assignSpp % assign species codes 

extractIBI % extract values from files struct

wt = [files(:).wt];
dur = [files(:).dur];
spp = [files(:).spp]; 

%% calculate species means
% tabulate -- but only where duration exceeds duration threshold 

for i = unique(spp(~isnan(spp)))
    sp_wt(i) = nanmean(wt(find(spp == i)));
    sp_sdwt(i) = nanstd(wt(find(spp == i)));
    sp_IBI(i) = nanmean(mnIBI(find(spp == i)));
    sp_mdIBI(i) = nanmean(mdIBI(find(spp == i)));
    sp_mnf(i) = nanmean(mnf(find(spp == i)));
    sp_mdf(i) = nanmean(mdf(find(spp == i)));
    sp_dur(i) = sum(dur(find(spp == i))); % total hours per species
    sp_ct(i) = length(find(spp == i)); 
    sp_ind(i) = find(spp == i,1); % find first appearance of species index
end

tab = table(sp_ct',sp_wt',sp_sdwt',sp_dur')

% order for plotting for phylogeny
phyorder = [16,15,12,11,14,9,13,8,10,7,5,6,1,2,3,4];
assignCol % assign colours to species

%% fit polynomial
[p3,s] = polyfit(log10(sp_wt), log10(sp_IBI), 1);
f3 = polyval(p3,log10(sp_wt)); 
resid = log10(sp_IBI)-f3; 
% retrieve parameters
b3 = p3(1); a3 = 10^(p3(2));

% plot residuals
figure(2), clf, hold on
plot(resid,phyorder,'o','LineWidth',1)
plot([0 0],[1 16],'k:')
adjustfigurefont('Helvetica',14)
grid on
set(gcf,'position',[730   119   293   670],'paperpositionmode','auto')
print -dpng -r300 BreathCounts_Resid_phylo

% plot
figure(12), clf
loglog(sp_wt, sp_IBI, '.', sp_wt, a3*sp_wt.^b3, '-'), hold on
loglog(sp_wt,sp_mdIBI,'x')

% evaluate fit
cf = fit(log10(sp_wt)',log10(sp_IBI)','poly1');
cf_coeff = coeffvalues(cf);
cf_confint = confint(cf);
a = 10.^(cf_coeff(2)); b = cf_coeff(1);
b_uncert = (cf_confint(2,1) - cf_confint(1,1))/2;
a_uncert = (cf_confint(2,2) - cf_confint(1,2))/2;
%%
figure(3), clf, hold on % have to fix errorbar on y axis
set(gca, 'XScale', 'log', 'yscale','log')
xlabel('Body Mass (kg)'), ylabel('IBI'), ylim([10 10^4])
set(gcf,'position',[360.3333  197.6667  913.3333  420.0000], 'paperpositionmode','auto')


% for i = 1:length(files)
%     h = scatter(repmat(files(i).wt,length(files(i).tdiff),1)+(randn(1,length(files(i).tdiff))*files(i).wt/50)',files(i).tdiff,'markeredgecolor',files(i).col); drawnow;
%     h.MarkerEdgeAlpha = 0.2;
%
%     % mindiff(i) = min(files(i).tdiff);
% end
odont = 1:10; 
orig = [1 2 5 6 7 8 10];
myst = 11:16; 

 % plot original 7 
for i = 1:length(files)
    if ismember(files(i).spp,orig) == 1 
    h = scatter(files(i).wt,mnIBI(i),'o','markerfacecolor',files(i).col,'markeredgecolor','k');
    h.MarkerEdgeAlpha = 0.5; h.MarkerFaceAlpha = 0.5;
    end
    % plot(files(i).wt,mdIBI(i),'kx')
end

% add species means
plot(sp_wt(orig),sp_IBI(orig),'ko','markersize',10,'markerfacecolor','k')
% add regression line from PGLS
PGLSa = 10^(0.836); PGLSb = 0.267;
plot(sp_wt, PGLSa*sp_wt.^PGLSb, 'k-') % USE PGLS VALUES
plot(sort(sp_wt), PGLSa*sort(sp_wt).^0.25,'k:')
plot(sort(sp_wt), PGLSa*sort(sp_wt).^0.5,'k:')
%plot(sp_wt([1 2 6 4]), a3*sp_wt([1 2 6 4]).^(0.25*1.03),'k:')
%plot(sp_wt([7 5 3]), a3*sp_wt([7 5 3]).^(0.25*0.91),'k:')
%plot(sp_wt([1 2 6 4]), a3*sp_wt([1 2 6 4]).^(0.5*1.03),'k:')
%plot(sp_wt([7 5 3]), a3*sp_wt([7 5 3]).^(0.5*0.91),'k:')
adjustfigurefont('Helvetica',14)
print -dpng -r300 BreathCounts_MeanMedian_Lines_n7
%%
figure(3), clf, hold on % have to fix errorbar on y axis
set(gca, 'XScale', 'log', 'yscale','log')
xlabel('Body Mass (kg)'), ylabel('IBI'), ylim([10 10^4])
set(gcf,'position',[360.3333  197.6667  913.3333  420.0000], 'paperpositionmode','auto')

% plot all odontocetes
for i = 1:length(files)
    if ismember(files(i).spp,odont) == 1 
    h = scatter(files(i).wt,mnIBI(i),'o','markerfacecolor',files(i).col,'markeredgecolor','k');
    h.MarkerEdgeAlpha = 0.5; h.MarkerFaceAlpha = 0.5;
    end
    % plot(files(i).wt,mdIBI(i),'kx')
end


plot(sp_wt(odont),sp_IBI(odont),'ko','markersize',10,'markerfacecolor','k')
cf = fit(log10(sp_wt(odont))',log10(sp_IBI(odont))','poly1');
cf_coeff = coeffvalues(cf);
cf_confint = confint(cf);
a = 10.^(cf_coeff(2)); b = cf_coeff(1);
plot(sort(sp_wt(odont)), a*sort(sp_wt(odont)).^b, 'k-') % USE PGLS VALUES
adjustfigurefont('Helvetica',14)
plot(sort(sp_wt), PGLSa*sort(sp_wt).^0.25,'k:')
plot(sort(sp_wt), PGLSa*sort(sp_wt).^0.5,'k:')

print -dpng -r300 BreathCounts_MeanMedian_Lines_odont
%% 
% plot all mysticetes
for i = 1:length(files)
    if ismember(files(i).spp,myst) == 1 
    h = scatter(files(i).wt,mnIBI(i),'o','markerfacecolor',files(i).col,'markeredgecolor','k');
    h.MarkerEdgeAlpha = 0.5; h.MarkerFaceAlpha = 0.5;
    end
    % plot(files(i).wt,mdIBI(i),'kx')
end

% plot(sp_wt, a3*sp_wt.^b3,'k--')

plot(sp_wt(myst),sp_IBI(myst),'ko','markersize',10)
cf = fit(log10(sp_wt)',log10(sp_IBI)','poly1');
cf_coeff = coeffvalues(cf);
cf_confint = confint(cf);
a = 10.^(cf_coeff(2)); b = cf_coeff(1);
plot(sort(sp_wt), a*sort(sp_wt).^b, 'k--') 

print -dpng -r300 BreathCounts_MeanMedian_Lines_all


%% CDF plot -- work on this to normalize duration
figure(8), clf, hold on
for i = 1:length(files)
    h = plot(sort(files(i).tdiff),(1:length(files(i).tdiff))/length(files(i).tdiff));
    set(h,'color',files(i).col)
end
set(gca,'xscale','log'),
set(gcf,'position',[360.3333  197.6667  913.3333  420.0000], 'paperpositionmode','auto')
ylabel('Proportion')
xlabel('IBI (sec)'), grid on
patch([1 10^4 10^4 1 1], [0 0 1 1 0],[1 1 1],'facealpha',0.6)

% Add species averages on top of semitransparent individual samples
% initialize new structure
for i = unique(spp(~isnan(spp))), spps(i).tdiff = NaN; end
% vertcat with all tdiff from each file based on spp
for fl = 1:length(files)
    if ~isnan(spp(fl)); 
    spps(spp(fl)).tdiff = vertcat(spps(spp(fl)).tdiff, files(fl).tdiff);
    spps(spp(fl)).col = files(fl).col;
    end
end

for i = unique(spp(~isnan(spp))),
    h = plot(sort(spps(i).tdiff),(1:length(spps(i).tdiff))/length(spps(i).tdiff));
    set(h,'color',spps(i).col,'linewidth',2)
end

print -dpng -r300 BreathCounts_CDF_all

figure(81), clf, hold on
for i = unique(spp(~isnan(spp))),
    h = plot(sort(spps(i).tdiff),(1:length(spps(i).tdiff))/length(spps(i).tdiff));
    set(h,'color',spps(i).col,'linewidth',2)
end
set(gca,'xscale','log'), grid on,
set(gcf,'position',[360.3333  197.6667  913.3333  420.0000], 'paperpositionmode','auto')
xlabel('Inter-Breath Interval (sec)'), ylabel('Proportion')

print -dpng -r300 BreathCounts_CDF_spp

% do for just gm11_148c for talk
% figure(82), clf
% i = 59;
% h = plot(sort(files(i).tdiff),(1:length(files(i).tdiff))/length(files(i).tdiff));
% set(h,'color',files(i).col,'LineWidth',2)
% set(gca,'xscale','log')
% ylabel('Proportion')
% xlabel('IBI (sec)'), grid on
% 

%% massage for output for R: weight, spp weight, IBI, species code, file index, weight+file code
% all = horzcat(repmat(files(1).wt,length(files(1).tdiff),1),... % individual weight
%     repmat(sp_wt(files(1).spp),length(files(1).tdiff),1),... % species average weight
%     files(1).tdiff,... % inter-breath intervals
%     repmat(files(1).spp,length(files(1).tdiff),1),... % species code
%     repmat(1,length(files(1).tdiff),1),... % file number
%     repmat(round(files(1).wt)+(1/10),length(files(1).tdiff),1)); % unique identifier
% for i = 2:length(files)
%     test = horzcat(repmat(files(i).wt,length(files(i).tdiff),1),... % individual weight
%         repmat(sp_wt(files(i).spp),length(files(i).tdiff),1),... % species weight
%         files(i).tdiff,... % inter-breath intervals
%         repmat(files(i).spp,length(files(i).tdiff),1),... % species code
%         repmat(i,length(files(i).tdiff),1),... % file number
%         repmat(round(files(i).wt)+(i/10),length(files(i).tdiff),1)); % unique identifier
%     all = vertcat(all,test);
% end
% 
% csvwrite('JoyPlotBreathDatTest_export.csv',all,1,0)

%% export for PGLS

csvwrite('BreathCounts_PGLSdata',[sp_wt; sp_IBI; sp_mdIBI; unique(spp(~isnan(spp))); sp_mnf; sp_ct]',1,0)

% csvwrite('BreathCounts_PGLSdata',[sp_wt; sp_IBI; sp_mdIBI; unique([files(:).spp]); sp_mnf; sp_TLC]',1,0)
