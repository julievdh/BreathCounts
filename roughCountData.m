clear, close all

% load count data from tags
load('countdata')

for i = 1:length(files)
    if isempty(files(i).tdiff) == 0
        % calculate duration of resp audit in hours
        files(i).dur = (files(i).resp(end) - files(i).resp(1))/3600;
        if files(i).dur >= 6
            mnIBI(i) = mean(files(i).tdiff); % mean IBI
            stdIBI(i) = std(files(i).tdiff); % STD
            mdIBI(i) = median(files(i).tdiff); % median IBI
            modeIBI(i) = mode(files(i).tdiff); % mode
            
            % if isfield(files(i),'swim_ct') == 1
            mnf(i) = 60/mnIBI(i); % mean frequency
            mdf(i) = 60/mdIBI(i); % median frequency
            % end
        else
            mnIBI(i) = NaN; % mean IBI
            stdIBI(i) = NaN; % STD
            mdIBI(i) = NaN; % median IBI
            modeIBI(i) = NaN; % mode
        end
        else files(i).dur = NaN;
    end
end

wt = [files(:).wt];
dur = [files(:).dur];

%% calculate species means
assignSpp % assign species codes
for i = unique([files(:).spp])
    sp_wt(i) = nanmean(wt(find([files(:).spp] == i)));
    sp_IBI(i) = nanmean(mnIBI(find([files(:).spp] == i)));
    sp_mdIBI(i) = nanmean(mdIBI(find([files(:).spp] == i)));
    sp_mnf(i) = nanmean(mnf(find([files(:).spp] == i)));
    sp_mdf(i) = nanmean(mdf(find([files(:).spp] == i)));
    sp_dur(i) = sum(dur(find([files(:).spp] == i))); % total hours per species
    sp_ct(i) = length(find([files(:).spp] == i));
    sp_ind(i) = find([files(:).spp] == i,1); % find first appearance of species index
end
assignCol % assign colours to species

%% fit polynomial
[p3,s] = polyfit(log10(sp_wt), log10(sp_IBI), 1);

% retrieve parameters
b3 = p3(1); a3 = 10^(p3(2));

% plot
figure(2), clf
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
for i = 1:length(files)
    h = scatter(files(i).wt,mnIBI(i),'o','markerfacecolor',files(i).col,'markeredgecolor','k');
    h.MarkerEdgeAlpha = 0.5; h.MarkerFaceAlpha = 0.5;
    
    % plot(files(i).wt,mdIBI(i),'kx')
end

% add species means
plot(sp_wt,sp_IBI,'ko','markersize',10,'markerfacecolor','k')
% add regression line from PGLS
% PGLSa = 10^(0.836); PGLSb = 0.267;
% plot(sp_wt, PGLSa*sp_wt.^PGLSb, 'k-') % USE PGLS VALUES
plot(sp_wt, a3*sp_wt.^0.25,'k:')
plot(sp_wt, a3*sp_wt.^0.5,'k:')
%plot(sp_wt([1 2 6 4]), a3*sp_wt([1 2 6 4]).^(0.25*1.03),'k:')
%plot(sp_wt([7 5 3]), a3*sp_wt([7 5 3]).^(0.25*0.91),'k:')
%plot(sp_wt([1 2 6 4]), a3*sp_wt([1 2 6 4]).^(0.5*1.03),'k:')
%plot(sp_wt([7 5 3]), a3*sp_wt([7 5 3]).^(0.5*0.91),'k:')
plot(sp_wt, a3*sp_wt.^b3,'k--')


print -dpng -r300 BreathCounts_MeanMedian_Lines

p3%% CDF plot -- work on this to normalize duration
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
for i = unique([files(:).spp]), spp(i).tdiff = NaN; end
% vertcat with all tdiff from each file based on spp
for fl = 1:length(files)
    species = files(fl).spp;
    spp(species).tdiff = vertcat(spp(species).tdiff, files(fl).tdiff);
    spp(species).col = files(fl).col;
end

for i = unique([files(:).spp]),
    h = plot(sort(spp(i).tdiff),(1:length(spp(i).tdiff))/length(spp(i).tdiff));
    set(h,'color',spp(i).col,'linewidth',2)
end

print -dpng -r300 BreathCounts_CDF_all

figure(81), clf, hold on
for i = unique([files(:).spp]),
    h = plot(sort(spp(i).tdiff),(1:length(spp(i).tdiff))/length(spp(i).tdiff));
    set(h,'color',spp(i).col,'linewidth',2)
end
set(gca,'xscale','log'), grid on,
set(gcf,'position',[360.3333  197.6667  913.3333  420.0000], 'paperpositionmode','auto')
xlabel('Inter-Breath Interval (sec)'), ylabel('Proportion')

print -dpng -r300 BreathCounts_CDF_spp

% do for just gm11_148c for talk
figure(82), clf
i = 59;
h = plot(sort(files(i).tdiff),(1:length(files(i).tdiff))/length(files(i).tdiff));
set(h,'color',files(i).col,'LineWidth',2)
set(gca,'xscale','log')
ylabel('Proportion')
xlabel('IBI (sec)'), grid on
%% massage for output for R: weight, spp weight, IBI, species code, file index, weight+file code
all = horzcat(repmat(files(1).wt,length(files(1).tdiff),1),... % individual weight
    repmat(sp_wt(files(1).spp),length(files(1).tdiff),1),... % species average weight
    files(1).tdiff,... % inter-breath intervals
    repmat(files(1).spp,length(files(1).tdiff),1),... % species code
    repmat(1,length(files(1).tdiff),1),... % file number
    repmat(round(files(1).wt)+(1/10),length(files(1).tdiff),1)); % unique identifier
for i = 2:length(files)
    test = horzcat(repmat(files(i).wt,length(files(i).tdiff),1),... % individual weight
        repmat(sp_wt(files(i).spp),length(files(i).tdiff),1),... % species weight
        files(i).tdiff,... % inter-breath intervals
        repmat(files(i).spp,length(files(i).tdiff),1),... % species code
        repmat(i,length(files(i).tdiff),1),... % file number
        repmat(round(files(i).wt)+(i/10),length(files(i).tdiff),1)); % unique identifier
    all = vertcat(all,test);
end

csvwrite('JoyPlotBreathDatTest_export.csv',all,1,0)

%% export for PGLS

% csvwrite('BreathCounts_PGLSdata',[sp_wt; sp_IBI; sp_mdIBI; unique([files(:).spp]); sp_mnf; sp_TLC]',1,0)
