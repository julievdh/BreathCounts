% Mortola data with ours

Mortola = load('MortolaData.mat');
applyLungMass

nonRum = find(Mortola.Ruminant == 0 & Mortola.Marine == 0 & Mortola.AquaticBird == 0 & Mortola.TerrestrialBird == 0);

figure(10), clf, hold on
plot(Mortola.M(nonRum),Mortola.f(nonRum),'ko')
set(gca,'xscale','log','yscale','log'), 
xlim([10^-2 10^6]), ylim([0.1 200])
p = polyfit(log10(Mortola.M(nonRum)), log10(Mortola.f(nonRum)), 1);
x = linspace(min(Mortola.M(nonRum)),max(Mortola.M(nonRum)),100);
plot(x,10^(p(2))*x.^(p(1)),'k-')

for i = 1:length(files)
    scatter(wt(i),mnf(i),50,'ko','MarkerFaceColor',files(i).col,'MarkerFaceAlpha',0.5)
end

plot(sp_wt,sp_mnf,'ko','markerfacecolor','k','markersize',9)

[pf] = polyfit(log10(sp_wt), log10(sp_mnf), 1);
plot(sp_wt,10^(pf(2))*sp_wt.^(pf(1)),'k-')

xlabel('Body Mass (kg)'), ylabel('Mean Ventilation Frequency (/min)')
adjustfigurefont('Helvetica',14)
print -dpng -r300 BreathCounts_MortolaJVDH_2019
% add labels
% for i = unique([files(:).spp])
%     text(sp_wt(i),10^(pf(2))*sp_wt(i).^(pf(1)),files(sp_ind(i)).tag(1:2))
% end

%% IBI Mortola

figure(13), clf, hold on
plot(Mortola.M(nonRum),60./Mortola.f(nonRum),'ko')
set(gca,'xscale','log','yscale','log'), 
xlim([10^-2 10^6]), %ylim([0.1 200])
p = polyfit(log10(Mortola.M(nonRum)), log10(60./Mortola.f(nonRum)), 1);
x = linspace(min(Mortola.M(nonRum)),max(Mortola.M(nonRum)),100);
plot(x,10^(p(2))*x.^(p(1)),'k-')

for i = 1:length(files)
    scatter(wt(i),mnIBI(i),'ko','MarkerFaceColor',files(i).col,'MarkerFaceAlpha',0.5)
end

plot(sp_wt,sp_IBI,'ko','markerfacecolor','k','MarkerSize',9)

[pf] = polyfit(log10(sp_wt), log10(sp_IBI), 1);
plot(sp_wt,10^(pf(2))*sp_wt.^(pf(1)),'k-')

xlabel('Body Mass (kg)'), ylabel('Mean Inter-Breath Interval(sec)')
adjustfigurefont('Helvetica',14)
print -dpng -r300 BreathCounts_MortolaJVDH_IBI_2019

%% VE Mortola

figure(11), clf, hold on
Mortola.VT(nonRum) = 0.0681.*Mortola.M(nonRum).^1.02;

plot(Mortola.M(nonRum),Mortola.f(nonRum).*Mortola.VT(nonRum)','ko') 
set(gca,'xscale','log','yscale','log'), 
xlim([10^-2 10^5]), % ylim([0.1 200])
p = polyfit(log10(Mortola.M(nonRum)), log10(Mortola.f(nonRum).*Mortola.VT(nonRum)'), 1);
x = linspace(min(Mortola.M(nonRum)),max(Mortola.M(nonRum)),100);
plot(x,10^(p(2))*x.^(p(1)),'k-')

for i = 1:length(files)
    scatter(wt(i),mnf(i)*TLC(i),'ko','MarkerFaceColor',files(i).col,'MarkerFaceAlpha',0.5)
end

plot(sp_wt,sp_mnf.*sp_TLC,'ko','markerfacecolor','k','markersize',9)

[pVE] = polyfit(log10(sp_wt), log10(sp_mnf.*sp_TLC), 1);
plot(sp_wt,10^(pVE(2))*sp_wt.^(pVE(1)),'k-')

xlabel('Body Mass (kg)'), ylabel('Ventilation (L/min)')
adjustfigurefont('Helvetica',14), axis equal
print -dpng -r300 BreathCounts_MortolaVE_2019

%% could you plot on the same figure?
figure(12), clf, hold on

plot(sp_wt,sp_mnf.*sp_TLC,'ko','markerfacecolor','k')
plot(sp_wt,sp_IBI,'ko','markerfacecolor','k')

plot(Mortola.M(nonRum),60./Mortola.f(nonRum),'ko')
plot(Mortola.M(nonRum),Mortola.f(nonRum)'.*Mortola.VT(nonRum),'ko') % 

set(gca,'xscale','log','yscale','log')