% Mortola data with ours

Mortola = load('MortolaData');

nonRum = find(Mortola.Ruminant == 0 & Mortola.Marine == 0 & Mortola.AquaticBird == 0 & Mortola.TerrestrialBird == 0);

figure(10), hold on
plot(Mortola.M(nonRum),Mortola.f(nonRum),'ko')
set(gca,'xscale','log','yscale','log'), 
xlim([10^-2 10^6]), ylim([0.1 200])
p = polyfit(log10(Mortola.M(nonRum)), log10(Mortola.f(nonRum)), 1);
x = linspace(min(Mortola.M(nonRum)),max(Mortola.M(nonRum)),100);
plot(x,10^(p(2))*x.^(p(1)),'k-')

for i = 1:length(files)
    plot(wt(i),mnf(i),'ko','MarkerFaceColor',files(i).col)
end

plot(sp_wt,sp_mnf,'ko','markerfacecolor','k')

[pf] = polyfit(log10(sp_wt), log10(sp_mnf), 1);
plot(sp_wt,10^(pf(2))*sp_wt.^(pf(1)),'k-')

xlabel('Body Mass (kg)'), ylabel('Mean Ventilation Frequency (/min)')
adjustfigurefont
print -dpng -r300 BreathCounts_MortolaJVDH

%% IBI Mortola

figure(13), hold on
plot(Mortola.M(nonRum),60./Mortola.f(nonRum),'ko')
set(gca,'xscale','log','yscale','log'), 
xlim([10^-2 10^6]), %ylim([0.1 200])
p = polyfit(log10(Mortola.M(nonRum)), log10(60./Mortola.f(nonRum)), 1);
x = linspace(min(Mortola.M(nonRum)),max(Mortola.M(nonRum)),100);
plot(x,10^(p(2))*x.^(p(1)),'k-')

for i = 1:length(files)
    plot(wt(i),mnIBI(i),'ko','MarkerFaceColor',files(i).col)
end

plot(sp_wt,sp_IBI,'ko','markerfacecolor','k')

[pf] = polyfit(log10(sp_wt), log10(sp_IBI), 1);
plot(sp_wt,10^(pf(2))*sp_wt.^(pf(1)),'k-')

xlabel('Body Mass (kg)'), ylabel('Mean Inter-Breath Interval(sec)')
adjustfigurefont
print -dpng -r300 BreathCounts_MortolaJVDH_IBI

%% VE Mortola

figure(11), clf, hold on
plot(Mortola.M(nonRum),Mortola.f(nonRum).*Mortola.M(nonRum)/10,'ko') % ASSUMPTION THAT VT = 1/10 mass in nonRuminants -- CHECK (we know it at least scales as ^1)
set(gca,'xscale','log','yscale','log'), 
xlim([10^-2 10^5]), % ylim([0.1 200])
p = polyfit(log10(Mortola.M(nonRum)), log10(Mortola.f(nonRum).*Mortola.M(nonRum)/10), 1);
x = linspace(min(Mortola.M(nonRum)),max(Mortola.M(nonRum)),100);
plot(x,10^(p(2))*x.^(p(1)),'k-')

for i = 1:length(files)
    plot(wt(i),mnf(i)*TLC(i),'ko','MarkerFaceColor',files(i).col)
end

plot(sp_wt,sp_mnf.*sp_TLC,'ko','markerfacecolor','k')

[pVE] = polyfit(log10(sp_wt), log10(sp_mnf.*sp_TLC), 1);
plot(sp_wt,10^(pVE(2))*sp_wt.^(pVE(1)),'k-')

xlabel('Body Mass (kg)'), ylabel('Minute Volume (L)')
adjustfigurefont, axis equal
print -dpng -r300 BreathCounts_MortolaVE 

%% could you plot on the same figure?
figure(12), clf, hold on

plot(sp_wt,sp_mnf.*sp_TLC,'ko','markerfacecolor','k')
plot(sp_wt,sp_IBI,'ko','markerfacecolor','k')

plot(Mortola.M(nonRum),60./Mortola.f(nonRum),'ko')
Mortola.VT(nonRum) = 7.96.*Mortola.M(nonRum).^1.04;
plot(Mortola.M(nonRum),Mortola.f(nonRum)'.*Mortola.VT(nonRum),'ko') % ASSUMPTION THAT VT = 1/10 mass in nonRuminants -- CHECK (we know it at least scales as ^1)

set(gca,'xscale','log','yscale','log')