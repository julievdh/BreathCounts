% just try this
applyLungMass

% M = 30:10:40000;
% figure(1), clf
% col = [1 0 0, 0 0 1];
% % PGLSa = 17.24; PGLSb = 0.15;
% for scaling = [0.75, 0.5]
% FMR = M.^scaling;
% FR = 10*M.^(scaling-1);
% IBI = M.^(1-scaling);
% % IBI2 = FR*60; % just check it scales the same way
% 
% 
% %loglog(M,FMR,'k'), hold on
% loglog(M,FR,'k'), hold on
% %loglog(M,IBI,'k')
% % loglog(M,IBI2,':')
% 
% end

% loglog(sp_wt,sp_IBI,'o');
figure(19), clf, hold on
for i = 1:length(files)
    scatter(wt(i),mnf(i),'ko','MarkerFaceColor',files(i).col,'markerfacealpha',0.5)
end

plot(sp_wt,sp_mnf,'ko','markerfacecolor','k','markerSize',9)
plot([10 100000],[1 1],':')
xlabel('Body Mass (kg)'), ylabel('Mean Ventilation Frequency (/min)')
set(gca,'xscale','log','yscale','log')
adjustfigurefont('Helvetica',14)
print -dpng -r300 BreathCounts_meanF

VE = mnf.*TLC; % fR * VT
figure(29),clf, hold on
for i = 1:length(files)
plot([files(i).wt],VE(i),'ko','markerfacecolor',files(i).col)
end
plot(sp_wt,sp_mnf.*sp_TLC,'ko','markersize',8,'markerfacecolor','k')
set(gca,'xscale','log','yscale','log')
xlabel('Body Mass (kg)'), ylabel('Minute Ventilation (L)')

[pVE,sVE] = polyfit(log10(sp_wt), log10(sp_mnf.*sp_TLC), 1);
plot(sp_wt,10^(pVE(2))*sp_wt.^(pVE(1)),'k-')
adjustfigurefont
print -dpng -r300 BreathCounts_VE

%% lung scaling
figure(12), clf, hold on
for i = 1:length(sp_wt)
plot(sp_wt(i),[files(sp_ind(i)).lungsc],'ko','markerfacecolor',files(sp_ind(i)).col,'markersize',8)
end
set(gca,'xscale','log')
xlabel('Body Mass (kg)'), ylabel('TLC/Body Mass')
adjustfigurefont
% print -dpng -r300 BreathCounts_LungScaling

%% plot yy with VE and IBI
figure(92), clf
loglog(sp_wt,sp_IBI,'o'), hold on
loglog(sp_wt,a3*sp_wt.^b3,'k-')
loglog(sp_wt,sp_mnf.*sp_TLC,'v')
loglog(sp_wt,10^(pVE(2))*sp_wt.^(pVE(1)),'k--')

% get fit
cf = fit(log10(sp_wt)',log10(sp_mnf.*sp_TLC)','poly1');
cf_coeff = coeffvalues(cf);
cf_confint = confint(cf);
a = 10.^(cf_coeff(2)); b = cf_coeff(1);
b_uncert = (cf_confint(2,1) - cf_confint(1,1))/2; 
a_uncert = (cf_confint(2,2) - cf_confint(1,2))/2;

return
%% resps per day
% resp frequency is per minute, 
figure, loglog(wt,mnf*60*24,'.')
hold on
loglog(sp_wt,sp_mnf*60*24,'ko','markerfacecolor','k')
loglog([10 100000],[1*60*24 1*60*24],':')
xlabel('Body Mass (kg)'), ylabel('Ventilations per day')
