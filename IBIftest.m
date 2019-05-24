ag
% test how IBI should scale with W 

% using data from Mortola 2006
% beluga, tursiops, fin whale, Hp, bottlenose whale, orca, pilot whale
% W = [1135 159 51700 52 6650 2883 1034];
% f = [1.3 2.3 0.8 1 1.1 0.9 2.5];

% use non-ruminant mammals so actually have proper relationship to be sure
% are testing the right thing
W = [11 12 500 63 125 145 5.2 37 85 39 77 274 2.3 4.9 2.0 1.2 3.1 3.8 5.0 0.4 620 278 205 235 3091 26 217 12 2611 4550];
f = [17 22 13 9 15 17 32 13 36 21 38 15 28 28 36 37 21 63 18 56 12 15 23 23 11 17 17 15 5 7]; 

% make hypothetical
Wrange = 10:10:6000;
fest = Wrange.^(-0.25);

figure(1), clf
loglog(W,f,'o'), hold on
loglog(Wrange,fest)

%% regress
p = polyfit(log(W), log(f), 1);

% retrieve original parameters
b = p(1);
a = exp(p(2));

% plot
loglog(W, f, '.', W, a*W.^b, 'r')
%%
p2 = polyfit(log(W), log(1./f), 1);

% retrieve original parameters
b2 = p2(1);
a2 = exp(p2(2));

% plot
figure(2)
loglog(W, 1./f, '.', W, a2*W.^b2, 'r'), hold on

loglog(W,1./f,'o')
loglog(Wrange,1./fest)


