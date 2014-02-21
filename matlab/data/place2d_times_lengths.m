

nodes =     [   833            3014           6514           12637          12752           15406       ];
ltot =      [   1.671122e-01   1.413727e-01   1.316508e-01   2.101602e-02   4.297591e-02    4.195974e-02];
ltot_pe =   [   1.282769e-01   2.635241e-01   1.701552e-01   2.435289e-02   9.081522e-02    6.586360e-02];
ltot_p =    [   1.221689e+00   2.981378e-01   2.853210e-01   1.025899e+00   1.111488e-01    1.280730e-01];

t =         [   2.404631e-01   9.700774e-01   1.207850e+01    1.401057e+01  7.375932e+00    1.247893e+01];
tpe =       [   2.533275e-01   9.178206e-01   1.665013e+01    1.738294e+01  8.126292e+00    1.223591e+01];
tp =        [   1.175052e-01   4.373538e-01   4.578858e+00    5.649411e+00  2.700431e+00    3.789010e+00];

tp_savings = tpe./tp;
lp_penalty = ltot_p./ltot_pe;

%%
figure(1)
clf
loglog(nodes,ltot,'k')
hold on
loglog(nodes,ltot_pe,'b')
loglog(nodes,ltot_p,'r')
xlabel('number of nodes')
ylabel('wirelength^2')
xlim([0.5*min(nodes) 2*max(nodes)])
ylim([0.5*min(ltot) 2*max(ltot_p)])
fixfigs(1,3,14,12)

figure(2)
clf
loglog(nodes,t,'k')
hold on
loglog(nodes,tpe,'b')
loglog(nodes,tp,'r')
xlabel('number of nodes')
ylabel('solution time (s)')
xlim([0.5*min(nodes) 2*max(nodes)])
ylim([0.5*min(tp) 2*max(tpe)])
fixfigs(2,3,14,12)

figure(3)
clf
loglog(nodes,lp_penalty,'bo','markerfacecolor','b')
xlabel('number of nodes')
ylabel('wirelength^2 (approx) / wirelength^2 (exact)')
xlim([0.5*min(nodes) 2*max(nodes)])
ylim([0.5*min(lp_penalty) 2*max(lp_penalty)])
fixfigs(3,3,14,12)

figure(4)
clf
semilogx(nodes,tp_savings,'bo','markerfacecolor','b')
xlabel('number of nodes')
ylabel('speedup factor')
xlim([0.5*min(nodes) 2*max(nodes)])
ylim([1 5])
fixfigs(4,3,14,12)
