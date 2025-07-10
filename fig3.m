load('results/optimal_fluctuation_windows.mat')

figure(31)
e=linspace(-1.5,3,50);
ax(1)=subplot(2,2,1)
histogram(log10(fluc.cv(:,1)),e,'EdgeColor','none')
box off; set(gca,'TickDir','out')
ax(2)=subplot(2,2,3)
histogram(log10(fluc.cv(:,2)),e,'EdgeColor','none')
box off; set(gca,'TickDir','out')
linkaxes(ax)

e=linspace(log10(60),log10(600),12);
axc(1)=subplot(2,2,2)
histogram(log10(fluc.optimal_window_size(:,1)),e,'EdgeColor','none')
box off; set(gca,'TickDir','out')
axc(2)=subplot(2,2,4)
histogram(log10(fluc.optimal_window_size(:,2)),e,'EdgeColor','none')
box off; set(gca,'TickDir','out')
linkaxes(axc)

prctile(fluc.optimal_window_size,[25 50 75])
prctile(fluc.cv,[25 50 75])

%%

figure(32)
subplot(2,3,1)
loglog(syn_table.syn_efficacy,fluc.cv,'o')
subplot(2,3,2)
loglog(syn_table.pre_firing_rate,fluc.cv,'o')
subplot(2,3,3)
loglog(syn_table.post_firing_rate,fluc.cv,'o')
subplot(2,3,4)
loglog(syn_table.syn_efficacy,fluc.optimal_window_size,'o')
subplot(2,3,5)
loglog(syn_table.pre_firing_rate,fluc.optimal_window_size,'o')
subplot(2,3,6)
loglog(syn_table.post_firing_rate,fluc.optimal_window_size,'o')