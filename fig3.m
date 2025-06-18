figure(31)
e=linspace(-1.5,3,50);
ax(1)=subplot(2,2,1)
histogram(log10(fluc.cv(:,1)),e,'EdgeColor','none')
box off; set(gca,'TickDir','out')
ax(2)=subplot(2,2,3)
histogram(log10(fluc.cv(:,2)),e,'EdgeColor','none')
box off; set(gca,'TickDir','out')
linkaxes(ax)

e=linspace(0,600,25);
axc(1)=subplot(2,2,2)
histogram(fluc.optimal_window_size(:,2),e,'EdgeColor','none')
box off; set(gca,'TickDir','out')
axc(2)=subplot(2,2,4)
histogram(fluc.optimal_window_size(:,1),e,'EdgeColor','none')
box off; set(gca,'TickDir','out')
linkaxes(axc)

prctile(fluc.optimal_window_size,[25 50 75])
prctile(fluc.cv,[25 50 75])