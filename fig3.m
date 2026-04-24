load("results\optimal_fluctuation_windows.mat")
load('results\syn_index.mat')

%% CV and optimal window size results for data and controls

prctile(fluc.optimal_window_size,[25 50 75])
prctile(fluc.cv,[25 50 75])

figure(31)
subplot(2,1,1)
e=linspace(-1.5,3,1024);
[f,x] = ksdensity(log10(fluc.cv(:,1)),e);
shadedErrorBar(x,f*0,[f; f*0],'b');
hold on
[fcontrol,x] = ksdensity(log10(fluc.cv(:,2)),e);
shadedErrorBar(x,fcontrol*0,[fcontrol; fcontrol*0],'r');
box off; set(gca,'TickDir','out')
hold off
xlabel('CV')

subplot(2,1,2)
xtick = [log10([60 120 240 480])];
xticklabels = arrayfun(@num2str,(10.^xtick)/60,'UniformOutput',false);
e=linspace(log10(60)-.1,log10(600)+.1,1024);
[f,x] = ksdensity(log10(fluc.optimal_window_size(:,1)),e,'Support',[min(e) max(e)],'BoundaryCorrection','reflection');
shadedErrorBar(x,f*0,[f; f*0],'b');
box off; set(gca,'TickDir','out')
hold on
[fcontrol,x] = ksdensity(log10(fluc.optimal_window_size(:,2)),e,'Support',[min(e) max(e)],'BoundaryCorrection','reflection');
shadedErrorBar(x,fcontrol*0,[fcontrol; fcontrol*0],'r');
box off; set(gca,'TickDir','out')
hold off
set(gca,'XTick',xtick)
set(gca,'XTickLabel',xticklabels)
xlabel(['Optimal Window Size' ...
    ''])
axis tight

%%

uws1 = min(diff(log10(unique(fluc.optimal_window_size(:,1)))));
uws2 = min(diff(log10(unique(fluc.optimal_window_size(:,2)))));
jit = (randn(size(fluc.optimal_window_size))).*[uws1 uws2];


figure(32)
axc(2)=subplot(3,1,1)
plot(log10(fluc.optimal_window_size)+jit,log10(syn_table.syn_efficacy),'.')
set(gca,'XTick',tick)
set(gca,'XTickLabel','')
ylabel('Efficacy')
box off; set(gca,'TickDir','out')
addLinTrend(log10(fluc.optimal_window_size(:,1)),log10(syn_table.syn_efficacy),'b')
addLinTrend(log10(fluc.optimal_window_size(:,2)),log10(syn_table.syn_efficacy),'r')

axc(3)=subplot(3,1,2)
plot(log10(fluc.optimal_window_size)+jit,log10(syn_table.pre_firing_rate),'.')
ylabel('Pre Firing Rate')
set(gca,'XTick',xtick)
set(gca,'XTickLabel','')
box off; set(gca,'TickDir','out')
addLinTrend(log10(fluc.optimal_window_size(:,1)),log10(syn_table.pre_firing_rate),'b')
addLinTrend(log10(fluc.optimal_window_size(:,2)),log10(syn_table.pre_firing_rate),'r')
ylim([-1 2])

axc(4)=subplot(3,1,3)
plot(log10(fluc.optimal_window_size)+jit,log10(syn_table.post_firing_rate),'.')
ylabel('Post Firing Rate')
xlabel('Optimal Window Size [min]')
set(gca,'XTick',xtick)
% set(gca,'XTickLabel',xticklabels)
box off; set(gca,'TickDir','out')
linkaxes(axc,'x')
xlim([log10(60) log10(600)])
addLinTrend(log10(fluc.optimal_window_size(:,1)),log10(syn_table.post_firing_rate),'b')
addLinTrend(log10(fluc.optimal_window_size(:,2)),log10(syn_table.post_firing_rate),'r')
ylim([-1 2])

%% flip axes


uws1 = min(diff(log10(unique(fluc.optimal_window_size(:,1)))));
uws2 = min(diff(log10(unique(fluc.optimal_window_size(:,2)))));
jit = (randn(size(fluc.optimal_window_size))).*[uws1 uws2];
etick = [log10([.01 .02 .04 .08 .16 .32])];
eticklabels = arrayfun(@num2str,(10.^etick)*100,'UniformOutput',false);
figure(32)
axc(2)=subplot(3,1,1)
plot(log10(syn_table.syn_efficacy),log10(fluc.optimal_window_size)+jit,'.')
axis tight
set(gca,'YTick',xtick)
set(gca,'YTickLabel',xticklabels)
set(gca,'XTick',etick)
set(gca,'XTickLabel',eticklabels)
xlabel('Efficacy')

box off; set(gca,'TickDir','out')
addLinTrend(log10(syn_table.syn_efficacy),log10(fluc.optimal_window_size(:,1)),'b')
addLinTrend(log10(syn_table.syn_efficacy),log10(fluc.optimal_window_size(:,2)),'r')

axc(3)=subplot(3,1,2)
plot(log10(syn_table.pre_firing_rate),log10(fluc.optimal_window_size)+jit,'.')
xlabel('Pre Firing Rate')
set(gca,'YTick',xtick)
set(gca,'YTickLabel',xticklabels)
box off; set(gca,'TickDir','out')
addLinTrend(log10(syn_table.pre_firing_rate),log10(fluc.optimal_window_size(:,1)),'b')
addLinTrend(log10(syn_table.pre_firing_rate),log10(fluc.optimal_window_size(:,2)),'r')
xlim([-.6 2])

axc(4)=subplot(3,1,3)
plot(log10(syn_table.post_firing_rate),log10(fluc.optimal_window_size)+jit,'.')
xlabel('Post Firing Rate')
ylabel('Optimal Window Size [min]')
set(gca,'YTick',xtick)
set(gca,'YTickLabel',xticklabels)
box off; set(gca,'TickDir','out')
linkaxes(axc,'y')
ylim([log10(60) log10(600)])
addLinTrend(log10(syn_table.post_firing_rate),log10(fluc.optimal_window_size(:,1)),'b')
addLinTrend(log10(syn_table.post_firing_rate),log10(fluc.optimal_window_size(:,2)),'r')
xlim([-.6 2])


%% example of optimal window size for shuffled data

% fast synaptic effect estimated from model
syn_v = fitout.atoms.syn1(exp(fitout.feat_params))>0.01;
syn_t = params.t(syn_v);
sparams = params;
sparams.interval = [min(syn_t) max(syn_t)];

% shuffle randomly (negative control)
Tshuff = shuffle_preserve_efficacy(Tlist,[1 2],sparams);
[opt_res_shuff,smooth_estimate_shuff] = optimizeWindowSize([1 2],Tshuff,params,fitout,20,0);

figure(33)
shadedErrorBar(smooth_estimate.tcenter/60,smooth_estimate.eff,smooth_estimate.effse,'b')
hold on
shadedErrorBar(smooth_estimate_shuff.tcenter/60,smooth_estimate_shuff.eff,smooth_estimate_shuff.effse,'r')
hold off
xlim([0 max(cellfun(@max,Tlist))/60])
line(xlim(),[0 0])
axis tight; box off;
set(gca,'TickDir','out')
xlabel('Time [min]')
ylabel('Efficacy')
tcenter = linspace(0,max(cellfun(@max,Tlist)),11);
rr = nanmedian(smooth_estimate.effse);
yl = ylim();
ylim([max(yl(1),min(smooth_estimate.eff)-rr*2) min(yl(2),max(smooth_estimate.eff)+rr*2)])



%%
figure(34)
opt_res.dev(~isfinite(opt_res.dev))=NaN;
b=nanmean(squeeze(opt_res.dev(1,2,:)));
% errorbar(opt_res.win_vec,mean(squeeze(opt_res.dev(:,2,:)),2)-b,std(squeeze(opt_res.dev(:,2,:)),1,2)/sqrt(5))
plot(opt_res.win_vec,nanmean(squeeze(opt_res.dev(:,2,:)),2)-b)
hold on
b=nanmean(squeeze(opt_res_shuff.dev(1,2,:)));
% errorbar(opt_res_shuff.win_vec,mean(squeeze(opt_res_shuff.dev(:,2,:)),2)-b,std(squeeze(opt_res_shuff.dev(:,2,:)),1,2)/sqrt(5))
plot(opt_res.win_vec,nanmean(squeeze(opt_res_shuff.dev(:,2,:)),2)-b)
hold off
box off
xlabel('Optimal Window Size')
ylabel('Deviance')