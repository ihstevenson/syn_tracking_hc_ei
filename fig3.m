load("results\optimal_fluctuation_windows.mat")


%% CV and optimal window size results for data and controls

prctile(fluc.optimal_window_size,[25 50 75])
prctile(fluc.cv,[25 50 75])

figure(31)
e=linspace(-1.5,3,1024);
[f,x] = ksdensity(log10(fluc.cv(:,1)),e);
shadedErrorBar(x,f*0,[f; f*0]);
hold on
[fcontrol,x] = ksdensity(log10(fluc.cv(:,2)),e);
shadedErrorBar(x,fcontrol*0,[fcontrol; fcontrol*0]);
box off; set(gca,'TickDir','out')
hold off
xlabel('CV')

xtick = [log10([60 120 240 480])];
xticklabels = arrayfun(@num2str,(10.^xtick)/60,'UniformOutput',false);
figure(32)
axc(1)=subplot(4,1,1)
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
set(gca,'XTickLabel','')

axc(2)=subplot(4,1,2)
plot(log10(fluc.optimal_window_size),syn_table.syn_efficacy,'.')
set(gca,'XTick',xtick)
set(gca,'XTickLabel','')
ylabel('Efficacy')
box off; set(gca,'TickDir','out')

axc(3)=subplot(4,1,3)
plot(log10(fluc.optimal_window_size),log10(syn_table.pre_firing_rate),'.')
ylabel('Pre Firing Rate')
set(gca,'XTick',xtick)
set(gca,'XTickLabel','')
box off; set(gca,'TickDir','out')

axc(4)=subplot(4,1,4)
plot(log10(fluc.optimal_window_size),log10(syn_table.post_firing_rate),'.')
ylabel('Post Firing Rate')
xlabel('Optimal Window Size [min]')
set(gca,'XTick',xtick)
set(gca,'XTickLabel',xticklabels)
box off; set(gca,'TickDir','out')
linkaxes(axc,'x')
xlim([log10(60) log10(600)])

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
b=mean(squeeze(opt_res.dev(1,2,:)));
% errorbar(opt_res.win_vec,mean(squeeze(opt_res.dev(:,2,:)),2)-b,std(squeeze(opt_res.dev(:,2,:)),1,2)/sqrt(5))
plot(opt_res.win_vec,mean(squeeze(opt_res.dev(:,2,:)),2)-b)
hold on
b=mean(squeeze(opt_res_shuff.dev(1,2,:)));
% errorbar(opt_res_shuff.win_vec,mean(squeeze(opt_res_shuff.dev(:,2,:)),2)-b,std(squeeze(opt_res_shuff.dev(:,2,:)),1,2)/sqrt(5))
plot(opt_res.win_vec,mean(squeeze(opt_res_shuff.dev(:,2,:)),2)-b)
hold off