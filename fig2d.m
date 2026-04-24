
load('results\subset_prepost_spike_times.mat')
load('results\syn_index')
%% 

% Data plotted in Fig 2 comes from pair_ids 1-6
% corresponding to pairs in syn_table(subset_index,:)
% running optimal window size estimation takes ~10min
pair_id = 5; % from data subset
Tlist(1) = Tlist_pre(pair_id);
Tlist(2) = Tlist_post(pair_id);

% convolutional CCG model parameters
params=[];
params.interval_centers = [-.02 .02];
params.num_bins = 201;
params.dt = range(params.interval_centers)/(params.num_bins-1);
params.interval = params.interval_centers+[-1 1]*params.dt/2;
params.t = linspace(params.interval(1),params.interval(2),params.num_bins+1);
params.t = params.t(1:end-1)+mean(diff(params.t))/2;

params.mask = ones(size(params.t));
params.mask(abs(params.t)<0.001)=0;
params.mask = params.mask>0;
params.nknots = 6;
params.asym=false;
params.doPlot=false;

% fit overall model
fitout = composite_model_fit([1 2],Tlist,params,'i>j');

cc = spk_xcorr(Tlist{1},Tlist{2},params.interval(1),params.interval(2),params.num_bins);
ac = spk_xcorr(Tlist{1},Tlist{1},params.interval(1),params.interval(2),params.num_bins);
ac(params.t==0)=NaN;
acp = spk_xcorr(Tlist{2},Tlist{2},params.interval(1),params.interval(2),params.num_bins);
acp(params.t==0)=NaN;

figure(24)
clf
subplot(4,10,11)
bar(params.t,cc,1,'EdgeColor','none')
hold on
plot(params.t,fitout.yhat*length(Tlist{1}))
hold off
box off; set(gca,'TickDir','out')
line([0 0],ylim())
subplot(4,10,20+1)
bar(params.t,ac,1,'EdgeColor','none')
box off; set(gca,'TickDir','out')
line([0 0],ylim())
subplot(4,10,30+1)
bar(params.t,acp,1,'EdgeColor','none')
box off; set(gca,'TickDir','out')
line([0 0],ylim())



% run optimal window size estimation (default is k=5 cross-validation folds)
[opt_res,smooth_estimate] = optimizeWindowSize([1 2],Tlist,params,fitout,20,0);


subplot(4,10,2:10)
shadedErrorBar(smooth_estimate.tcenter/60,smooth_estimate.eff,smooth_estimate.effse)
xlim([0 max(cellfun(@max,Tlist))/60])
line(xlim(),[0 0])
axis tight; box off;
set(gca,'TickDir','out')
xlabel('Time [min]')
ylabel('Efficacy')
tcenter = linspace(0,max(cellfun(@max,Tlist)),11);
hold on
plot(tcenter/60,interp1(smooth_estimate.tcenter,smooth_estimate.eff,tcenter)-interp1(smooth_estimate.tcenter,smooth_estimate.effse,tcenter),'^')
hold off
rr = nanmedian(smooth_estimate.effse);
yl = ylim();
ylim([max(yl(1),min(smooth_estimate.eff)-rr*2) min(yl(2),max(smooth_estimate.eff)+rr*2)])


% Get local windows showing CCG and ACGs
ww = 10*60; % window for plotting local ccgs and acgs
cc(:,2:10)=NaN;
ac(:,2:10)=NaN;
acp(:,2:10)=NaN;

for i=2:10
    tid = (Tlist{1}>(tcenter(i)-ww/2) & Tlist{1}<(tcenter(i)+ww/2));
    if sum(tid)>0
        cc(:,i) = spk_xcorr(Tlist{1}(tid),Tlist{2},params.interval(1),params.interval(2),params.num_bins);
        ac(:,i) = spk_xcorr(Tlist{1}(tid),Tlist{1},params.interval(1),params.interval(2),params.num_bins);
        ac(params.t==0,i)=NaN;
        cc(:,i) = cc(:,i)/sum(tid);
        ac(:,i) = ac(:,i)/sum(tid);
    end
    tidp = (Tlist{2}>(tcenter(i)-ww/2) & Tlist{2}<(tcenter(i)+ww/2));
    if sum(tidp)>0
        acp(:,i) = spk_xcorr(Tlist{2}(tidp),Tlist{2},params.interval(1),params.interval(2),params.num_bins);
        acp(params.t==0,i)=NaN;
        acp(:,i) = acp(:,i)/sum(tidp);
    end
end

for i=2:10
    subplot(4,10,10+i)
    bar(params.t,cc(:,i),1,'EdgeColor','none')
    box off; set(gca,'TickDir','out')
    ylim([0 max(max(cc(:,2:end)))])
    xlim([-.01 .01])
    subplot(4,10,20+i)
    bar(params.t,ac(:,i),1,'EdgeColor','none')
    box off; set(gca,'TickDir','out')
    ylim([0 max(max(ac(:,2:end)))])
    subplot(4,10,30+i)
    bar(params.t,acp(:,i),1,'EdgeColor','none')
    box off; set(gca,'TickDir','out')
    ylim([0 max(max(acp(:,2:end)))])
end
