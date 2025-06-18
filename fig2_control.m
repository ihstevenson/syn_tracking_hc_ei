res_dir = 'C:\Users\istevenson\Documents\GitHub\fluctuation_project\detection\';

glist = dir([res_dir '*_detected_final.mat']);
rec_num=29
cnx = 26

% Fig 2D examples are...
% 1-13, 35-3
% 51-8, 42-10
% 29-26, 25-25

prefix = strsplit(glist(rec_num).name,'_');
prefix = prefix{1};
load([res_dir prefix '_detected_final'])
load([res_dir prefix '_opt_window'])

[r,c] = find(result.detected_cnx);
fit_res.ssi = result.ssi(sub2ind(size(result.ssi),r,c))';
fit_res.kidx = (fit_res.cvsmoo<.1) & fit_res.ssi<0.2 & (fit_res.rsq>.5) & (fit_res.dur<1) & (fit_res.lat>0.75) & (fit_res.lat<10) & (abs(fit_res.eff)<.5) & (abs(fit_res.eff)>.01) & fit_res.ex;
ijlist = [r(fit_res.kidx) c(fit_res.kidx)];
ij = ijlist(cnx,:);


%% optimize window size for all classified putative synapses

fidx = find(r==ijlist(cnx,1) & c==ijlist(cnx,2));
fsyn = fitout(fidx);
sparams=params;
syn_v = fsyn.atoms.syn1(exp(fsyn.feat_params))>0.01;
syn_t = params.t(syn_v);
sparams.interval = [min(syn_t) max(syn_t)];

% Shuffle randomly (negative control)
Tshuff = shuffle_preserve_efficacy(Tlist,ij,sparams);

% Shuffle weighted (positive control)
%get baseline fraction
pfrac = 1-sum(fsyn.lamsmoo(syn_v))/sum(fsyn.yhat(syn_v));
Tshuff = shuffle_manip_efficacy(Tlist,ij,sparams,pfrac,1/(60*40));

params.nknots=6;
params.doPlot = false;

figure(2)
[opt_res_shuff(cnx),smooth_estimate_shuff(cnx)] = optimizeWindowSize(ijlist(cnx,:),Tshuff,params,fsyn,20,1);


%% original
figure(31)

subplot(4,10,2:10)
shadedErrorBar(smooth_estimate(cnx).tcenter/60,smooth_estimate(cnx).eff,smooth_estimate(cnx).effse)
xlim([0 max(cellfun(@max,Tlist))/60])
line(xlim(),[0 0])
axis tight; box off;
set(gca,'TickDir','out')
xlabel('Time [min]')
ylabel('Efficacy')
tcenter = linspace(0,max(cellfun(@max,Tlist)),11);
hold on
plot(tcenter/60,interp1(smooth_estimate(cnx).tcenter,smooth_estimate(cnx).eff,tcenter)-interp1(smooth_estimate(cnx).tcenter,smooth_estimate(cnx).effse,tcenter),'^')
hold off
rr = nanmedian(smooth_estimate(cnx).effse);
yl = ylim();
ylim([max(yl(1),min(smooth_estimate(cnx).eff)-rr*2) min(yl(2),max(smooth_estimate(cnx).eff)+rr*2)])

fidx=find(r==ijlist(cnx,1) & c==ijlist(cnx,2));
fitout(fidx)
smooth_estimate(cnx)

title(strjoin([prefix ' ' num2str(tab_sub.unit_id(r(fidx))) ' ' num2str(tab_sub.unit_id(c(fidx))) tab_sub.ecephys_structure_acronym(r(fidx)) tab_sub.ecephys_structure_acronym(c(fidx))]))

%
% figure(32)
cc = spk_xcorr(Tlist{ijlist(cnx,1)},Tlist{ijlist(cnx,2)},params.interval(1),params.interval(2),params.num_bins);
ac = spk_xcorr(Tlist{ijlist(cnx,1)},Tlist{ijlist(cnx,1)},params.interval(1),params.interval(2),params.num_bins);
ac(params.t==0)=NaN;
acp = spk_xcorr(Tlist{ijlist(cnx,2)},Tlist{ijlist(cnx,2)},params.interval(1),params.interval(2),params.num_bins);
acp(params.t==0)=NaN;
% clf
subplot(4,10,11)
bar(params.t,cc,1,'EdgeColor','none')
hold on
plot(params.t,fitout(fidx).yhat*length(Tlist{ijlist(cnx,1)}))
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

ww = smooth_estimate(cnx).window;
ww = 10*60;
cc(:,2:10)=NaN;
ac(:,2:10)=NaN;
acp(:,2:10)=NaN;

for i=2:10
    tid = (Tlist{ijlist(cnx,1)}>(tcenter(i)-ww/2) & Tlist{ijlist(cnx,1)}<(tcenter(i)+ww/2));
    if sum(tid)>0
        cc(:,i) = spk_xcorr(Tlist{ijlist(cnx,1)}(tid),Tlist{ijlist(cnx,2)},params.interval(1),params.interval(2),params.num_bins);
        ac(:,i) = spk_xcorr(Tlist{ijlist(cnx,1)}(tid),Tlist{ijlist(cnx,1)},params.interval(1),params.interval(2),params.num_bins);
        ac(params.t==0,i)=NaN;
        cc(:,i) = cc(:,i)/sum(tid);
        ac(:,i) = ac(:,i)/sum(tid);
    end
    tidp = (Tlist{ijlist(cnx,2)}>(tcenter(i)-ww/2) & Tlist{ijlist(cnx,2)}<(tcenter(i)+ww/2));
    if sum(tidp)>0
        acp(:,i) = spk_xcorr(Tlist{ijlist(cnx,2)}(tidp),Tlist{ijlist(cnx,2)},params.interval(1),params.interval(2),params.num_bins);
        acp(params.t==0,i)=NaN;
        acp(:,i) = acp(:,i)/sum(tidp);
    end
end

for i=2:10
    subplot(4,10,10+i)
    bar(params.t,cc(:,i),1,'EdgeColor','none')
    box off; set(gca,'TickDir','out')
    %                 line([0 0],ylim())
    ylim([0 max(max(cc(:,2:end)))])
    xlim([-.01 .01])
    subplot(4,10,20+i)
    bar(params.t,ac(:,i),1,'EdgeColor','none')
    box off; set(gca,'TickDir','out')
    %                 line([0 0],ylim())
    ylim([0 max(max(ac(:,2:end)))])
    subplot(4,10,30+i)
    bar(params.t,acp(:,i),1,'EdgeColor','none')
    box off; set(gca,'TickDir','out')
    %                 line([0 0],ylim())
    ylim([0 max(max(acp(:,2:end)))])
end

%%

figure(32)

subplot(4,10,2:10)
shadedErrorBar(smooth_estimate_shuff(cnx).tcenter/60,smooth_estimate_shuff(cnx).eff,smooth_estimate_shuff(cnx).effse)
xlim([0 max(cellfun(@max,Tlist))/60])
line(xlim(),[0 0])
axis tight; box off;
set(gca,'TickDir','out')
xlabel('Time [min]')
ylabel('Efficacy')
tcenter = linspace(0,max(cellfun(@max,Tlist)),11);
hold on
plot(tcenter/60,interp1(smooth_estimate_shuff(cnx).tcenter,smooth_estimate_shuff(cnx).eff,tcenter)-interp1(smooth_estimate_shuff(cnx).tcenter,smooth_estimate_shuff(cnx).effse,tcenter),'^')
hold off
rr = nanmedian(smooth_estimate(cnx).effse);
%yl = ylim();
ylim([max(yl(1),min(smooth_estimate(cnx).eff)-rr*2) min(yl(2),max(smooth_estimate(cnx).eff)+rr*2)])

fidx=find(r==ijlist(cnx,1) & c==ijlist(cnx,2));
fitout(fidx)
smooth_estimate(cnx)

title(strjoin([prefix ' ' num2str(tab_sub.unit_id(r(fidx))) ' ' num2str(tab_sub.unit_id(c(fidx))) tab_sub.ecephys_structure_acronym(r(fidx)) tab_sub.ecephys_structure_acronym(c(fidx))]))

%
% figure(32)
cc = spk_xcorr(Tshuff{ijlist(cnx,1)},Tshuff{ijlist(cnx,2)},params.interval(1),params.interval(2),params.num_bins);
ac = spk_xcorr(Tshuff{ijlist(cnx,1)},Tshuff{ijlist(cnx,1)},params.interval(1),params.interval(2),params.num_bins);
ac(params.t==0)=NaN;
acp = spk_xcorr(Tlist{ijlist(cnx,2)},Tlist{ijlist(cnx,2)},params.interval(1),params.interval(2),params.num_bins);
acp(params.t==0)=NaN;
% clf
subplot(4,10,11)
bar(params.t,cc,1,'EdgeColor','none')
hold on
plot(params.t,fitout(fidx).yhat*length(Tlist{ijlist(cnx,1)}))
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

ww = smooth_estimate_shuff(cnx).window;
ww = 10*60;
cc(:,2:10)=NaN;
ac(:,2:10)=NaN;
acp(:,2:10)=NaN;

for i=2:10
    tid = (Tshuff{ijlist(cnx,1)}>(tcenter(i)-ww/2) & Tshuff{ijlist(cnx,1)}<(tcenter(i)+ww/2));
    if sum(tid)>0
        cc(:,i) = spk_xcorr(Tshuff{ijlist(cnx,1)}(tid),Tshuff{ijlist(cnx,2)},params.interval(1),params.interval(2),params.num_bins);
        ac(:,i) = spk_xcorr(Tshuff{ijlist(cnx,1)}(tid),Tshuff{ijlist(cnx,1)},params.interval(1),params.interval(2),params.num_bins);
        ac(params.t==0,i)=NaN;
        cc(:,i) = cc(:,i)/sum(tid);
        ac(:,i) = ac(:,i)/sum(tid);
    end
    tidp = (Tshuff{ijlist(cnx,2)}>(tcenter(i)-ww/2) & Tshuff{ijlist(cnx,2)}<(tcenter(i)+ww/2));
    if sum(tidp)>0
        acp(:,i) = spk_xcorr(Tshuff{ijlist(cnx,2)}(tidp),Tshuff{ijlist(cnx,2)},params.interval(1),params.interval(2),params.num_bins);
        acp(params.t==0,i)=NaN;
        acp(:,i) = acp(:,i)/sum(tidp);
    end
end

for i=2:10
    subplot(4,10,10+i)
    bar(params.t,cc(:,i),1,'EdgeColor','none')
    box off; set(gca,'TickDir','out')
    %                 line([0 0],ylim())
    ylim([0 max(max(cc(:,2:end)))])
    xlim([-.01 .01])
    subplot(4,10,20+i)
    bar(params.t,ac(:,i),1,'EdgeColor','none')
    box off; set(gca,'TickDir','out')
    %                 line([0 0],ylim())
    ylim([0 max(max(ac(:,2:end)))])
    subplot(4,10,30+i)
    bar(params.t,acp(:,i),1,'EdgeColor','none')
    box off; set(gca,'TickDir','out')
    %                 line([0 0],ylim())
    ylim([0 max(max(acp(:,2:end)))])
end

%% get results for all pairs...

glist = dir([res_dir '*_detected_final.mat']);
for rec_num=23:length(glist)
    prefix = strsplit(glist(rec_num).name,'_');
    prefix = prefix{1};
    load([res_dir prefix '_detected_final'])
    load([res_dir prefix '_opt_window'])
    [r,c] = find(result.detected_cnx);
    fit_res.ssi = result.ssi(sub2ind(size(result.ssi),r,c))';
    fit_res.kidx = (fit_res.cvsmoo<.1) & fit_res.ssi<0.2 & (fit_res.rsq>.5) & (fit_res.dur<1) & (fit_res.lat>0.75) & (fit_res.lat<10) & (abs(fit_res.eff)<.5) & (abs(fit_res.eff)>.01) & fit_res.ex;
    ijlist = [r(fit_res.kidx) c(fit_res.kidx)];

    clear opt_res_shuff smooth_estimate_shuff
    if ~isempty(ijlist)

        for cnx = 1:size(ijlist,1)

            ij = ijlist(cnx,:);

            fidx = find(r==ijlist(cnx,1) & c==ijlist(cnx,2));
            fsyn = fitout(fidx);
            sparams=params;
            syn_v = fsyn.atoms.syn1(exp(fsyn.feat_params))>0.01;
            syn_t = params.t(syn_v);
            sparams.interval = [min(syn_t) max(syn_t)];

            % Shuffle randomly (negative control)
            % Tshuff = shuffle_preserve_efficacy(Tlist,ij,sparams);

            % Shuffle weighted (positive control)
            %get baseline fraction
            pfrac = 1-sum(fsyn.lamsmoo(syn_v))/sum(fsyn.yhat(syn_v));
            Tshuff = shuffle_manip_efficacy(Tlist,ij,sparams,pfrac,1/(60*40));

            params.nknots=6;
            params.doPlot = false;

            figure(2)
            [opt_res_shuff(cnx),smooth_estimate_shuff(cnx)] = optimizeWindowSize(ijlist(cnx,:),Tshuff,params,fsyn,20,1);
        end
        save(['E:\fluctuation_control_pos\' prefix '_opt_window_shuff.mat'],'opt_res_shuff','smooth_estimate_shuff')
    end
end


%% summarize results for all pairs...

optimal_window_size=[];
fluc_cv=[];
contrib=[];
npost=[];
T=[];
cc=1;

glist = dir([res_dir '*_detected_final.mat']);
for rec_num=1:length(glist)
    prefix = strsplit(glist(rec_num).name,'_');
    prefix = prefix{1};
    load([res_dir prefix '_detected_final'])
    [r,c] = find(result.detected_cnx);
    fit_res.ssi = result.ssi(sub2ind(size(result.ssi),r,c))';
    fit_res.kidx = (fit_res.cvsmoo<.1) & fit_res.ssi<0.2 & (fit_res.rsq>.5) & (fit_res.dur<1) & (fit_res.lat>0.75) & (fit_res.lat<10) & (abs(fit_res.eff)<.5) & (abs(fit_res.eff)>.01) & fit_res.ex;
    ijlist = [r(fit_res.kidx) c(fit_res.kidx)];

    if ~isempty(ijlist)
        load([res_dir prefix '_opt_window'])
        load(['E:\fluctuation_control\' prefix '_opt_window_shuff.mat'])
        for cnx = 1:size(ijlist,1)
            fidx = find(r==ijlist(cnx,1) & c==ijlist(cnx,2));
            fsyn = fitout(fidx);

            optimal_window_size(cc,1) = smooth_estimate(cnx).window;
            fluc_cv(cc,1) = nanstd(smooth_estimate(cnx).eff)/fsyn.eff;

            optimal_window_size(cc,2) = smooth_estimate_shuff(cnx).window;
            fluc_cv(cc,2) = nanstd(smooth_estimate_shuff(cnx).eff)/fsyn.eff;

            contrib(cc,1) = length(Tlist{ijlist(cnx,1)})*fsyn.eff/length(Tlist{ijlist(cnx,2)});
            npost(cc,1) = length(Tlist{ijlist(cnx,2)});
            T(cc,1) = max(cellfun(@max,Tlist));
            cc=cc+1;
        end
    end
end

%%

figure(101)
e=linspace(-1.5,3,100);
ax(1)=subplot(2,2,1)
histogram(log10(fluc_cv(:,1)),e,'EdgeColor','none')
box off; set(gca,'TickDir','out')
ax(2)=subplot(2,2,3)
histogram(log10(fluc_cv(:,2)),e,'EdgeColor','none')
box off; set(gca,'TickDir','out')
linkaxes(ax)

e=linspace(0,600,10);
axc(1)=subplot(2,2,2)
histogram(optimal_window_size(:,1),e,'EdgeColor','none')
box off; set(gca,'TickDir','out')
axc(2)=subplot(2,2,4)
histogram(optimal_window_size(:,2),e,'EdgeColor','none')
box off; set(gca,'TickDir','out')
linkaxes(axc)

prctile(optimal_window_size,[25 50 75])
prctile(fluc_cv,[25 50 75])