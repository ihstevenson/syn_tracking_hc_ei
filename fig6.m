
sr = 1.25e3;
bpFilt = designfilt('bandpassfir','FilterOrder',20, ...
    'CutoffFrequency1',100,'CutoffFrequency2',250, ...
    'SampleRate',sr);
ripple_params.zscore_cutoff = 2.5;
ripple_params.zscore_edge = 2;
ripple_params.minimum_duration = 15; % ms
ripple_params.minimum_interval = 15; % ms


ccc=1;
ripple_all=cell(0);
for rec_num=1:length(glist)
    prefix = strsplit(glist(rec_num).name,'_');
    prefix = prefix{1};
    fprintf('%s...\n',prefix)

    load([res_dir prefix '_detected_final'])

    [r,c] = find(result.detected_cnx);
    fit_res.ssi = result.ssi(sub2ind(size(result.ssi),r,c))';
    fit_res.kidx = (fit_res.cvsmoo<.1) & fit_res.ssi<0.2 & (fit_res.rsq>.5) & (fit_res.dur<1) & (fit_res.lat>0.75) & (fit_res.lat<10) & (abs(fit_res.eff)<.5) & (abs(fit_res.eff)>.01) & fit_res.ex;
    ijlist = [r(fit_res.kidx) c(fit_res.kidx)];
    T = max(cellfun(@max,Tlist));

    for cnx=1:size(ijlist,1)
        fprintf('%02i...\n',cnx)
        pre = ijlist(cnx,1);
        post = ijlist(cnx,2);
        fidx = find(r==pre & c==post);

        probe_id = tab_sub.probe_id(pre);
        if isfile([lfp_dir prefix '_' num2str(probe_id) '_lfp.mat'])
            load([lfp_dir prefix '_' num2str(probe_id) '_lfp.mat'])
        else
            ripples = struct();
            ripple_all{ccc} = ripples;
            fprintf('no LFP')
            ccc=ccc+1;
            continue
        end

        chans = double(chans);
        chan = tab_sub.peak_channel_id(pre);
        [~,chan_index]=min(abs(chans-chan));
        chan = chans(chan_index);

        lfp = double(lfp(:,chan_index));
        lfp(isnan(lfp))=0;
        lfp_filt = filtfilt(bpFilt,lfp);

        [ripples,zscore] = ripple_detection_v2(lfp_filt,sr,time,ripple_params);
        ripple_all{ccc} = ripples;
        ccc=ccc+1;
    end
end

%% visualize single example with ripples

addpath('C:\Users\istevenson\Dropbox\glm_dist\')
addpath('C:\Users\istevenson\Dropbox\__matlab_toolboxes\own_stuff')
rec_num=19;
prefix = strsplit(glist(rec_num).name,'_');
prefix = prefix{1};
fprintf('%s...\n',prefix)

load([res_dir prefix '_detected_final'])
cnx_id = cnx_list(:,1)==rec_num;
% rate_estimate_sub = rate_estimate_all(cnx_id);
% clear rate_estimate_all

[r,c] = find(result.detected_cnx);
fit_res.ssi = result.ssi(sub2ind(size(result.ssi),r,c))';
fit_res.kidx = (fit_res.cvsmoo<.1) & fit_res.ssi<0.2 & (fit_res.rsq>.5) & (fit_res.dur<1) & (fit_res.lat>0.75) & (fit_res.lat<10) & (abs(fit_res.eff)<.5) & (abs(fit_res.eff)>.01) & fit_res.ex;
ijlist = [r(fit_res.kidx) c(fit_res.kidx)];
T = max(cellfun(@max,Tlist));

unique_cells = unique([ijlist(:,1); ijlist(:,2)]);
[~,sidx]=sort(tab_sub{unique_cells,'ecephys_structure_acronym'});
unique_cells = unique_cells(sidx);
tab_sub{unique_cells,'ecephys_structure_acronym'}
[~,~,ia]=unique(tab_sub{unique_cells,'ecephys_structure_acronym'});

probet = tabulate(categorical(tab_sub.probe_id(unique_cells)));
[~,i]=max(cellfun(@(x) x,probet(:,2)));
probe_id=probet{i,1};
% probe_id=probet{2,1};
load([lfp_dir prefix '_' probe_id '_lfp.mat'])

% limit to cells on same probe
same_probe = tab_sub{unique_cells,'probe_id'}==str2num(probe_id);
[~,~,ia]=unique(tab_sub{unique_cells(same_probe),'ecephys_structure_acronym'});

neuron_subset = unique_cells(same_probe);
closest_lfp_id=[];
ispre=[];
for i=1:length(neuron_subset)
    [~,closest_lfp_id(i)]=min(abs(chans-tab_sub{neuron_subset(i),'peak_channel_id'}));
    ispre(i) = any(ijlist(:,1)==neuron_subset(i));
end
[~,ispre_srt]=sort(ispre);

%%
unique(closest_lfp_id(ispre>0))
cnxf = find(cnx_id);


rip_num = 390;
rwin=.4;
trange = [ripple_all{cnxf(1)}(rip_num).start_time-rwin ripple_all{cnxf(1)}(rip_num).end_time+rwin];

figure(1)
clf
subplot(3,1,1)
tid = time>trange(1) & time<trange(2);
plot(time(tid),lfp(tid,18))
xlim(trange)
ylim([-1 1]*max(abs(ylim())))
ylim([-2 2]*8e-4)
% line([trange(1)+rwin trange(2)-rwin],[1 1]*.9*max(ylim()))

st1=arrayfun(@(x) x.start_time,ripple_all{cnxf(1)});
st2=arrayfun(@(x) x.end_time,ripple_all{cnxf(1)});
multirips = st2>trange(1) & st1<trange(2);
line([st1(multirips)' st2(multirips)']',repmat([1 1],sum(multirips),1)'*.9*max(ylim()),'Color','k')


box off; set(gca,'TickDir','out')
subplot(3,1,2)
% drawRaster(Tlist(unique_cells(same_probe)),1,ia,[trange(1) trange(2)])
drawRaster(Tlist(neuron_subset(ispre_srt)),1,ispre(ispre_srt),[trange(1) trange(2)])
box off; set(gca,'TickDir','out')

%% xcorr - conditioning on ripples

addpath('C:\Users\istevenson\Documents\GitHub\syngalong_v0-1\code')

fparams=[];
fparams.trange = [-0.025 0.025];
fparams.tn = 102;
fparams.t = linspace(fparams.trange(1),fparams.trange(2),fparams.tn);
fparams.t = fparams.t(1:end-1)+mean(diff(fparams.t))/2;

sparams=[];
sparams.trange = [-0.25 0.25];
sparams.tn = 102;
sparams.t = linspace(sparams.trange(1),sparams.trange(2),sparams.tn);
sparams.t = sparams.t(1:end-1)+mean(diff(sparams.t))/2;


%%
cnx=28;
cnx_id = cnx_list(:,1)==rec_num;
cnxf = find(cnx_id);

ripple_starts = arrayfun(@(x)x.start_time,ripple_all{cnxf(cnx)});
ripple_stops = arrayfun(@(x)x.end_time,ripple_all{cnxf(cnx)});

[c,deltaT] = spk_xcorr(Tlist{ijlist(cnx,1)},Tlist{ijlist(cnx,2)},fparams.trange(1),fparams.trange(2),fparams.tn-1);
subplot(2,3,1)
bar(params.t,c,1,'EdgeColor','none')
box off; set(gca,'TickDir','out')


[ra,deltaT] = spk_xcorr(ripple_starts,ripple_starts,sparams.trange(1),sparams.trange(2),sparams.tn-1);
[sc1,deltaT1] = spk_xcorr(ripple_starts,Tlist{ijlist(cnx,1)},sparams.trange(1),sparams.trange(2),sparams.tn-1);
[sc2,deltaT2] = spk_xcorr(ripple_starts,Tlist{ijlist(cnx,2)},sparams.trange(1),sparams.trange(2),sparams.tn-1);
ra(sparams.t==0)=NaN;

subplot(2,3,2)
bar(sparams.t,sc1,1,'EdgeColor','none')
box off; set(gca,'TickDir','out')
subplot(2,3,3)
bar(sparams.t,sc2,1,'EdgeColor','none')
box off; set(gca,'TickDir','out')
subplot(2,3,6)
bar(sparams.t,ra,1,'EdgeColor','none')
box off; set(gca,'TickDir','out')


%t1_ripp = unique(deltaT1(:,3));
t1_ripp = unique(deltaT1(deltaT1(:,1)>0 & deltaT1(:,1)<ripple_stops(deltaT1(:,2))',3));
t2_ripp = unique(deltaT2(deltaT2(:,1)>0 & deltaT2(:,1)<ripple_stops(deltaT2(:,2))',3));
[c1,deltaT] = spk_xcorr(Tlist{ijlist(cnx,1)}(t1_ripp),Tlist{ijlist(cnx,2)},fparams.trange(1),fparams.trange(2),fparams.tn-1);
subplot(2,3,4)
bar(params.t,c1,1,'EdgeColor','none')
box off; set(gca,'TickDir','out')

t1_noripp = setdiff(1:length(Tlist{ijlist(cnx,1)}),t1_ripp);
[c2,deltaT] = spk_xcorr(Tlist{ijlist(cnx,1)}(t1_noripp),Tlist{ijlist(cnx,2)},fparams.trange(1),fparams.trange(2),fparams.tn-1);
subplot(2,3,5)
bar(params.t,c2,1,'EdgeColor','none')
box off; set(gca,'TickDir','out')


%% collect ripple results over all pairs

cnxc=1;
ripp_pair_stats=[];
ripple_act_pre=[];
ripple_act_post=[];
ripple_act_auto=[];
xcorr_ripple = [];
xcorr_noripple = [];

for rec_num=1:length(glist)
    prefix = strsplit(glist(rec_num).name,'_');
    prefix = prefix{1};
    fprintf('%s...\n',prefix)

    load([res_dir prefix '_detected_final'])


    [r,c] = find(result.detected_cnx);
    fit_res.ssi = result.ssi(sub2ind(size(result.ssi),r,c))';
    fit_res.kidx = (fit_res.cvsmoo<.1) & fit_res.ssi<0.2 & (fit_res.rsq>.5) & (fit_res.dur<1) & (fit_res.lat>0.75) & (fit_res.lat<10) & (abs(fit_res.eff)<.5) & (abs(fit_res.eff)>.01) & fit_res.ex;
    ijlist = [r(fit_res.kidx) c(fit_res.kidx)];
    T = max(cellfun(@max,Tlist));

    cnx_id = cnx_list(:,1)==rec_num;
    cnxf = find(cnx_id);


    for cnx=1:size(ijlist,1)
        fprintf('%02i...\n',cnx)

        if ~isfield(ripple_all{cnxf(cnx)},'start_time')
            ripple_starts=[];
            ripple_stops=[];
        else
            ripple_starts = arrayfun(@(x)x.start_time,ripple_all{cnxf(cnx)});
            ripple_stops = arrayfun(@(x)x.end_time,ripple_all{cnxf(cnx)});
        end
        ripp_len = ripple_stops-ripple_starts;

        if ~isempty(ripple_starts)
            [c,deltaT] = spk_xcorr(Tlist{ijlist(cnx,1)},Tlist{ijlist(cnx,2)},fparams.trange(1),fparams.trange(2),fparams.tn-1);
            [ra,deltaT] = spk_xcorr(ripple_starts,ripple_starts,sparams.trange(1),sparams.trange(2),sparams.tn-1);
            [sc1,deltaT1] = spk_xcorr(ripple_starts,Tlist{ijlist(cnx,1)},sparams.trange(1),sparams.trange(2),sparams.tn-1);
            [sc2,deltaT2] = spk_xcorr(ripple_starts,Tlist{ijlist(cnx,2)},sparams.trange(1),sparams.trange(2),sparams.tn-1);
            t1_ripp = unique(deltaT1(deltaT1(:,1)>0 & deltaT1(:,1)<ripp_len(deltaT1(:,2))',3));
            t2_ripp = unique(deltaT2(deltaT2(:,1)>0 & deltaT2(:,1)<ripp_len(deltaT2(:,2))',3));
            [c1,deltaT] = spk_xcorr(Tlist{ijlist(cnx,1)}(t1_ripp),Tlist{ijlist(cnx,2)},fparams.trange(1),fparams.trange(2),fparams.tn-1);
            t1_noripp = setdiff(1:length(Tlist{ijlist(cnx,1)}),t1_ripp);
            [c2,deltaT] = spk_xcorr(Tlist{ijlist(cnx,1)}(t1_noripp),Tlist{ijlist(cnx,2)},fparams.trange(1),fparams.trange(2),fparams.tn-1);

            % ripple activity and conditional pair xcorrs
            ripple_act_pre(cnxc,:)=sc1';
            ripple_act_post(cnxc,:)=sc2';
            ripple_act_auto(cnxc,:)=ra';
            xcorr_ripple(cnxc,:) = c1';
            xcorr_noripple(cnxc,:) = c2';
        end
        % summary stats
        T=max(cellfun(@max,Tlist));
        rippT = sum(ripple_stops-ripple_starts);

        ripp_pair_stats(cnxc,1)=length(ripple_starts); % number of ripples
        ripp_pair_stats(cnxc,2)=mean(ripple_stops-ripple_starts); % ripple length
        ripp_pair_stats(cnxc,3)=length(t1_ripp)/length(Tlist{ijlist(cnx,1)}); % fraction of pre spike in ripples
        ripp_pair_stats(cnxc,4)=length(t1_ripp)/rippT; % pre-rate during ripples
        ripp_pair_stats(cnxc,5)=(length(Tlist{ijlist(cnx,1)})-length(t1_ripp))/(T-rippT);% pre-rate not during ripples
        ripp_pair_stats(cnxc,6)=length(t2_ripp)/rippT; % post-rate during ripples
        ripp_pair_stats(cnxc,7)=(length(Tlist{ijlist(cnx,2)})-length(t2_ripp))/(T-rippT); % post-rate not during ripples
        ripp_pair_stats(cnxc,8)=T; % recording time
        cnxc=cnxc+1;
    end
end

%% plot ripple bulk rate statistics

ripp_pair_stats(ripp_pair_stats(:,1)==0,1)=NaN;
[nanmean(ripp_pair_stats(:,1)) nanstd(ripp_pair_stats(:,1))] % total ripples
[nanmean(ripp_pair_stats(:,1)./ripp_pair_stats(:,8)) nanstd(ripp_pair_stats(:,1)./ripp_pair_stats(:,8))]*60 % ripples per min
[mean(ripp_pair_stats(ripp_pair_stats(:,1)>0,2))] % avg ripp duration



inclus = ripp_pair_stats(:,1)>100 & isfinite(ripp_pair_stats(:,1));
figure(3)
clf
subplot(1,2,1)
loglog(ripp_pair_stats(inclus,5),ripp_pair_stats(inclus,4),'.')
box off; set(gca,'TickDir','out')
line(xlim(),xlim())
% subplot(1,3,2)
hold on
loglog(ripp_pair_stats(inclus,7),ripp_pair_stats(inclus,6),'.')
hold off
box off; set(gca,'TickDir','out')
xlabel('Rate not during ripples [Hz]')
ylabel('Rate during ripples [Hz]')
% line(xlim(),xlim())
subplot(1,2,2)
% no strong correlation between rate pre- and post- changes
loglog(ripp_pair_stats(inclus,4)./ripp_pair_stats(inclus,5),ripp_pair_stats(inclus,6)./ripp_pair_stats(inclus,7),'.')
xlabel('Pre- Ripple Change [%]')
ylabel('Post- Ripple Change [%]')
box off; set(gca,'TickDir','out')
[r,p]=corrcoef(ripp_pair_stats(inclus,4)./ripp_pair_stats(inclus,5),ripp_pair_stats(inclus,6)./ripp_pair_stats(inclus,7))
line([1 1],ylim())
line(xlim(),[1 1])

% rate changes during ripples are similar on average for pre and post
mean(ripp_pair_stats(inclus,4)./ripp_pair_stats(inclus,5))
mean(ripp_pair_stats(inclus,6)./ripp_pair_stats(inclus,7))

% percent of units that decrease with ripple
sum((ripp_pair_stats(inclus,4)./ripp_pair_stats(inclus,5))<1)/sum(inclus)
sum((ripp_pair_stats(inclus,6)./ripp_pair_stats(inclus,7))<1)/sum(inclus)

%% plot ripple-locked rates
%
% ripple_act_pre=[];
% ripple_act_post=[];
% ripple_act_auto=[];
% xcorr_ripple = [];
% xcorr_noripple = [];

figure(4)
% subplot(1,2,1)
% imagesc(sparams.t,1:sum(inclus),ripple_act_pre(inclus,:))
% subplot(1,2,2)
% imagesc(sparams.t,1:sum(inclus),ripple_act_post(inclus,:))


r1n = ripple_act_pre(inclus,:)./ripp_pair_stats(inclus,1)/mean(diff(sparams.t));
r2n = ripple_act_post(inclus,:)./ripp_pair_stats(inclus,1)/mean(diff(sparams.t));

% [~,sidx] = sort(sqrt(ripp_pair_stats(inclus,5).*ripp_pair_stats(inclus,7)));
% [~,sidx] = sort(ripp_pair_stats(inclus,7));
[~,sidx] = sort(sqrt(ripp_pair_stats(inclus,4).*ripp_pair_stats(inclus,6)));
r1n=r1n(sidx,:);
r2n=r2n(sidx,:);
[~,sidx2] = sort(sqrt(ripp_pair_stats(:,4).*ripp_pair_stats(:,6)));

cl = max(max(r1n(:)),max(r2n(:)));
subplot(1,2,1)
imagesc(sparams.t,1:sum(inclus),sqrt(r1n))
caxis([0 sqrt(cl)])
colorbar
line([0 0],ylim())
box off; set(gca,'TickDir','out')
xlabel('Time following ripple onset [s]')
subplot(1,2,2)
imagesc(sparams.t,1:sum(inclus),sqrt(r2n))
caxis([0 sqrt(cl)])
colorbar
line([0 0],ylim())
box off; set(gca,'TickDir','out')
xlabel('Time following ripple onset [s]')


r1n = ripple_act_pre(inclus,:)./ripp_pair_stats(inclus,1)/mean(diff(sparams.t));
r2n = ripple_act_post(inclus,:)./ripp_pair_stats(inclus,1)/mean(diff(sparams.t));
figure(5)
subplot(2,1,1)
plot(sparams.t,mean(r1n))
hold on
plot(sparams.t,mean(r2n))
hold off
line([0 0],ylim())
box off; set(gca,'TickDir','out')
xlabel('Time following ripple onset [s]')
ylabel('Rate [Hz]')
xlim([-.25 .25])
subplot(2,1,2)
plot(sparams.t,mean(r1n./ripp_pair_stats(inclus,5)))
hold on
plot(sparams.t,mean(r2n./ripp_pair_stats(inclus,7)))
hold off
line([0 0],ylim())
box off; set(gca,'TickDir','out')
ylabel('Normalized Rate')
line(xlim(),[1 1])
xlabel('Time following ripple onset [s]')
xlim([-.25 .25])
%%
% single pair...
% joint quantile hist, overall and during ripples...
% xcorr before and after ripples


% Fig 3A [14-14]
% Fig 3B [25-9]
% Fig 3C [1-6,29-3,25-28]
% Fig 3D [51-10,25-1,40-10]
% Fig 3E [2-27,32-9,34-9]
% rec_num=32; cnx=4; % example weakening w/ ripples
% rec_num=2; cnx=13; % example strengthening w/ ripples, ripple quiet
rec_num= 14; cnx=14; % fig 6d top - 
% rec_num= 25; cnx=9; % fig 6d top - cnxf=223
rec_num=49; cnx=21; % test
rec_num=27; cnx=14; % test
rec_num=49; cnx=27; % test
prefix = strsplit(glist(rec_num).name,'_');
prefix = prefix{1};
fprintf('%s...\n',prefix)

load([res_dir prefix '_detected_final'])

[r,c] = find(result.detected_cnx);
fit_res.ssi = result.ssi(sub2ind(size(result.ssi),r,c))';
fit_res.kidx = (fit_res.cvsmoo<.1) & fit_res.ssi<0.2 & (fit_res.rsq>.5) & (fit_res.dur<1) & (fit_res.lat>0.75) & (fit_res.lat<10) & (abs(fit_res.eff)<.5) & (abs(fit_res.eff)>.01) & fit_res.ex;
ijlist = [r(fit_res.kidx) c(fit_res.kidx)];
cnx_id = cnx_list(:,1)==rec_num;
cnxf = find(cnx_id);

% get results from initial xcorr
dwdir = 'E:\fluctuation_backup\dw\';
load([dwdir prefix '_joint_res_300'])

pre = ijlist(cnx,1);
post = ijlist(cnx,2);
[r,c] = find(result.detected_cnx);
fidx = find(r==pre & c==post);

T = max(cellfun(@max,Tlist));
dt=1e-1;

ripple_starts = arrayfun(@(x)x.start_time,ripple_all{cnxf(cnx)});
ripple_stops = arrayfun(@(x)x.end_time,ripple_all{cnxf(cnx)});

% get ripple locations
tin = linspace(0,T,ceil(T/dt));
pattern_ids = zeros(1,ceil(T/dt));
for r=1:size(ripple_starts,2)
    tid = tin>(ripple_starts(r)-.1) & tin<(ripple_stops(r)+.1);
    pattern_ids(tid) = (ripple_stops(r)-ripple_starts(r))/sum(tid);
end

% rates
[y_pre,t_pre,optw,W1,C1] = sskernel(Tlist{pre},tin);
pre_fr = y_pre/sum(y_pre)*length(Tlist{pre});
[y_post,t_post,optw2,W12,C12] = sskernel(Tlist{post},tin);
pos_fr = y_post/sum(y_post)*length(Tlist{post});

% joint rate statistics
pedges = linspace(1,length(pre_fr),10);
pre_rank = tiedrank(pre_fr);
pos_rank = tiedrank(pos_fr);
[N,Xedges,Yedges,binX,binY] = histcounts2(pre_rank,pos_rank,pedges,pedges);

[Nripp,Xedges,Yedges,binX,binY] = histcounts2(pre_rank(pattern_ids>0),pos_rank(pattern_ids>0),pedges,pedges);

pattern_ids = pattern_ids/sum(pattern_ids)*length(Tlist{ijlist(cnx,1)});
fitw = spk_xcorr_fwdrev(Tlist,ijlist(cnx,:),tin,pattern_ids,fitout(fidx).feat_params,params);


figure(8)
subplot(2,1,1)
plot(tin/60,fitw.fwd_pattern_mask,tin/60,fitw.rev_pattern_mask)
axis tight; box off; set(gca,'TickDir','out')
subplot(2,1,2)
bar(tin/60,pattern_ids,1,'EdgeColor','none')
box off; set(gca,'TickDir','out')

figure(9)
ax(1)=subplot(1,2,1)
bar(params.t,fitw.fwd.c/max(fitw.fwd.a),1,'EdgeColor','none')
hold on
plot(params.t,fitw.fwd.yhat)
hold off
box off; set(gca,'TickDir','out')
title('Fwd')
ax(2)=subplot(1,2,2)
bar(params.t,fitw.rev.c/max(fitw.rev.a),1,'EdgeColor','none')
hold on
plot(params.t,fitw.rev.yhat)
hold off
box off; set(gca,'TickDir','out')
title('Rev')
linkaxes(ax)
ylim([0 1.2*max(max(fitw.fwd.c/max(fitw.fwd.a)),max(fitw.rev.c/max(fitw.rev.a)))])

figure(10)
subplot(1,2,1)
imagesc(N)
axis image off; set(gca,'TickDir','out'); set(gca,'YDir','normal')
subplot(1,2,2)
imagesc(Nripp)
axis image off; set(gca,'TickDir','out'); set(gca,'YDir','normal')


%% save plots for examples
figure(8)
exportgraphics(gcf,[prefix '_' num2str(cnx) '_ripple_events.pdf'],'ContentType','vector')
figure(9)
exportgraphics(gcf,[prefix '_' num2str(cnx) '_ripple_xcorr.pdf'],'ContentType','vector')
figure(10)
exportgraphics(gcf,[prefix '_' num2str(cnx) '_ripple_jointrate.pdf'],'ContentType','vector')

%% loop to collect results on efficacy changes before and after ripples
cnxc=1;
freff=[]; freffse=[];
Nall=[]; Nrippall=[]; cellfeat=[];
riphist = [];
for rec_num=1:length(glist)
    prefix = strsplit(glist(rec_num).name,'_');
    prefix = prefix{1};
    fprintf('%s...\n',prefix)

    load([res_dir prefix '_detected_final'])

    [r,c] = find(result.detected_cnx);
    fit_res.ssi = result.ssi(sub2ind(size(result.ssi),r,c))';
    fit_res.kidx = (fit_res.cvsmoo<.1) & fit_res.ssi<0.2 & (fit_res.rsq>.5) & (fit_res.dur<1) & (fit_res.lat>0.75) & (fit_res.lat<10) & (abs(fit_res.eff)<.5) & (abs(fit_res.eff)>.01) & fit_res.ex;
    ijlist = [r(fit_res.kidx) c(fit_res.kidx)];
    cnx_id = cnx_list(:,1)==rec_num;
    cnxf = find(cnx_id);

    if ~isempty(ijlist)
        % get results from initial xcorr
        dwdir = 'E:\fluctuation_backup\dw\';
        load([dwdir prefix '_joint_res_300'])
    end

    for cnx=1:size(ijlist,1)
        fprintf('%02i...\n',cnx)
        pre = ijlist(cnx,1);
        post = ijlist(cnx,2);
        [r,c] = find(result.detected_cnx);
        fidx = find(r==pre & c==post);

        cellfeat(cnxc,1:2) = tab_sub{ijlist(cnx,:),'waveform_duration'}';
        cellfeat(cnxc,3:4) = tab_sub{ijlist(cnx,:),'firing_rate'}';
        

        T = max(cellfun(@max,Tlist));
        dt=1e-1;

        if ~isfield(ripple_all{cnxf(cnx)},'start_time')
            cnxc=cnxc+1;
            continue;
        end
        ripple_starts = arrayfun(@(x)x.start_time,ripple_all{cnxf(cnx)});
        ripple_stops = arrayfun(@(x)x.end_time,ripple_all{cnxf(cnx)});
        if isempty(ripple_starts)
            cnxc=cnxc+1;
            continue;
        end

        % get ripple locations
        tin = linspace(0,T,ceil(T/dt));
        pattern_ids = zeros(1,ceil(T/dt));
        for r=1:size(ripple_starts,2)
            tid = tin>(ripple_starts(r)-.1) & tin<(ripple_stops(r)+.1);
            pattern_ids(tid) = (ripple_stops(r)-ripple_starts(r))/sum(tid);
        end

        % ripple histogram
        riphist(cnxc,:) = histcounts(ripple_starts,0:(60):(180*60));
        riphist(cnxc,floor(T/60):end)=NaN;

        % rates
        [y_pre,t_pre,optw,W1,C1] = sskernel(Tlist{pre},tin);
        pre_fr = y_pre/sum(y_pre)*length(Tlist{pre});
        [y_post,t_post,optw2,W12,C12] = sskernel(Tlist{post},tin);
        pos_fr = y_post/sum(y_post)*length(Tlist{post});

        % joint rate statistics
        pedges = linspace(1,length(pre_fr),10);
        pre_rank = tiedrank(pre_fr);
        pos_rank = tiedrank(pos_fr);
        [N,Xedges,Yedges,binX,binY] = histcounts2(pre_rank,pos_rank,pedges,pedges);

        [Nripp,Xedges,Yedges,binX,binY] = histcounts2(pre_rank(pattern_ids>0),pos_rank(pattern_ids>0),pedges,pedges);

        pattern_ids = pattern_ids;
        fitw = spk_xcorr_fwdrev(Tlist,ijlist(cnx,:),tin,pattern_ids,fitout(fidx).feat_params,params);
        freff(cnxc,:) = [fitw.fwd.eff fitw.rev.eff];
        freffse(cnxc,:) = [fitw.fwd.effse fitw.rev.effse];
        Nall(:,:,cnxc)=N;
        Nrippall(:,:,cnxc)=Nripp;
        cnxc=cnxc+1;
    end
end

%%

figure(11)
subplot(1,2,1)
imagesc(mean(Nall(:,:,inclus),3))
axis image off; set(gca,'TickDir','out'); set(gca,'YDir','normal')
subplot(1,2,2)
imagesc(mean(Nrippall(:,:,inclus),3))
axis image off; set(gca,'TickDir','out'); set(gca,'YDir','normal')


figure(12)
clf
% errorbar(freff(inclus,1),freff(inclus,2),freffse(inclus,2),freffse(inclus,2),freffse(inclus,1),freffse(inclus,1),'o')
ny = [abs(log10(freff(inclus,2))-log10(freff(inclus,2)-freffse(inclus,2))) abs(log10(freff(inclus,2))-log10(freff(inclus,2)+freffse(inclus,2)))];
nx = [abs(log10(freff(inclus,1))-log10(freff(inclus,1)-freffse(inclus,1))) abs(log10(freff(inclus,1))-log10(freff(inclus,1)+freffse(inclus,1)))];
errorbar(log10(freff(inclus,2)),log10(freff(inclus,1)),nx(:,1),nx(:,2),ny(:,1),ny(:,2),'.','CapSize',0)
box off; set(gca,'TickDir','out')
line(xlim(),xlim())
xlim([-2.5 -.25])
ylim([-2.5 -.25])

hold on
pickids = [124; 223]; % examples from 6d
% pickids = [637]; % test outlier from 6d - heavily biased by ripples at beginning only
% pickids = [264]; % test another outlier
pickids = [643]; % test another outlier
plot(log10(freff(pickids,2)),log10(freff(pickids,1)),'ro')
hold off


set(gca,'XTick',log10([0.002:0.001:0.01 0.02:0.01:.1 .2:.1:.5]))
set(gca,'XTickLabel',arrayfun(@num2str,[0.002:0.001:0.01 0.02:0.01:.1 .2:.1:.5],'UniformOutput',false))
set(gca,'YTick',log10([0.002:0.001:0.01 0.02:0.01:.1 .2:.1:.5]))
set(gca,'YTickLabel',arrayfun(@num2str,[0.002:0.001:0.01 0.02:0.01:.1 .2:.1:.5],'UniformOutput',false))

xlabel('Before'); ylabel('After')

figure(13)
% histogram(freff(inclus,1)-freff(inclus,2),'EdgeColor','none')
histogram((freff(inclus,1)-freff(inclus,2))./freff(inclus,2),100,'EdgeColor','none')
% histogram(log2(freff(inclus,2))-log2(freff(inclus,1)),'EdgeColor','none')
box off; set(gca,'TickDir','out')
line([0 0],ylim())
% line([0 0]+mean(freff(inclus,1)-freff(inclus,2)),ylim())
xlim([-1 1])
xlabel('Efficacy Change [%]')


%%
figure(14)
subplot(2,2,1)
% plot(log(ripp_pair_stats(inclus,4)./ripp_pair_stats(inclus,5)),freff(inclus,1)-freff(inclus,2),'o')
plot(log(ripp_pair_stats(inclus,4)./ripp_pair_stats(inclus,5)),log(freff(inclus,1)./freff(inclus,2)),'o')
corrcoef(log(ripp_pair_stats(inclus,4)./ripp_pair_stats(inclus,5)),freff(inclus,1)-freff(inclus,2))
corrcoef(log(ripp_pair_stats(inclus,4)./ripp_pair_stats(inclus,5)),log(freff(inclus,1)./freff(inclus,2)))
subplot(2,2,2)
% plot(log(ripp_pair_stats(inclus,6)./ripp_pair_stats(inclus,7)),freff(inclus,1)-freff(inclus,2),'o')
plot(log(ripp_pair_stats(inclus,6)./ripp_pair_stats(inclus,7)),log(freff(inclus,1)./freff(inclus,2)),'o')
corrcoef(log(ripp_pair_stats(inclus,6)./ripp_pair_stats(inclus,7)),freff(inclus,1)-freff(inclus,2))
subplot(2,2,3)
plot(log(freff(inclus,1)),log(freff(inclus,1)./freff(inclus,2)),'o')
corrcoef(log(freff(inclus,1)),freff(inclus,1)-freff(inclus,2))

%%

figure(15)
clf

nq=10;
deltaw_id = tiedrank((freff(inclus,1)-freff(inclus,2))./freffse(inclus,1));
deltaw_id = ceil(deltaw_id/max(deltaw_id)*nq);
nain = sum(sum(Nall(:,:,inclus),1),2);
nrin = sum(sum(Nrippall(:,:,inclus),1),2);
Nai = Nall(:,:,inclus)./nain;
Nri = Nrippall(:,:,inclus)./nrin;

for i=1:max(deltaw_id)
    subplot(2,nq,i)
    imagesc(mean(Nai(:,:,deltaw_id==i),3))
    axis image
    set(gca,'YDir','normal')
    subplot(2,nq,i+nq)
    imagesc(mean(Nri(:,:,deltaw_id==i),3))
    axis image
    set(gca,'YDir','normal')
end

%% ripple pattern summary

a=nanmean(riphist,2).^2./nanmean(riphist.^2,2);
S=(1-a)./(1-1./sum(isfinite(riphist),2));
