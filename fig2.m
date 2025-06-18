

dparams=[];
dparams.trange = [-0.025 0.025];
dparams.tn = 102;
dparams.t = linspace(dparams.trange(1),dparams.trange(2),dparams.tn);
dparams.t = dparams.t(1:end-1)+mean(diff(dparams.t))/2;
dparams.jit_min_t = 0.0011;
dparams.ssi_min_t = 0.00049;
dparams.ssi_max = 0.4;
dparams.alph = 0.1;
dparams.jit_conv_window = 2;
dparams.mask=abs(dparams.t)>0.0000001;


%% evaluate single recording

res_dir = 'C:\Users\istevenson\Documents\GitHub\fluctuation_project\detection\';

glist = dir([res_dir '*_detected_final.mat']);
rec_num = 1
prefix = strsplit(glist(rec_num).name,'_');
prefix = prefix{1};
load([res_dir prefix '_detected_final'])

[r,c] = find(result.detected_cnx);

fit_res.ssi = result.ssi(sub2ind(size(result.ssi),r,c))';
fit_res.kidx = (fit_res.cvsmoo<.1) & fit_res.ssi<0.2 & (fit_res.rsq>.5) & (fit_res.dur<1) & (fit_res.lat>0.75) & (fit_res.lat<10) & (abs(fit_res.eff)<.5) & (abs(fit_res.eff)>.01) & fit_res.ex;
ijlist = [r(fit_res.kidx) c(fit_res.kidx)];

sidx = find(result.detected_cnx(:));
sidx = sidx(fit_res.kidx);
figure(1)
offset=1;
plotPairs(result,dparams,sidx,offset)

[sum(fit_res.kidx) length(fit_res.kidx)]

%%
[r c fit_res.cvsmoo' fit_res.rsq' fit_res.lat' fit_res.dur' fit_res.eff' fit_res.kidx' tab_sub.waveform_duration(r) fit_res.ex']


%% Load all detection results


res_dir = 'C:\Users\istevenson\Documents\GitHub\fluctuation_project\detection\';

glist = dir([res_dir '*_detected_final.mat']);

area_count = zeros(1,5);
det_count = zeros(5,5);
kid_count = zeros(5,5);
pair_count = zeros(5,5);

wave_all = [];
kid_params_all = [];
all_anat=[];

for rec_num = 1:length(glist)
    % for rec_num = 1:21
    prefix = strsplit(glist(rec_num).name,'_');
    prefix = prefix{1};

    load([res_dir prefix '_detected_final'])
    %     load([prefix '_detected_2knot'])

    all_anat = [all_anat; data.anatomy];

    rparams=[];
    rparams.trange = [-0.025 0.025];
    rparams.tn = 102;
    rparams.t = linspace(rparams.trange(1),rparams.trange(2),rparams.tn);
    rparams.t = rparams.t(1:end-1)+mean(diff(rparams.t))/2;
    rparams.jit_min_t = 0.0011;
    rparams.ssi_min_t = 0.00049;
    rparams.ssi_max = 0.4;
    rparams.alph = 0.1;
    rparams.jit_conv_window = 2;
    rparams.mask=abs(rparams.t)>0.0000001;

    amat = data.area_mat(data.ids,:);

    area_count = area_count+sum(data.area_mat);
    [r,c] = find(result.detected_cnx);
    for i=1:length(r)
        rid = find(amat(r(i),:));
        cid = find(amat(c(i),:));
        det_count(rid,cid)=det_count(rid,cid)+1;
    end
    [r,c] = find(result.detected_cnx);
    pair_count = pair_count + max(0,sum(data.area_mat)'*sum(data.area_mat)-diag(sum(data.area_mat)));


    fit_res.ssi = result.ssi(sub2ind(size(result.ssi),r,c))';
    fit_res.kidx = (fit_res.cvsmoo<.1) & fit_res.ssi<0.2 & (fit_res.rsq>.5) & (fit_res.dur<1) & (fit_res.lat>0.75) & (fit_res.lat<10) & (abs(fit_res.eff)<.5) & (abs(fit_res.eff)>.01) & fit_res.ex;
    ijlist = [r(fit_res.kidx) c(fit_res.kidx)];

    for i=1:size(ijlist,1)
        rid = find(amat(ijlist(i,1),:));
        cid = find(amat(ijlist(i,2),:));
        kid_count(rid,cid)=kid_count(rid,cid)+1;
    end
    size(ijlist,1)

    pre_vec = zeros(size(tab_sub,1),1);
    pos_vec = zeros(size(tab_sub,1),1);
    if ~isempty(ijlist)
        pre_vec(ijlist(:,1)) = 1;
        pos_vec(ijlist(:,2)) = 1;
    end
    wave_all = [wave_all; [tab_sub.waveform_halfwidth tab_sub.waveform_duration tab_sub.firing_rate pre_vec pos_vec]];

    kid_params_all = [kid_params_all; [fit_res.rsq(fit_res.kidx)' fit_res.lat(fit_res.kidx)' fit_res.dur(fit_res.kidx)' fit_res.eff(fit_res.kidx)']];
end

area_count
sum(pair_count(:)) % pairs
sum(det_count(:)) % initial pass
sum(kid_count(:)) % model pass

%% nice waveform plot

figure(12)
[n1,x1e,y1e] = histcounts2(log10(wave_all(:,1)),log10(wave_all(:,2)),256);
[n2,x2e,y2e] = histcounts2(log10(wave_all(:,3)),log10(wave_all(:,2)),256);
x1 = x1e(1:end-1)+diff(x1e)/2;
y1 = y1e(1:end-1)+diff(y1e)/2;
x2 = x2e(1:end-1)+diff(x2e)/2;
y2 = y2e(1:end-1)+diff(y2e)/2;
gwin = gausswin(64);
gwin = gwin*gwin';
s1 = conv2(n1,gwin,'same');
s2 = conv2(n2,gwin,'same');

subplot(3,5,[6 7 11 12])
% imagesc(x1,y1,s1')
imagesc(x1,y1,log(s1'+50))
set(gca,'YDir','normal')
set(gca,'TickDir','out')
y1r = unique(round((10.^y1)*10)/10);
set(gca,'YTick',log10(y1r))
set(gca,'YTickLabels',y1r)
x1r = unique(round((10.^x1)*10)/10);
set(gca,'XTick',log10(x1r))
set(gca,'XTickLabels',x1r)
box off
xlabel('Waveform Half-width [ms]')
ylabel('Waveform Duration [ms]')
xlim([min(x1) log10(0.6)])
ylim([min(y1) log10(1.4)])


cmap = lines(2);
% avoid overplottting
jit_wave = wave_all(:,1:2);
jit_wave(:,1) = jit_wave(:,1)+(rand(size(jit_wave(:,1)))-.5)*0.0137;
jit_wave(:,2) = jit_wave(:,2)+(rand(size(jit_wave(:,2)))-.5)*0.0137;

hold on
scatter(log10(jit_wave(wave_all(:,4)>0,1)),log10(jit_wave(wave_all(:,4)>0,2)),10,'MarkerFaceColor',cmap(1,:),'MarkerEdgeColor','none')
scatter(log10(jit_wave(wave_all(:,5)>0,1)),log10(jit_wave(wave_all(:,5)>0,2)),10,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','none')
hold off

subplot(3,5,[8 9 13 14])
% imagesc(x2,y2,s2')
imagesc(x2,y2,log(s2'+50))
set(gca,'YDir','normal')
set(gca,'TickDir','out')
set(gca,'YDir','normal')
set(gca,'TickDir','out')
y2r = unique(round((10.^y2)*10)/10);
set(gca,'YTick',log10(y2r))
set(gca,'YTickLabels',y2r)
% x2r = unique(round((10.^x2)));
x2r = 2.^[-3:6];
set(gca,'XTick',log10(x2r))
set(gca,'XTickLabels',x2r)
box off
colormap(1-gray(256))
xlabel('Firing Rate [Hz]')
ylabel('Waveform Duration [ms]')
% xlim([min(x1) log10(0.6)])
ylim([min(y1) log10(1.4)])
xl=xlim();

hold on
scatter(log10(wave_all(wave_all(:,4)>0,3)),log10(jit_wave(wave_all(:,4)>0,2)),10,'MarkerFaceColor',cmap(1,:),'MarkerEdgeColor','none')
scatter(log10(wave_all(wave_all(:,5)>0,3)),log10(jit_wave(wave_all(:,5)>0,2)),10,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','none')
hold off

subplot(3,5,1:2)
edges = linspace(log10(0.06),log10(0.65),256);
[fi,x]=ksdensity(log10(jit_wave(wave_all(:,4)>0,1)),edges,'bandwidth',range(edges)/15);
[fe,x]=ksdensity(log10(jit_wave(wave_all(:,5)>0,1)),edges,'bandwidth',range(edges)/15);
% plot(x,fi,x,fe)
shadedErrorBar(x,fi*0,[fi; fi*0])
hold on
shadedErrorBar(x,fe*0,[fe; fe*0])
hold off
box off; set(gca,'TickDir','out')
xlim([min(x1) log10(0.6)])

subplot(3,5,3:4)
edges = linspace(log10(0.1),log10(100),256);
[fi,x]=ksdensity(log10(wave_all(wave_all(:,4)>0,3)),edges,'bandwidth',range(edges)/15);
[fe,x]=ksdensity(log10(wave_all(wave_all(:,5)>0,3)),edges,'bandwidth',range(edges)/15);
% plot(x,fi,x,fe)
shadedErrorBar(x,fi*0,[fi; fi*0])
hold on
shadedErrorBar(x,fe*0,[fe; fe*0])
hold off
box off; set(gca,'TickDir','out')
xlim(xl)

subplot(3,5,[10 15])
edges = linspace(log10(0.08),log10(1.5),256);
[fi,x]=ksdensity(log10(wave_all(wave_all(:,4)>0,2)),edges,'bandwidth',range(edges)/15);
[fe,x]=ksdensity(log10(wave_all(wave_all(:,5)>0,2)),edges,'bandwidth',range(edges)/15);
% plot(x,fi,x,fe)
shadedErrorBar(x,fi*0,[fi; fi*0])
hold on
shadedErrorBar(x,fe*0,[fe; fe*0])
hold off
box off; set(gca,'TickDir','out')
camroll(90)
set(gca,'XDir','reverse')
xlim([min(y1) log10(1.4)])

%% putative synapse statistics

figure(13)
subplot(3,1,1)
histogram(kid_params_all(:,2),25,'EdgeColor','none')
xlabel('Latency [ms]')
box off; set(gca,'TickDir','out')
subplot(3,1,2)
histogram(kid_params_all(:,3),25,'EdgeColor','none')
xlabel('Time Constant [ms]')
box off; set(gca,'TickDir','out')
subplot(3,1,3)
tmp = kid_params_all(:,4)*100;
tmp(tmp>23)=23;
histogram(tmp,25,'EdgeColor','none')
xlabel('Efficacy [%]')
box off; set(gca,'TickDir','out')
xlim([0 23])

%% area counts

figure(14)
subplot(5,1,1)
imagesc(log(area_count))
text(1:5,ones(1,5)-.4,data.area_list,'HorizontalAlignment','center')
text(1:5,ones(1,5),arrayfun(@num2str,area_count,'UniformOutput',false),'HorizontalAlignment','center')
box off; axis off

subplot(5,1,2:5)
imagesc(log(pair_count))
for i=1:size(kid_count,1)
    text(1:5,ones(1,5)+i-1,arrayfun(@num2str,kid_count(i,:),'UniformOutput',false),'HorizontalAlignment','center')
end
axis equal
box off; axis off
(kid_count./pair_count)*100
colormap(1-gray(256))
colormap(cbrew(256,'blues'))

%%

ij = [126 134]
fidx = find(r==ij(1) & c==ij(2));
fsyn = fitout(fidx);

params.nknots=2;
params.doPlot = false;
figure(2)
[opt_res,smooth_estimate] = optimizeWindowSize(ij,Tlist,params,fsyn,20,1);

%% optimize window size for all classified putative synapses

for cnx = 1:size(ijlist,1)
    fidx = find(r==ij(cnx) & c==ij(cnx));
    fsyn = fitout(fidx);

    params.nknots=2;
    params.doPlot = false;
    figure(2)
    [opt_res(cnx),smooth_estimate(cnx)] = optimizeWindowSize(ijlist(cnx,:),Tlist,params,fsyn,20,1);
end

%%

figure(31)
cnx=6;
subplot(3,10,2:10)
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
% clf
subplot(3,10,11)
bar(params.t,cc,1,'EdgeColor','none')
hold on
plot(params.t,fitout(fidx).yhat*length(Tlist{ijlist(cnx,1)}))
hold off
box off; set(gca,'TickDir','out')
line([0 0],ylim())
subplot(3,10,20+1)
bar(params.t,ac,1,'EdgeColor','none')
box off; set(gca,'TickDir','out')
line([0 0],ylim())


ww = smooth_estimate(cnx).window;
ww = 10*60;
for i=2:10
    tid = (Tlist{ijlist(cnx,1)}>(tcenter(i)-ww/2) & Tlist{ijlist(cnx,1)}<(tcenter(i)+ww/2));
    cc(:,i) = spk_xcorr(Tlist{ijlist(cnx,1)}(tid),Tlist{ijlist(cnx,2)},params.interval(1),params.interval(2),params.num_bins);
    ac(:,i) = spk_xcorr(Tlist{ijlist(cnx,1)}(tid),Tlist{ijlist(cnx,1)},params.interval(1),params.interval(2),params.num_bins);
    ac(params.t==0,i)=NaN;
    cc(:,i) = cc(:,i)/sum(tid);
    ac(:,i) = ac(:,i)/sum(tid);
end

for i=2:10
    subplot(3,10,10+i)
    bar(params.t,cc(:,i),1,'EdgeColor','none')
    box off; set(gca,'TickDir','out')
    line([0 0],ylim())
    ylim([0 max(max(cc(:,2:end)))])
    subplot(3,10,20+i)
    bar(params.t,ac(:,i),1,'EdgeColor','none')
    box off; set(gca,'TickDir','out')
    line([0 0],ylim())
    ylim([0 max(max(ac(:,2:end)))])
end


%%

res_dir = 'C:\Users\istevenson\Documents\GitHub\fluctuation_project\detection\';

glist = dir([res_dir '*_detected_final.mat']);

% opt_win=[];
% for rec_num = 1:length(glist)
for rec_num=19
    prefix = strsplit(glist(rec_num).name,'_');
    prefix = prefix{1};

    load([res_dir prefix '_detected_final'])

    [r,c] = find(result.detected_cnx);
    fit_res.ssi = result.ssi(sub2ind(size(result.ssi),r,c))';
    fit_res.kidx = (fit_res.cvsmoo<.1) & fit_res.ssi<0.2 & (fit_res.rsq>.5) & (fit_res.dur<1) & (fit_res.lat>0.75) & (fit_res.lat<10) & (abs(fit_res.eff)<.5) & (abs(fit_res.eff)>.01) & fit_res.ex;
    ijlist = [r(fit_res.kidx) c(fit_res.kidx)];


    load([res_dir prefix '_opt_window'])
    %     if ~isempty(dir([res_dir prefix '_opt_window.mat']))
    %         load([res_dir prefix '_opt_window'])
    %         opt_win = [opt_win arrayfun(@(x) x.window,smooth_estimate)];
    %     end
end


%%

res_dir = 'C:\Users\istevenson\Documents\GitHub\fluctuation_project\detection\';

glist = dir([res_dir '*_detected_final.mat']);

% opt_win=[];
% for rec_num = 1:length(glist)
for rec_num=[1:length(glist)]
    prefix = strsplit(glist(rec_num).name,'_');
    prefix = prefix{1};

    load([res_dir prefix '_detected_final'])

    [r,c] = find(result.detected_cnx);
    fit_res.ssi = result.ssi(sub2ind(size(result.ssi),r,c))';
    fit_res.kidx = (fit_res.cvsmoo<.1) & fit_res.ssi<0.2 & (fit_res.rsq>.5) & (fit_res.dur<1) & (fit_res.lat>0.75) & (fit_res.lat<10) & (abs(fit_res.eff)<.5) & (abs(fit_res.eff)>.01) & fit_res.ex;
    ijlist = [r(fit_res.kidx) c(fit_res.kidx)];


    %     load([res_dir prefix '_opt_window'])
    if ~isempty(dir([res_dir prefix '_opt_window.mat']))
        load([res_dir prefix '_opt_window'])
        opt_win = [opt_win arrayfun(@(x) x.window,smooth_estimate)];


        for cnx=1:size(ijlist,1)
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

            exportgraphics(gcf,sprintf('%03i_%03i_fluctuation.pdf',rec_num,cnx),'ContentType','vector')
        end
    end
end