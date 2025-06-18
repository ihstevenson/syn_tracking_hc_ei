
function [opt_res,smooth_estimate] = optimizeWindowSize(ij,Tlist,params,fitout,vecn,doPlot)

% fit overall correlogram
% params.nknots = 4;
% params.doPlot = true;
% params.mask = ones(size(params.t));
% params.mask(abs(params.t)<0.0001)=0;
% params.mask = params.mask>0;
% params.asym=false;

pre = ij(1);
post = ij(2);

if nargin<4
    figure(22)
    fitout = composite_model_fit([pre post],Tlist,params,'i>j');
end
if nargin<5
    vecn=10;
end
if nargin<6
    doPlot=true;
end


% efficacy over time with different window sizes
% params.nknots = 2;
% params.doPlot = false;
win_vec = logspace(log10(60),log10(600),vecn); % 10s to 10min timescale

test_window = 10; % resolution to use for testing
cvk = 5;
T = max(cellfun(@max,Tlist));

eff_multi = cell(0);
total_f=[];
fitout_local_all = cell(0);

edges = 0:test_window:(ceil(T/test_window)*test_window);
[~,~,bin] = histcounts(Tlist{pre},edges);
cv = cvpartition(length(Tlist{pre}),'KFold',cvk);

tic
for k=1:cvk
    for w=1:length(win_vec)
        fprintf('window %0.2f...',win_vec(w))
        window = win_vec(w);
        stride = test_window;

        tcenter=[]; ydat=[];
        tstart=0;
        i=1;
        total_f(w,1:2,k) = 0;


        while (tstart+window)<T
            pre_spk_id = (Tlist{pre}>tstart & Tlist{pre}<(tstart+window)) & cv.training(k);
            Zlist = Tlist;
            Zlist{post} = Zlist{post}(Zlist{post}>tstart & Zlist{post}<(tstart+window));
            tcenter(i) = tstart+window/2;
            if sum(pre_spk_id)>1 & length(Zlist{post})>1
                fitout_local_all{w,i} = composite_model_fit([pre post],Zlist,params,'i>j',pre_spk_id,fitout.feat_params);
                ydat(i,1) = fitout_local_all{w,i}.eff;
                ydat(i,2) = fitout_local_all{w,i}.effse;
                total_f(w,1,k) = total_f(w,1,k)+fitout_local_all{w,i}.f;
            else
                fitout_local_all{w,i} = [];
                ydat(i,1) = NaN;
                ydat(i,2) = NaN;
            end
            tstart=tstart+stride;
            i=i+1;
        end
        eff_multi{w} = [tcenter' ydat];

        % test set - fixed short window...
        for i=1:max(bin)
            pre_spk_id = bin==i & cv.test(k);
            if sum(pre_spk_id)>1
                [c,~] = spk_xcorr(Tlist{pre}(pre_spk_id),Tlist{post},params.interval(1),params.interval(2),params.num_bins);
                [~,model_i] = min(abs(tcenter-(i-1)*test_window+test_window/2));
                if ~isempty(fitout_local_all{w,model_i})
                    p=fitout_local_all{w,model_i}.yhat;
                    p(p==0)=eps;
                    N=sum(pre_spk_id);
%                     total_f(w,2,k) = total_f(w,2,k)-2*sum(c.*log(p)+(N-c).*log(1-p)+gammaln(N+1)-gammaln(c+1)-gammaln(N-c+1));
                    total_f(w,2,k) = total_f(w,2,k)-2*nansum(c.*log(p)+(N-c).*log(1-p));
                else
                    total_f(w,2,k) = NaN;
                end
            end
        end

        toc
    end
end

opt_res.win_vec = win_vec;
opt_res.dev = total_f;
opt_res.eff_multi = eff_multi;


[~,min_dev] = min(nanmean(total_f(:,2,:),3));
if doPlot
    % plot optimization result
    figure(24)
    semilogx(win_vec,squeeze(total_f(:,2,:)),'ko')
    hold on
    semilogx(win_vec,nanmean(total_f(:,2,:),3),'rx')
    box off; set(gca,'TickDir','out')
    xlabel('Window Size [s]')
    ylabel('Cross-Validated Deviance')
    
%     yyy = nanmean(total_f(:,2,:),3);
%     nkeep = isfinite(yyy);
%     p = polyfit(log(win_vec(nkeep)),yyy(nkeep),5);
%     x0 = logspace(log10(min(win_vec(nkeep))),log10(max(win_vec(nkeep))),512);
%     yhat = polyval(p,log(x0));
%     hold on
%     plot(x0,yhat)
%     hold off
%     [~,ix0]=min(yhat);
%     optw = x0(ix0)

    % plot smoothing at different levels
    figure(25)
    ivec = [1 min_dev vecn];
    for i=1:3
        id = ivec(i);
        ax(i)=subplot(3,1,i)
        plot(eff_multi{id}(:,1)/60,eff_multi{id}(:,2))
        hold on
        plot(eff_multi{id}(:,1)/60,eff_multi{id}(:,2)+eff_multi{id}(:,3),'k')
        plot(eff_multi{id}(:,1)/60,eff_multi{id}(:,2)-eff_multi{id}(:,3),'k')
        hold off
        box off; set(gca,'TickDir','out')
        ylabel(sprintf('Eff [win =%.1fs]',win_vec(id)))
    end
    xlabel('Time [min]')
    linkaxes(ax)
    ylim([-0.1 0.5])
end


% run the optimized window size, short stride, full data
fprintf('running optimized window size...\n')

window = win_vec(min_dev);
stride = test_window;

fitout_local_opt = cell(0);
tcenter=[]; ydat=[];
tstart=0; call=[];
i=1;
while (tstart+window)<T
    pre_spk_id = (Tlist{pre}>tstart & Tlist{pre}<(tstart+window));
    Zlist = Tlist;
    Zlist{post} = Zlist{post}(Zlist{post}>tstart & Zlist{post}<(tstart+window));
    tcenter(i) = tstart+window/2;
    if sum(pre_spk_id)>0 & length(Zlist{post})>0
        fitout_local_opt{i} = composite_model_fit([pre post],Zlist,params,'i>j',pre_spk_id,fitout.feat_params);
        ydat(i,1) = fitout_local_opt{i}.eff;
        ydat(i,2) = fitout_local_opt{i}.effse;

        [c,~] = spk_xcorr(Tlist{pre}(pre_spk_id),Tlist{post},params.interval(1),params.interval(2),params.num_bins);
        call(i,:) = c./sum(pre_spk_id);
    else
        ydat(i,1) = NaN;
        ydat(i,2) = NaN;
    end

    tstart=tstart+stride;
    i=i+1;
end
smooth_estimate.tcenter= tcenter;
smooth_estimate.window = window;
smooth_estimate.eff = ydat(:,1);
smooth_estimate.effse = ydat(:,2);



if doPlot
    % show optimized window result
    tvec = [11:20:160]*60;
    figure(26)
    subplot(3,1,1:2)
    plot(tcenter/60,ydat(:,1))
    hold on
    plot(tcenter/60,ydat(:,1)+ydat(:,2),'k')
    plot(tcenter/60,ydat(:,1)-ydat(:,2),'k')
    yl=ylim();
    plot(tvec/60,yl(2)*.9,'r^')
    hold off
    box off; set(gca,'TickDir','out')

    xlim([0 max(tcenter/60)])
    subplot(3,1,3)
    imagesc(tcenter/60,params.t*1000,call')
    box off; set(gca,'TickDir','out')

    % pick a few time points to show
    ax=[];
    figure(27)
    for i=1:length(tvec)
        [~,tid] = min(abs(tcenter-tvec(i)));
        ax(i)=subplot(1,length(tvec),i);
        bar(params.t,call(tid,:),1,'EdgeColor','flat')
        box off; set(gca,'TickDir','out')
        title(sprintf('%0.2f min',tcenter(tid)/60))
    end
    linkaxes(ax)
    for i=1:length(tvec)
        subplot(1,length(tvec),i)
        line([0 0],ylim())
    end
end