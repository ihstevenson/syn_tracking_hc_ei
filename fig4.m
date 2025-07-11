
%%

% estimate firing rates
pair_id = 8; % from data subset
Tlist(1) = Tlist_pre(pair_id);
Tlist(2) = Tlist_post(pair_id);

dt=0.1; % time resolution for rate effects [s]
tin = linspace(0,T(pair_id),ceil(T(pair_id)/dt));
[y_pre,t_pre,optw,W1,C1] = sskernel(Tlist{1},tin);
pre_fr = y_pre/sum(y_pre)*length(Tlist{1});
[y_post,t_post,optw2,W12,C12] = sskernel(Tlist{2},tin);
pos_fr = y_post/sum(y_post)*length(Tlist{2});

% joint rate statistics
pedges = linspace(1,length(pre_fr),10);
pre_rank = tiedrank(pre_fr);
pos_rank = tiedrank(pos_fr);
[N,Xedges,Yedges,binX,binY] = histcounts2(pre_rank,pos_rank,pedges,pedges);

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

fitout = composite_model_fit([1 2],Tlist,params,'i>j');

% Fit forward and reverse weighted models relative to pre-post rate pairings
nsplit=10; % matrix resolution will be [nsplit x nsplit]
dfac = length(pre_rank)/nsplit;
pedges = (1:(nsplit)+1)-.5;
[Nlo,~,~,binXlo,binYlo] = histcounts2(ceil(pre_rank/dfac),ceil(pos_rank/dfac),pedges,pedges);
for i=1:nsplit
    for j=1:nsplit
        pattern_ids = binXlo==i & binYlo==j;
        fitw(i,j) = spk_xcorr_fwdrev(Tlist,[1 2],tin,pattern_ids,fitout.feat_params,params);
        fprintf('.')
    end
    fprintf('\n')
end


%%
figure(41)
clf
subplot(5,5,1)
imagesc(N)
set(gca,'TickDir','out'); set(gca,'YDir','normal'); box off; axis equal
axis image
ylabel('Pre')
xlabel('Post')

subplot(5,5,3:5)
plot(tin/60,pre_fr/dt,tin/60,pos_fr/dt)
box off; set(gca,'TickDir','out')
axis tight

% show weighted ccgs for corner pre-post rate combinations
c=1; ax=[];
for i=[1 nsplit]
    for j=[1 nsplit]
        fprintf('.')
        figure(41)
        ax(2*c-1)=subplot(5,5,c*5+1)
        box off; set(gca,'TickDir','out')
        bar(params.t,fitw(i,j).rev.c,1,'EdgeColor','flat')
        hold on
        plot(params.t,fitw(i,j).rev.yhat*fitw(i,j).rev.N)
        hold off

        ax(2*c)=subplot(5,5,c*5+2)
        bar(params.t,fitw(i,j).fwd.c,1,'EdgeColor','flat')
        hold on
        plot(params.t,fitw(i,j).fwd.yhat*fitw(i,j).fwd.N)
        hold off
        box off; set(gca,'TickDir','out')

        subplot(5,5,c*5+(3:5))
        plot(tin/60,fitw(i,j).fwd_pattern_mask/max(fitw(i,j).fwd_pattern_mask)/2,tin/60,fitw(i,j).rev_pattern_mask/max(fitw(i,j).rev_pattern_mask)/2)
        hold on
        hedges=linspace(0,T(pair_id)/60,256);
        histogram(find(fitw(i,j).pattern_ids)*dt/60,hedges,'Normalization','pdf','EdgeColor','none')
        hold off
        box off; set(gca,'TickDir','out')
        axis tight
        c=c+1;
    end
end
linkaxes(ax,'y')

%% Efficacy matrix and BCM fit

w=[];
for i=1:size(fitw,1)
    for j=1:size(fitw,2)
        w(i,j,1) = fitw(i,j).fwd.eff;
        w(i,j,2) = fitw(i,j).rev.eff;
    end
end
dw = (w(:,:,1)-w(:,:,2))*100;

% Fit BCM to observed efficacy change matrix
cent = linspace(0,100,11);
cent = cent(1:end-1)+diff(cent)/2;
pre_pct = prctile(pre_fr,cent);
pos_pct = prctile(pos_fr,cent);

[rpos,rpre] = meshgrid(pos_pct,pre_pct);
X = [rpos(:).^2.*rpre(:) -rpos(:).*rpre(:) rpre(:) rpos(:)];
b = glmfit(X,dw(:),'normal');
dwhat = glmval(b,X,'identity');
dwhat = reshape(dwhat,size(fitw));


% Visualize efficacy matrices
figure(42)
set(gca,'CLim',[0 max(get(gca,'CLim'))])

subplot(2,2,1)
imagesc(w(:,:,2))
set(gca,'YDir','normal')
cl(1,:)=get(gca,'CLim');
colorbar
title('Efficacy Before')

subplot(2,2,2)
imagesc(w(:,:,1))
set(gca,'YDir','normal')
colorbar
title('Efficacy After')
cl(2,:)=get(gca,'CLim');
set(gca,'CLim',[min(cl(:,1)) max(cl(:,2))])
subplot(2,2,1)
set(gca,'CLim',[min(cl(:,1)) max(cl(:,2))])

subplot(2,2,3)
imagesc(dw)
set(gca,'YDir','normal')
cl=get(gca,'CLim');
set(gca,'Clim',[-1 1]*max(abs(cl)));
colorbar
title('Observed \DeltaEfficacy')

subplot(2,2,4)
imagesc(dwhat)
set(gca,'YDir','normal')
set(gca,'Clim',[-1 1]*max(abs(cl)));
colorbar
axis equal tight
title('Model \DeltaEfficacy')
xlabel('')