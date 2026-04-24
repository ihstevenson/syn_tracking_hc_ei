
% takes [opt_res,smooth_estimate] result from optimal window (see fig2d.m)

% optimal bandwidth rate estimates at dt resolution
dt=0.1; % time resolution for rate effects [s]
tin = linspace(0,T(pair_id),ceil(T(pair_id)/dt));
[y_pre,t_pre,optw,W1,C1] = sskernel(Tlist{1},tin);
pre_fr = y_pre/sum(y_pre)*length(Tlist{1});
[y_post,t_post,optw2,W12,C12] = sskernel(Tlist{2},tin);
pos_fr = y_post/sum(y_post)*length(Tlist{2});

% get (anti-aliased) rate, that matches efficacy windows
% smooth_estimate using smooth_estimate.window windows w/ 10 stride
tid = interp1(tin,1:length(tin),smooth_estimate.tcenter,'linear');
[fpre,ty] = resample(pre_fr,tin,0.1,1,100,'linear');
fpre10 = interp1(ty,fpre,smooth_estimate.tcenter,'linear')';
[fpos,ty] = resample(pos_fr,tin,0.1,1,100,'linear');
fpos10 = interp1(ty,fpos,smooth_estimate.tcenter,'linear')';
fpre = smoothdata(fpre10,"movmean",smooth_estimate.window/10);
fpos = smoothdata(fpos10,"movmean",smooth_estimate.window/10);

% fit STP model
y = smooth_estimate.eff;
yse = smooth_estimate.effse.^2;
X = zscore(real([log(fpre) log(fpos)]));
validid=all(X>-10,2) & all(isfinite(X),2) & isfinite(y);
[b,dev,stat] = glmfit(X(validid,:),y(validid),'normal','weights',1./yse(validid));
effhat = glmval(b,X,'identity');
aic(1) = 2*length(b)+dev;

% fit STP+BCM models
% get covariance and 3rd order terms from high-resolution rate estimates
fpp_hi = pre_fr.*pos_fr;
fp2p_hi = pre_fr.*pos_fr.^2;

% filter and downsample to resolution of efficacy observations
tau = 30;
twin = linspace(0,tau*10,tau*10/0.1+1);
wfun = twin/tau.*exp(1-twin/tau);
wfun = wfun/sum(wfun);
fpp_hi_filt = conv(fpp_hi,wfun);
fpp_hi_filt = fpp_hi_filt(1:length(fpp_hi));
fpp_filt = interp1(tin,fpp_hi_filt,smooth_estimate.tcenter,'linear')';
fp2p_hi_filt = conv(fp2p_hi,wfun);
fp2p_hi_filt = fp2p_hi_filt(1:length(fp2p_hi));
fp2p_filt = interp1(tin,fp2p_hi_filt,smooth_estimate.tcenter,'linear')';
fpp_filt = smoothdata(fpp_filt,"movmean",smooth_estimate.window/10);
fp2p_filt = smoothdata(fp2p_filt,"movmean",smooth_estimate.window/10);

% STP+covariance
X = zscore(real([log(fpre) log(fpos) fpp_filt]));
[b2,dev,stat] = glmfit(X(validid,:),y(validid),'normal','weights',1./yse(validid));
effhat(:,2) = glmval(b2,X,'identity');
aic(2) = 2*length(b2)+dev;

% STP+BCM
X = zscore(real([log(fpre) log(fpos) fpp_filt fp2p_filt]));
[b3,dev,stat] = glmfit(X(validid,:),y(validid),'normal','weights',1./yse(validid));
effhat(:,3) = glmval(b3,X,'identity');
aic(3) = 2*length(b3)+dev;

% goodness-of-fit
tmp = corrcoef([y(validid) effhat(validid,:)]);
r2 = tmp(1,2:end).^2;

% show traces
figure(51)
clf
subplot(2,1,1)
shadedErrorBar(smooth_estimate.tcenter/60,smooth_estimate.eff,smooth_estimate.effse)
line(xlim(),[0 0])
axis tight; box off;
set(gca,'TickDir','out')
xlabel('Time [min]')
ylabel('Efficacy')
hold on
plot(smooth_estimate.tcenter/60,effhat)
hold off
subplot(2,1,2)
plot(smooth_estimate.tcenter/60,zscore(fpre10))
hold on
plot(smooth_estimate.tcenter/60,zscore(fpos10))
plot(smooth_estimate.tcenter/60,zscore(fpp_filt))
hold off
ylabel('Features')

%% visualize learning rule as delta-efficacy matrix

cent = linspace(0,100,11);
cent = cent(1:end-1)+diff(cent)/2;
pre_pct = prctile(pre_fr,cent);
pos_pct = prctile(pos_fr,cent);
[rpos,rpre] = meshgrid(pos_pct,pre_pct);

zzm = mean([fpp_filt fp2p_filt]);
zzs = std([fpp_filt fp2p_filt]);

figure(52)
imagesc(b3(4)*(rpre.*rpos-zzm(1))./zzs(1)+b3(5)*(rpre.*rpos.^2-zzm(2))/zzs(2))
set(gca,'YDir','normal')
cl=get(gca,'CLim')
set(gca,'CLim',[-1 1]*max(abs(cl)))
colormap turbo
colorbar


