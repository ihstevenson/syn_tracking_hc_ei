
addpath('C:\Users\istevenson\Dropbox\__matlab_toolboxes')

res_dir = 'C:\Users\istevenson\Documents\GitHub\fluctuation_project\detection\';
bas_dir = 'D:\allen_neuropixel\';
glist = dir([res_dir '*_detected_final.mat']);
load('C:\Users\istevenson\Documents\GitHub\fluctuation_project\all_smooth_estimate')
load('C:\Users\istevenson\Documents\GitHub\fluctuation_project\all_smooth_estimate_extra')

rec_num=19; cnx=14;
prefix = strsplit(glist(rec_num).name,'_');
prefix = prefix{1};
ccc = find(cnx_list(:,1)==rec_num & cnx_list(:,2)==cnx);

load([res_dir prefix '_detected_final'])


%%
load('Z:\data\fluctuation_project_fastslow\761418226_fastslow_014.mat')
load('E:\761418226_fastslow_014.mat')

%%

% refit glm to get errorbars on stp effects
W_short = repmat(fit.wt_long.*fit.Xc, 1, fit.stp_Nq).*fit.stp_X;
offset = fit.beta0 + log(data.dt) + fit.hist*fit.hist_beta;
[~,~,stat] = glmfit(double(W_short),double(data.post_spk_vec(~data.tin_nan))','poisson','Offset',offset,'constant','off');
fit.covB=stat.covb;
clear W_short offset

isi_t = 1:min(10000,size(fit.stp_basis,2));
% modification function
se_mod_fn = sqrt(diag(fit.stp_basis(:,isi_t)'*fit.covB*fit.stp_basis(:,isi_t)));
figure(11)
subplot(2,1,1)
semilogx(1 + [fit.stp_basis(:,isi_t)'*fit.wt_short_param], 'r', 'LineWidth', 2)
hold on
plot(1 + fit.stp_basis(:,isi_t)'*fit.wt_short_param + se_mod_fn,'r:')
plot(1 + fit.stp_basis(:,isi_t)'*fit.wt_short_param - se_mod_fn,'r:')
hold off
yline(1, 'k--', 'LineWidth', 2);
box off
xl=xlim();
set(gca,'TickDir','out')

% histogram for ISI
subplot(2,1,2)
histogram(log10(diff(data.pre_spk_times/data.dt)),'EdgeColor','none')
xlim(log10(xl))
% xlim([log10(2) log10(ep)])
% xticks(log10(tickMark))
% xticklabels(arrayfun(@num2str, tickMark, 'UniformOutput', 0))
set(gca,'TickDir','out')
box off

%% match optimal window size for efficacy and get short/long effects



t = (1:length(data.tin_nan))*data.dt;
lam_full=t*NaN;
lam_nocnx=t*NaN;
lam_full(~data.tin_nan) = exp(fit.beta0 +...
    fit.wt_long.*fit.wt_short.*fit.Xc +...
    fit.hist*fit.hist_beta)*data.dt;
lam_nocnx(~data.tin_nan) = exp(fit.beta0 +fit.hist*fit.hist_beta)*data.dt;

eff_detail=[];
wt_s=t*NaN;
wt_l=t*NaN;
wt_s(~data.tin_nan)=fit.wt_short;
wt_l(~data.tin_nan)=fit.wt_long;
wt_swin=[];
wt_lwin=[];
w=smooth_estimate_all(ccc).window;
for i=1:length(smooth_estimate_all(ccc).tcenter)
    wid=t>(smooth_estimate_all(ccc).tcenter(i)-w/2) & t<(smooth_estimate_all(ccc).tcenter(i)+w/2);
    eff_detail(i)=nansum(lam_full(wid)-lam_nocnx(wid))./sum(data.pre_spk_vec(wid));
    wt_swin(i) = nanmean(wt_s(wid));
    wt_lwin(i) = nanmean(wt_l(wid));
end
%%

fpp_hi = rate_estimate_all(ccc).pre_fr.*rate_estimate_all(ccc).pos_fr;
[fpp,ty] = resample(fpp_hi,rate_estimate_all(ccc).tin,0.1,1,100,'linear');
fpp = interp1(ty,fpp,smooth_estimate_all(ccc).tcenter,'linear')';


fp2p_hi = rate_estimate_all(ccc).pre_fr.*rate_estimate_all(ccc).pos_fr.^2;
[fp2p,ty] = resample(fp2p_hi,rate_estimate_all(ccc).tin,0.1,1,100,'linear');
fp2p = interp1(ty,fp2p,smooth_estimate_all(ccc).tcenter,'linear')';

tau = 60;
twin = linspace(0,tau*10,tau*10/0.1+1);
wfun = twin/tau.*exp(1-twin/tau);
wfun = wfun/sum(wfun);

fpp_hi_filt = conv(fpp_hi,wfun);
fpp_hi_filt = fpp_hi_filt(1:length(fpp_hi));
fpp_filt = interp1(rate_estimate_all(ccc).tin,fpp_hi_filt,smooth_estimate_all(ccc).tcenter,'linear')';

fp2p_hi_filt = conv(fp2p_hi,wfun);
fp2p_hi_filt = fp2p_hi_filt(1:length(fp2p_hi));
fp2p_filt = interp1(rate_estimate_all(ccc).tin,fp2p_hi_filt,smooth_estimate_all(ccc).tcenter,'linear')';

fpp_filt = smoothdata(fpp_filt,"movmean",smooth_estimate_all(ccc).window/10);
fp2p_filt = smoothdata(fp2p_filt,"movmean",smooth_estimate_all(ccc).window/10);

%%
figure(10)
subplot(4,1,1)
shadedErrorBar(smooth_estimate_all(ccc).tcenter/60,smooth_estimate_all(ccc).eff,smooth_estimate_all(ccc).effse)
line(xlim(),[0 0])
axis tight; box off;
set(gca,'TickDir','out')
xlabel('Time [min]')
ylabel('Efficacy')
hold on
plot(smooth_estimate_all(ccc).tcenter/60,eff_detail)
hold off
subplot(4,1,2)
plot(smooth_estimate_all(ccc).tcenter/60,wt_swin)
hold on
plot(smooth_estimate_all(ccc).tcenter/60,wt_lwin)
hold off

subplot(4,1,3)
plot(rate_estimate_all(ccc).tin/60,rate_estimate_all(ccc).pre_fr)
hold on
plot(rate_estimate_all(ccc).tin/60,rate_estimate_all(ccc).pos_fr)
hold off

subplot(4,1,4)
plot(fpp_filt)
hold on
plot(fp2p_filt)
hold off

b=glmfit([fpp_filt fp2p_filt],wt_lwin');
wt_lwin_hat = b(1)+[fpp_filt fp2p_filt]*b(2:3);