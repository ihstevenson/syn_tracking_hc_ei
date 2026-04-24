c=1;
aic_stp=[]; aic_cov=[]; aic_bcm=[];
r2_stp=[]; r2_cov=[]; r2_bcm=[];
tau_vec = logspace(log10(10),log10(360),11);
glist = dir([output_dir '*_detection.mat']);
for rec_num=1:length(glist)
    prefix = strsplit(glist(rec_num).name,'_');
    prefix = prefix{1};
    load([output_dir prefix '_detection'])
    load([output_dir prefix '_opt_window'])

    T = max(cellfun(@max,data.Tlist));

    if ~isempty(detection.ijlist)
        for cnx = 1:size(detection.ijlist,1)
            ij = detection.ijlist(cnx,:);

            % optimal bandwidth rate estimates at dt resolution
            dt=0.1; % time resolution for rate effects [s]
            tin = linspace(0,T,ceil(T/dt));
            [y_pre,t_pre,optw,W1,C1] = sskernel(data.Tlist{ij(1)},tin);
            pre_fr = y_pre/sum(y_pre)*length(data.Tlist{ij(1)});
            [y_post,t_post,optw2,W12,C12] = sskernel(data.Tlist{ij(2)},tin);
            pos_fr = y_post/sum(y_post)*length(data.Tlist{ij(2)});

            % get (anti-aliased) rate, that matches efficacy windows
            % smooth_estimate using smooth_estimate.window windows w/ 10 stride
            tid = interp1(tin,1:length(tin),smooth_estimate(cnx).tcenter,'linear');
            [fpre,ty] = resample(pre_fr,tin,0.1,1,100,'linear');
            fpre10 = interp1(ty,fpre,smooth_estimate(cnx).tcenter,'linear')';
            [fpos,ty] = resample(pos_fr,tin,0.1,1,100,'linear');
            fpos10 = interp1(ty,fpos,smooth_estimate(cnx).tcenter,'linear')';
            fpre = smoothdata(fpre10,"movmean",smooth_estimate(cnx).window/10);
            fpos = smoothdata(fpos10,"movmean",smooth_estimate(cnx).window/10);

            % fit STP model
            y = smooth_estimate(cnx).eff;
            yse = smooth_estimate(cnx).effse.^2;
            X = zscore(real([log(fpre) log(fpos)]));
            validid=all(X>-10,2) & all(isfinite(X),2) & isfinite(y);
            [b,dev,stat] = glmfit(X(validid,:),y(validid),'normal','weights',1./yse(validid));
            effhat = glmval(b,X,'identity');
            aic_stp(c) = 2*length(b)+dev;

            % fit STP+BCM models
            % get covariance and 3rd order terms from high-resolution rate estimates
            fpp_hi = pre_fr.*pos_fr;
            fp2p_hi = pre_fr.*pos_fr.^2;

            % filter and downsample to resolution of efficacy observations
            for taui = 1:length(tau_vec)
                tau = tau_vec(taui);
                twin = linspace(0,tau*10,tau*10/0.1+1);
                wfun = twin/tau.*exp(1-twin/tau);
                wfun = wfun/sum(wfun);
                fpp_hi_filt = conv(fpp_hi,wfun);
                fpp_hi_filt = fpp_hi_filt(1:length(fpp_hi));
                fpp_filt = interp1(tin,fpp_hi_filt,smooth_estimate(cnx).tcenter,'linear')';
                fp2p_hi_filt = conv(fp2p_hi,wfun);
                fp2p_hi_filt = fp2p_hi_filt(1:length(fp2p_hi));
                fp2p_filt = interp1(tin,fp2p_hi_filt,smooth_estimate(cnx).tcenter,'linear')';
                fpp_filt = smoothdata(fpp_filt,"movmean",smooth_estimate(cnx).window/10);
                fp2p_filt = smoothdata(fp2p_filt,"movmean",smooth_estimate(cnx).window/10);

                % STP+covariance
                X = zscore(real([log(fpre) log(fpos) fpp_filt]));
                [b2,dev,stat] = glmfit(X(validid,:),y(validid),'normal','weights',1./yse(validid));
                effhat(:,2) = glmval(b2,X,'identity');
                aic_cov(c,taui) = 2*length(b2)+dev;

                % STP+BCM
                X = zscore(real([log(fpre) log(fpos) fpp_filt fp2p_filt]));
                [b3,dev,stat] = glmfit(X(validid,:),y(validid),'normal','weights',1./yse(validid));
                effhat(:,3) = glmval(b3,X,'identity');
                aic_bcm(c,taui) = 2*length(b3)+dev;

                % goodness-of-fit
                tmp = corrcoef([y(validid) effhat(validid,:)]);
                r2 = tmp(1,2:end).^2;
                r2_stp(c) = r2(1);
                r2_cov(c,taui) = r2(2);
                r2_bcm(c,taui) = r2(3);
            end
            c=c+1;
        end
    end
end

save('batch_linear_result','aic_stp','aic_cov','aic_bcm','r2_stp','r2_cov','r2_bcm')
%%

