
res_dir = 'C:\Users\istevenson\Documents\GitHub\fluctuation_project\detection\';
glist = dir([res_dir '*_detected_final.mat']);


ccc=1;
winsize=[]; cv=[];
% smooth_estimate_all=[]; rate_estimate_all=[]; cnx_list=[];
for rec_num=1:length(glist)
    prefix = strsplit(glist(rec_num).name,'_');
    prefix = prefix{1};

    load([res_dir prefix '_detected_final'])

    [r,c] = find(result.detected_cnx);
    fit_res.ssi = result.ssi(sub2ind(size(result.ssi),r,c))';
    fit_res.kidx = (fit_res.cvsmoo<.1) & fit_res.ssi<0.2 & (fit_res.rsq>.5) & (fit_res.dur<1) & (fit_res.lat>0.75) & (fit_res.lat<10) & (abs(fit_res.eff)<.5) & (abs(fit_res.eff)>.01) & fit_res.ex;
    ijlist = [r(fit_res.kidx) c(fit_res.kidx)];
    T = max(cellfun(@max,Tlist));

    fprintf('%03i...\n',rec_num)
    if isempty(ijlist)
        continue;
    end
    load([res_dir prefix '_opt_window'])

    % if ccc==1,
    %     smooth_estimate_all = smooth_estimate(1);
    % end

    for cnx=1:size(ijlist,1)
        cnx_list(ccc,:)=[rec_num cnx];
        winsize(ccc)=smooth_estimate(cnx).window;
        cv(ccc)=std(smooth_estimate(cnx).eff)/mean(smooth_estimate(cnx).eff);

        pre = ijlist(cnx,1);
        post = ijlist(cnx,2);
        fidx = find(r==pre & c==post);

        overall_eff(ccc,1)=fitout(fidx).eff;
        overall_eff(ccc,2)=fitout(fidx).effse;
        cnx_anat(ccc,1) = tab_sub.ecephys_structure_acronym(pre);
        cnx_anat(ccc,2) = tab_sub.ecephys_structure_acronym(post);
        overall_fr(ccc,1) = length(Tlist{pre})/T;
        overall_fr(ccc,2) = length(Tlist{post})/T;
        overall_fr(ccc,3) = tab_sub.waveform_duration(pre);
        overall_fr(ccc,4) = tab_sub.waveform_duration(post);

        % dt=1e-1;
        % tin = linspace(0,T,ceil(T/dt));
        % [y_pre,t_pre,optw,W1,C1] = sskernel(Tlist{pre},tin);
        % pre_fr = y_pre/sum(y_pre)*length(Tlist{pre});
        % [y_post,t_post,optw2,W12,C12] = sskernel(Tlist{post},tin);
        % pos_fr = y_post/sum(y_post)*length(Tlist{post});
        %
        % smooth_estimate_all(ccc) = smooth_estimate(cnx);
        % rate_estimate_all(ccc).tin = tin;
        % rate_estimate_all(ccc).pre_fr = pre_fr;
        % rate_estimate_all(ccc).pos_fr = pos_fr;

        ccc=ccc+1;
    end
end

% save('all_smooth_estimate','cnx_list','smooth_estimate_all','rate_estimate_all')
save('all_smooth_estimate_extra','cv','winsize','overall_eff','cnx_anat','overall_fr')

%%

doPlot=true;
% coeffs=cell(0);
% for ccc=1:length(rate_estimate_all)
    % for ccc=252 % 27-2
    for ccc=174 % 19-14 - use for fig 4B
    % for ccc=100 % 13-7 - use for fig 4A
    % for ccc=90 % 12-3
    % for ccc=400

    % get (anti-aliased) rate, that matches efficacy windows
    tid = interp1(rate_estimate_all(ccc).tin,1:length(rate_estimate_all(ccc).tin),smooth_estimate_all(ccc).tcenter,'linear');
    [fpre,ty] = resample(rate_estimate_all(ccc).pre_fr,rate_estimate_all(ccc).tin,0.1,1,100,'linear');
    fpre = interp1(ty,fpre,smooth_estimate_all(ccc).tcenter,'linear')';
    [fpos,ty] = resample(rate_estimate_all(ccc).pos_fr,rate_estimate_all(ccc).tin,0.1,1,100,'linear');
    fpos = interp1(ty,fpos,smooth_estimate_all(ccc).tcenter,'linear')';

    fpre10=fpre;
    fpos10=fpos;
    fpre = smoothdata(fpre,"movmean",smooth_estimate_all(ccc).window/10);
    fpos = smoothdata(fpos,"movmean",smooth_estimate_all(ccc).window/10);

    % naixin's jneuro model
    y = smooth_estimate_all(ccc).eff;
    yse = smooth_estimate_all(ccc).effse.^2;
    % X = [fpre fpos];
    % validid=all(X>10e-6,2);

    X = zscore(real([log(fpre) log(fpos)]));
    validid=all(X>-10,2) & all(isfinite(X),2) & isfinite(y);
    % [b,dev,stat] = glmfit(X(validid,:),y(validid),'normal');
    [b,dev,stat] = glmfit(X(validid,:),y(validid),'normal','weights',1./yse(validid));
    effhat = glmval(b,X,'identity');

    aic(ccc,1) = 2*length(b)+dev;

    % LTP/LTD model
    fpp_hi = rate_estimate_all(ccc).pre_fr.*rate_estimate_all(ccc).pos_fr;
    % fpp_hi = zscore(rate_estimate_all(ccc).pre_fr).*zscore(rate_estimate_all(ccc).pos_fr);
    % fpp_hi = log(rate_estimate_all(ccc).pre_fr).*log(rate_estimate_all(ccc).pos_fr);
    [fpp,ty] = resample(fpp_hi,rate_estimate_all(ccc).tin,0.1,1,100,'linear');
    fpp = interp1(ty,fpp,smooth_estimate_all(ccc).tcenter,'linear')';

    fp2p_hi = rate_estimate_all(ccc).pre_fr.*rate_estimate_all(ccc).pos_fr.^2;
    % fp2p_hi = zscore(rate_estimate_all(ccc).pre_fr).*zscore(rate_estimate_all(ccc).pos_fr).^2;
    % fp2p_hi = log(rate_estimate_all(ccc).pre_fr).*log(rate_estimate_all(ccc).pos_fr).^2;
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




    X = zscore(real([log(fpre) log(fpos) fpp_filt]));
    % [b2,dev,stat] = glmfit(X(validid,:),y(validid),'normal');
    [b2,dev,stat] = glmfit(X(validid,:),y(validid),'normal','weights',1./yse(validid));
    effhat(:,2) = glmval(b2,X,'identity');
    aic(ccc,2) = 2*length(b2)+dev;

    X = zscore(real([log(fpre) log(fpos) fpp_filt fp2p_filt]));
    % [b3,dev,stat] = glmfit(X(validid,:),y(validid),'normal');
    [b3,dev,stat] = glmfit(X(validid,:),y(validid),'normal','weights',1./yse(validid));
    effhat(:,3) = glmval(b3,X,'identity');
    aic(ccc,3) = 2*length(b3)+dev;
    % r2(ccc,:)=1-var(y(validid)-effhat(validid,:))./var(y(validid));
    tmp = corrcoef([y(validid) effhat(validid,:)]);
    r2(ccc,:) = tmp(1,2:end).^2;

    zzm = mean([fpp_filt fp2p_filt]);
    zzs = std([fpp_filt fp2p_filt]);

    coeffs{1}(:,ccc)=b;
    coeffs{2}(:,ccc)=b2;
    coeffs{3}(:,ccc)=b3;

    avg_eff(ccc) = mean(y(validid));

    if doPlot
        r2(ccc,:)
        cnx_anat(ccc,:)
        figure(28)
        clf
        subplot(2,1,1)
        shadedErrorBar(smooth_estimate_all(ccc).tcenter/60,smooth_estimate_all(ccc).eff,smooth_estimate_all(ccc).effse)
        line(xlim(),[0 0])
        axis tight; box off;
        set(gca,'TickDir','out')
        xlabel('Time [min]')
        ylabel('Efficacy')
        hold on
        plot(smooth_estimate_all(ccc).tcenter/60,effhat)
        hold off
        subplot(2,1,2)
        plot(smooth_estimate_all(ccc).tcenter/60,zscore(fpre10))
        hold on
        plot(smooth_estimate_all(ccc).tcenter/60,zscore(fpos10))
        plot(smooth_estimate_all(ccc).tcenter/60,zscore(fpp_filt))
        % plot(smooth_estimate_all(ccc).tcenter/60,zscore(fp2p_filt))
        hold off
    end
end

% exportgraphics(gcf,sprintf('%03i_%03i_dyn_lin_fit.pdf',cnx_list(ccc,1),cnx_list(ccc,2)),'ContentType','vector')

%% rule

% after running fig3 on same pair
figure(8)
imagesc(b3(4)*(rpre.*rpos-zzm(1))./zzs(1)+b3(5)*(rpre.*rpos.^2-zzm(2))/zzs(2))
set(gca,'YDir','normal')
cl=get(gca,'CLim')
set(gca,'CLim',[-1 1]*max(abs(cl)))
colormap turbo
colorbar
% exportgraphics(gcf,[prefix sprintf('_%03i_dw_dyn_lin.pdf',cnx)],'ContentType','vector')

%%

figure(29)
% plot(coeffs{1}(2,:),coeffs{1}(3,:),'.')
% plot(coeffs{2}(2,:),coeffs{2}(4,:),'.')
% scatter(coeffs{2}(3,:),coeffs{2}(4,:),r2(:,2)*100)
% scatter(coeffs{2}(2,:),coeffs{2}(4,:),r2(:,2)*100)
scatter(coeffs{3}(2,:),coeffs{3}(5,:),r2(:,3)*100)

sz=max(overall_eff(:,1),10e-12)*500;
% sz=r2(:,3)*200;
% sz = overall_fr(:,2)*2;
% sz = overall_fr(:,4)*200;
[uf,~,ua] = unique(cnx_anat(:,1));

cmap=lines(5);

yl=[-0.5 2];
subplot(3,3,1)
scatter(coeffs{3}(3,:),coeffs{3}(2,:),sz,cmap(ua,:),'filled','MarkerFaceAlpha',.5)
ylim(yl)
line([0 0],ylim())
line(xlim(),[0 0])
subplot(3,3,2)
scatter(coeffs{3}(4,:),coeffs{3}(2,:),sz,cmap(ua,:),'filled','MarkerFaceAlpha',.5)
ylim(yl)
line(xlim(),[0 0])
line([0 0],ylim())
subplot(3,3,3)
scatter(coeffs{3}(5,:),coeffs{3}(2,:),sz,cmap(ua,:),'filled','MarkerFaceAlpha',.5)
ylim(yl)
line(xlim(),[0 0])
line([0 0],ylim())

subplot(3,3,5)
scatter(coeffs{3}(4,:),coeffs{3}(3,:),sz)
subplot(3,3,6)
scatter(coeffs{3}(5,:),coeffs{3}(3,:),sz)

subplot(3,3,9)
scatter(coeffs{3}(5,:),coeffs{3}(4,:),sz,cmap(ua,:),'filled','MarkerFaceAlpha',.5)

% subplot(3,3,7)
% scatter(overall_fr(:,1),coeffs{3}(4,:)-coeffs{3}(5,:),sz,cmap(ua,:),'filled','MarkerFaceAlpha',.5)
% subplot(3,3,8)
% scatter(overall_fr(:,2),coeffs{3}(4,:)-coeffs{3}(5,:),sz,cmap(ua,:),'filled','MarkerFaceAlpha',.5)

%%
figure(30)
for i=1:5
    subplot(5,1,i)
    % plot(overall_eff(:,1),coeffs{3}(i,:),'o')
    plot(overall_fr(:,4),coeffs{3}(i,:),'o')
end

figure(30)
clf
plot(overall_fr(:,4),coeffs{3}(5,:)./coeffs{3}(4,:),'o')

figure(31)
clf
for i=1:5
    x=coeffs{3}(i,:);
    % x=coeffs{2}(i,:);
    % x=coeffs{1}(i,:);
    x(abs(x)>1)=NaN;
    edges = linspace(min(x),max(x),1024);
    for j=1:3
        subplot(2,5,i)
        % histogram(coeffs{3}(i,ua==j))
        [f,xi]= ksdensity(x(ua==j),edges);
        f=f/sum(f);
        plot(xi,f)
        hold on
    end
    hold off
    % xlim([-1 1]*max(abs(xlim()))/2)
    xlim([-1 1]*.05)
    box off; set(gca,'TickDir','out')
    subplot(2,5,i+5)
    boxplot(x,ua,'Orientation','horizontal','Whisker',100)
    xlim([-1 1]*.05)
    line([0 0],ylim())
    set(gca,'YTickLabel',uf)
end
% exportgraphics(gcf,'fig4_dyn_lin_coefs.pdf','ContentType','vector')

%%

figure(4)
clf
for m=1:3
    p0=linspace(0,100,512);
    plot(prctile(r2(:,m),p0),p0)
    hold on
    set(gca,'TickDir','out')
    box off
    xlabel('Variance Explained')
    ylabel('CDF')
end
hold off

%%

f=find(cnx_list(:,1)==25);
effmat=[];
for i=1:length(f)
    tcomm = linspace(0,max(smooth_estimate_all(f(i)).tcenter)+600,1024);
    validid=isfinite(smooth_estimate_all(f(i)).eff);
    effmat(i,:) = interp1(smooth_estimate_all(f(i)).tcenter(validid),smooth_estimate_all(f(i)).eff(validid),tcomm,'nearest',true);
end
imagesc(zscore(effmat(:,30:900)')')
set(gca,'CLim',[-1 1]*4)
% imagesc(effmat(:,30:900)./mean(effmat(:,30:900),2))
% set(gca,'CLim',[-0.5 3])
colorbar

%%

tauvec = logspace(1,log10(600),20);

for tvj = 1:length(tauvec)
    doPlot=false;
    % coeffs=cell(0);
    for ccc=1:length(rate_estimate_all)
        % for ccc=252 % 27-2
        %for ccc=174 % 19-14 - use for fig 4B
        % for ccc=100 % 13-7 - use for fig 4A
        % for ccc=90 % 12-3

        % get (anti-aliased) rate, that matches efficacy windows
        tid = interp1(rate_estimate_all(ccc).tin,1:length(rate_estimate_all(ccc).tin),smooth_estimate_all(ccc).tcenter,'linear');
        [fpre,ty] = resample(rate_estimate_all(ccc).pre_fr,rate_estimate_all(ccc).tin,0.1,1,100,'linear');
        fpre = interp1(ty,fpre,smooth_estimate_all(ccc).tcenter,'linear')';
        [fpos,ty] = resample(rate_estimate_all(ccc).pos_fr,rate_estimate_all(ccc).tin,0.1,1,100,'linear');
        fpos = interp1(ty,fpos,smooth_estimate_all(ccc).tcenter,'linear')';

        fpre10=fpre;
        fpos10=fpos;
        fpre = smoothdata(fpre,"movmean",smooth_estimate_all(ccc).window/10);
        fpos = smoothdata(fpos,"movmean",smooth_estimate_all(ccc).window/10);

        % naixin's jneuro model
        y = smooth_estimate_all(ccc).eff;
        yse = smooth_estimate_all(ccc).effse.^2;
        % X = [fpre fpos];
        % validid=all(X>10e-6,2);

        X = zscore(real([log(fpre) log(fpos)]));
        validid=all(X>-10,2) & all(isfinite(X),2) & isfinite(y);
        % [b,dev,stat] = glmfit(X(validid,:),y(validid),'normal');
        [b,dev,stat] = glmfit(X(validid,:),y(validid),'normal','weights',1./yse(validid));
        effhat = glmval(b,X,'identity');

        aic(ccc,1) = 2*length(b)+dev;

        % LTP/LTD model
        fpp_hi = rate_estimate_all(ccc).pre_fr.*rate_estimate_all(ccc).pos_fr;
        % fpp_hi = zscore(rate_estimate_all(ccc).pre_fr).*zscore(rate_estimate_all(ccc).pos_fr);
        % fpp_hi = log(rate_estimate_all(ccc).pre_fr).*log(rate_estimate_all(ccc).pos_fr);
        [fpp,ty] = resample(fpp_hi,rate_estimate_all(ccc).tin,0.1,1,100,'linear');
        fpp = interp1(ty,fpp,smooth_estimate_all(ccc).tcenter,'linear')';

        fp2p_hi = rate_estimate_all(ccc).pre_fr.*rate_estimate_all(ccc).pos_fr.^2;
        % fp2p_hi = zscore(rate_estimate_all(ccc).pre_fr).*zscore(rate_estimate_all(ccc).pos_fr).^2;
        % fp2p_hi = log(rate_estimate_all(ccc).pre_fr).*log(rate_estimate_all(ccc).pos_fr).^2;
        [fp2p,ty] = resample(fp2p_hi,rate_estimate_all(ccc).tin,0.1,1,100,'linear');
        fp2p = interp1(ty,fp2p,smooth_estimate_all(ccc).tcenter,'linear')';

        tau = tauvec(tvj);
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




        X = zscore(real([log(fpre) log(fpos) fpp_filt]));
        % [b2,dev,stat] = glmfit(X(validid,:),y(validid),'normal');
        [b2,dev,stat] = glmfit(X(validid,:),y(validid),'normal','weights',1./yse(validid));
        effhat(:,2) = glmval(b2,X,'identity');
        aic(ccc,2) = 2*length(b2)+dev;

        X = zscore(real([log(fpre) log(fpos) fpp_filt fp2p_filt]));
        % [b3,dev,stat] = glmfit(X(validid,:),y(validid),'normal');
        [b3,dev,stat] = glmfit(X(validid,:),y(validid),'normal','weights',1./yse(validid));
        effhat(:,3) = glmval(b3,X,'identity');
        aic(ccc,3) = 2*length(b3)+dev;
        % r2(ccc,:)=1-var(y(validid)-effhat(validid,:))./var(y(validid));
        tmp = corrcoef([y(validid) effhat(validid,:)]);
        r2(ccc,:) = tmp(1,2:end).^2;

        zzm = mean([fpp_filt fp2p_filt]);
        zzs = std([fpp_filt fp2p_filt]);

        coeffs{1}(:,ccc)=b;
        coeffs{2}(:,ccc)=b2;
        coeffs{3}(:,ccc)=b3;

        avg_eff(ccc) = mean(y(validid));

        
    end
    fprintf('.')
    r2all(:,:,tvj)=r2;
end
fprintf('\n')

%%

semilogx(tauvec,squeeze(nanmean(r2all))'*100)
box off; set(gca,'TickDir','out')
xlabel('Integration Time [s]')
ylabel('Variance Explained [%]')
set(gca,'XTick',[10 30 90 270])
xlim([10 640])

%%

p0=[25 50 75];
bar(reshape([squeeze(nanmedian(r2all(:,:,1)))' squeeze(nanmedian(r2all(:,:,9)))' squeeze(nanmedian(r2all(:,:,18)))' squeeze(nanmedian(max(r2all,[],3)))'],1,[]))
P=[prctile(r2all(:,:,1),p0) prctile(r2all(:,:,9),p0) prctile(r2all(:,:,18),p0) prctile(max(r2all,[],3),p0)];
hold on
errorbar(1:12,P(2,:),P(2,:)-P(1,:),P(3,:)-P(2,:),'o','CapSize',0)
hold off
box off; set(gca,'TickDir','out')