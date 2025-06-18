
res_dir = 'C:\Users\istevenson\Documents\GitHub\fluctuation_project\detection\';
dwdir = 'E:\fluctuation_backup\dw\';
glist = dir([res_dir '*_detected_final.mat']);

% Fig 3A [14-14]
% Fig 3B [25-9]

% [35-5, 25-13]

% Fig 3C [1-6,29-3,25-28]
% Fig 3D [51-10,25-1,40-10] 
% Fig 3E [2-27,32-9,34-9]
tmp=strsplit(glist(14).name,'_');
prefix = tmp{1};
cnx=14;


load([res_dir prefix '_detected_final'])
[r,c] = find(result.detected_cnx);
fit_res.ssi = result.ssi(sub2ind(size(result.ssi),r,c))';
fit_res.kidx = (fit_res.cvsmoo<.1) & fit_res.ssi<0.2 & (fit_res.rsq>.5) & (fit_res.dur<1) & (fit_res.lat>0.75) & (fit_res.lat<10) & (abs(fit_res.eff)<.5) & (abs(fit_res.eff)>.01) & fit_res.ex;
ijlist = [r(fit_res.kidx) c(fit_res.kidx)];
T = max(cellfun(@max,Tlist));
load([dwdir prefix '_joint_res_300.mat'])


% rates
pre = ijlist(cnx,1);
post = ijlist(cnx,2);
fidx = find(r==pre & c==post);

dt=1e-1;
tin = linspace(0,T,ceil(T/dt));
[y_pre,t_pre,optw,W1,C1] = sskernel(Tlist{pre},tin);
pre_fr = y_pre/sum(y_pre)*length(Tlist{pre});
[y_post,t_post,optw2,W12,C12] = sskernel(Tlist{post},tin);
pos_fr = y_post/sum(y_post)*length(Tlist{post});

% joint rate statistics
pedges = linspace(1,length(pre_fr),10);
pre_rank = tiedrank(pre_fr);
pos_rank = tiedrank(pos_fr);
[N,Xedges,Yedges,binX,binY] = histcounts2(pre_rank,pos_rank,pedges,pedges);


% xcorr weighted following joint rate conditions
[c,cmat] = spk_xcorr(Tlist{pre},Tlist{post},params.interval(1),params.interval(2),params.num_bins);
cedges = linspace(params.interval(1),params.interval(2),params.num_bins+1);
[~,~,cbins] = histcounts(cmat(:,1),cedges);
C = zeros(length(cbins),max(cbins));
for i=1:max(cbins)
    C(:,i) = cbins==i;
end

nsplit=10;
clear axjoint
dfac = length(pre_rank)/nsplit;
pedges = (1:(nsplit)+1)-.5;
[Nlo,~,~,binXlo,binYlo] = histcounts2(ceil(pre_rank/dfac),ceil(pos_rank/dfac),pedges,pedges);

pax = linspace(0,100,nsplit+1);
pax = pax(1:end-1)+diff(pax)/2;
pre_ax = prctile(pre_fr,pax)/dt;
pos_ax = prctile(pos_fr,pax)/dt;

% alpha fn weights
tau = 300;
twin = linspace(0,tau*10,tau*10/dt+1);
wfun = twin/tau.*exp(1-twin/tau);
wfun = wfun/sum(wfun);
% 
% figure(51)
% clf
% subplot(5,5,1)
% imagesc(N)
% set(gca,'TickDir','out'); set(gca,'YDir','normal'); box off; axis equal
% axis image
% ylabel('Pre')
% xlabel('Post')
% 
% subplot(5,5,3:5)
% plot(tin/60,pre_fr/dt,tin/60,pos_fr/dt)
% box off; set(gca,'TickDir','out')
% axis tight
% 
% c=1;
% for i=[1 10]
%     for j=[1 10]
%         fprintf('.')
%         if sum(binXlo==i & binYlo==j)>0
%             pattern_ids = binXlo==i & binYlo==j;
%             pattern_mask = conv(pattern_ids,wfun);
%             pattern_mask = pattern_mask(1:length(tin));
%             pre_w = interp1(tin,pattern_mask,Tlist{pre},'nearest')';
% 
%             figure(51)
%             subplot(5,5,c*5+1)
%             fitjoint_post(i,j) = composite_model_fit_weighted([pre post],Tlist,params,'i>j',pre_w,fitout(fidx).feat_params);
%             yl=ylim();
% 
%             pattern_mask2 = conv(fliplr(pattern_ids),wfun);
%             pattern_mask2 = fliplr(pattern_mask2(1:length(tin)));
%             pre_w2 = interp1(tin,pattern_mask2,Tlist{pre},'nearest')';
% 
%             subplot(5,5,c*5+2)
%             fitjoint_pre(i,j) = composite_model_fit_weighted([pre post],Tlist,params,'i>j',pre_w2,fitout(fidx).feat_params);
%             yl2=ylim();
% 
%             yl = [0 max(yl(2),yl2(2))];
%             ylim(yl);  box off; set(gca,'TickDir','out')
%             subplot(5,5,c*5+1)
%             ylim(yl); box off; set(gca,'TickDir','out')
% 
% 
%             subplot(5,5,c*5+(3:5))
%             plot(tin/60,pattern_mask/max(pattern_mask)/2,tin/60,pattern_mask2/max(pattern_mask2)/2)
%             hold on
%             hedges=linspace(0,T/60,256);
%             histogram(find(pattern_ids)*dt/60,hedges,'Normalization','pdf','EdgeColor','none')
%             hold off
%             box off; set(gca,'TickDir','out')
%             axis tight
%         end
%         c=c+1;
%     end
% 
%     fprintf('\n')
% end

% dw map and bcm fit
showPrePostDw(joint_res(cnx).fitjoint_pre,joint_res(cnx).fitjoint_post,joint_res(cnx).Nlo)
figure(24)
subplot(1,4,1)
set(gca,'CLim',[0 max(get(gca,'CLim'))])
subplot(1,4,2)
cl(1,:)=get(gca,'CLim');
subplot(1,4,3)
cl(2,:)=get(gca,'CLim');
set(gca,'CLim',[min(cl(:,1)) max(cl(:,2))])
subplot(1,4,2)
set(gca,'CLim',[min(cl(:,1)) max(cl(:,2))])

% bcm fit
cent = linspace(0,100,11);
cent = cent(1:end-1)+diff(cent)/2;
pre_pct = prctile(pre_fr,cent);
pos_pct = prctile(pos_fr,cent);

w=[];
for i=1:size(joint_res(cnx).fitjoint_post,1)
    for j=1:size(joint_res(cnx).fitjoint_post,2)
        w(i,j,1) = joint_res(cnx).fitjoint_post(i,j).eff;
        w(i,j,2) = joint_res(cnx).fitjoint_pre(i,j).eff;
    end
end
dw = (w(:,:,1)-w(:,:,2))*100;



[rpos,rpre] = meshgrid(pos_pct,pre_pct);
% X = [rpos(:).^2.*rpre(:) -rpos(:).*rpre(:)];
X = [rpos(:).^2.*rpre(:) -rpos(:).*rpre(:) rpre(:) rpos(:)];
% X = [rpos(:) rpre(:) rpos(:).*rpre(:)];
b = glmfit(X,dw(:),'normal');
dwhat = glmval(b,X,'identity');

% b = glmfit(X,dw(:),'normal','constant','off');
% dwhat = glmval(b,X,'identity','constant','off');
dwhat = reshape(dwhat,10,10);
figure(25)
subplot(2,2,3)
imagesc(dw)
set(gca,'YDir','normal')
subplot(2,2,4)
imagesc(dwhat)
set(gca,'YDir','normal')

subplot(2,2,3)
cl(1,:)=get(gca,'CLim');
subplot(2,2,4)
cl(2,:)=get(gca,'CLim');
set(gca,'CLim',[min(cl(:,1)) max(cl(:,2))])
subplot(2,2,3)
set(gca,'CLim',[min(cl(:,1)) max(cl(:,2))])

x0=[1:10]-.5;
figure(25)
ax(1)=subplot(2,2,1)
hold on
plot(x0,mean(dwhat(1:3,:)),':')
plot(x0,mean(dwhat(4:6,:)),':')
plot(x0,mean(dwhat(7:10,:)),':')
hold off
xlabel('Post Rate [%]')
box off; set(gca,'TickDir','out')
ax(2)=subplot(2,2,2)
hold on
plot(x0,mean(dwhat(:,1:3),2),':')
plot(x0,mean(dwhat(:,4:6),2),':')
plot(x0,mean(dwhat(:,7:10),2),':')
hold off
xlabel('Pre Rate [%]')
box off; set(gca,'TickDir','out')

linkaxes(ax)
subplot(2,2,3)
cl=get(gca,'CLim');
set(gca,'Clim',[-1 1]*max(abs(cl)));
axis equal tight
colorbar
subplot(2,2,4)
set(gca,'Clim',[-1 1]*max(abs(cl)));
colorbar
axis equal tight
% 
% % %%
% % figure(51)
% % exportgraphics(gcf,[prefix '_' sprintf('%03i',cnx) '_before_after_detail.pdf'],'ContentType','vector')
% figure(24)
% colormap parula
% exportgraphics(gcf,[prefix sprintf('_%03i_dw.pdf',cnx)],'ContentType','vector')
% figure(25)
% colormap turbo
% exportgraphics(gcf,[prefix sprintf('_%03i_dw_fit_v3.pdf',cnx)],'ContentType','vector')
