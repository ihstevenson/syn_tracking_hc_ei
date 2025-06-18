
% Process Allen Institute - Visual Coding Neuropixels recordings
% to detect putative synapses using an initial hypothesis test, followed
% by a convolutional model.

function [data,detection] = syn_detection(dat_dir,prefix)

load([dat_dir prefix '_spikes.mat'])
tab = readtable([dat_dir prefix '_units.csv']);
data.invalid = readtable([dat_dir prefix '_invalid.csv']);

% add sub-resolution jitter to prevent aliasing in the ccg
data.Tlist = cellfun(@(x) (x+(rand(size(x))-.5)*0.0001)',Tlist,'UniformOutput',false)';

data.anatomy = tab.ecephys_structure_acronym;
data.xorig = tab.probe_horizontal_position;
data.y = tab.probe_vertical_position;
data.probe = tab.probe_id;
[~,~,data.probeu]=unique(data.probe);
data.x = data.xorig+data.probeu*500-35;

data.area_list = {'CA1';'CA3';'DG';'SUB';'ProS'};
for area=1:length(data.area_list)
    data.area_mat(:,area) = strcmp(data.anatomy,data.area_list{area});
end

% Limit neurons so those with high SNR, >1k spikes, in relevant area
data.ids = tab.snr>2 & cellfun(@length,data.Tlist)>1000 & any(data.area_mat')';

data.Tlist=data.Tlist(data.ids);
data.tab_sub = tab(data.ids,:);

if isempty(Tlist)
    % no units matched the criteria
    return
end

% run hypothesis tests
nhst_params=[];
nhst_params.trange = [-0.025 0.025];
nhst_params.tn = 102;
nhst_params.t = linspace(nhst_params.trange(1),nhst_params.trange(2),nhst_params.tn);
nhst_params.t = nhst_params.t(1:end-1)+mean(diff(nhst_params.t))/2;
nhst_params.jit_min_t = 0.0011;
nhst_params.ssi_min_t = 0.00049;
nhst_params.ssi_max = 0.4;
nhst_params.alph = 0.1;
nhst_params.jit_conv_window = 2;
nhst_params.mask=abs(nhst_params.t)>0.0000001;


%
detection.result = nhst_detection(data.Tlist,nhst_params);

% fit pairs passing hypothesis test
[r,c] = find(detection.result.detected_cnx);

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

% cut-off at 200k spikes
Rlist=cellfun(@(x) x(1:min(length(x),200000)),data.Tlist,'UniformOutput',false);

figure(2)
clf
clear fitout
for i=1:length(r)
    fprintf('%03i/%03i...\n',i,length(r))
    fitout(i) = composite_model_fit([r(i) c(i)],Rlist,params,'i>j');
    drawnow

    % repair fits with nonphysiological effects...
    if fitout(i).eff<-1
        rparams = params;
        rparams.nknots=2;
        fitout2 = composite_model_fit([r(i) c(i)],Rlist,rparams,'i>j');
        fitout(i) = composite_model_fit([r(i) c(i)],Rlist,params,'i>j',[],fitout2.feat_params);
    end
end

% screen
fit_res=[];
% pseudo-Rsquared
rsq=[];
for i=1:length(fitout)
    fit_res.cvsmoo(i) = std(fitout(i).lamsmoo)/mean(fitout(i).lamsmoo);
    N = length(Tlist{r(i)});
    y = fitout(i).yhat(params.mask)*N+fitout(i).stat.resid;
    yhat = fitout(i).yhat(params.mask);
    llm = sum(y.*log(yhat)+(N-y).*log(1-yhat));
    ll0 = sum(y.*log(mean(y)/N)+(N-y).*log(1-mean(y)/N));
    lls = nansum(y.*log(y./N)+(N-y).*log(1-y./N));
    fit_res.rsq(i)=(llm-ll0)./(lls-ll0);
end

% latency and duration [ms] and efficacy
fit_res.lat = exp(arrayfun(@(x) x.feat_params(1),fitout))*1000;
fit_res.dur = exp(arrayfun(@(x) x.feat_params(2),fitout))*1000;
fit_res.eff = arrayfun(@(x) x.eff,fitout);
fit_res.ex = arrayfun(@(x) x.eff>0,fitout);
fit_res.in = arrayfun(@(x) x.eff<0,fitout);
fit_res.ssi = detection.result.ssi(sub2ind(size(detection.result.ssi),r,c))';
    
% select subset of candidate pairs
fit_res.kidx = (fit_res.cvsmoo<.2) & (fit_res.rsq>.7) & (fit_res.dur<1) & (fit_res.lat>.1) & (fit_res.lat<10) & (abs(fit_res.eff)<.5) & (abs(fit_res.eff)>.01) & fit_res.ex;
detection.ijlist = [r(fit_res.kidx) c(fit_res.kidx)];
detection.fidx = find(fit_res.kidx);

detection.nhst_params = nhst_params;
detection.params = params;
detection.fitout = fitout;
detection.fit_res = fit_res;