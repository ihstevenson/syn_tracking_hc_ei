% get weighted correlograms and fits using the events following (tin,pattern_ids)

function fitw = spk_xcorr_fwdrev(Tlist,ij,tin,pattern_ids,syn_init,params,wfun)

if nargin<7
    % alpha fn weights
    dt = mean(diff(tin));
    tau = 300;
    twin = linspace(0,tau*10,tau*10/dt+1);
    wfun = twin/tau.*exp(1-twin/tau);
    wfun = wfun/sum(wfun);
end

fitw.pattern_ids = pattern_ids;
fitw.fwd_pattern_mask = conv(pattern_ids,wfun);
fitw.fwd_pattern_mask = fitw.fwd_pattern_mask(1:length(tin));
fwd_spk_weights = interp1(tin,fitw.fwd_pattern_mask,Tlist{ij(1)},'nearest')';
fwd_spk_weights = fwd_spk_weights/mean(fwd_spk_weights); % total weight is same as unweighted result
fitw.fwd = composite_model_fit_weighted(ij,Tlist,params,'i>j',fwd_spk_weights,syn_init);

fitw.rev_pattern_mask = conv(fliplr(pattern_ids),wfun);
fitw.rev_pattern_mask = fliplr(fitw.rev_pattern_mask(1:length(tin)));
rev_spk_weights = interp1(tin,fitw.rev_pattern_mask,Tlist{ij(1)},'nearest')';
rev_spk_weights = rev_spk_weights/mean(rev_spk_weights);
fitw.rev = composite_model_fit_weighted(ij,Tlist,params,'i>j',rev_spk_weights,syn_init);
