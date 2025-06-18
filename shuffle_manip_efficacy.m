
% returns post-synaptic spike times for a null-model
% where efficacy does not fluctuate over the recording

function Tshuff = shuffle_manip_efficacy(Tlist,ij,params,shift_frac,w)

[c,dT] = spk_xcorr(Tlist{ij(1)},Tlist{ij(2)},params.interval(1),params.interval(2),params.num_bins);

keepidx=zeros(size(dT,1),1);
ua = unique(dT(:,3));
for i=1:length(ua)
    r = find(dT(:,3)==ua(i));
    keepidx(datasample(r,1))=1;
end
n = length(ua);
ShuffPost = Tlist{ij(2)};
ShuffPost = Tlist{ij(2)}(setdiff(1:length(ShuffPost),unique(dT(:,3)))); % keep spikes that aren't in the xcorr
wt = (sin(2*pi*w*Tlist{ij(1)})+1)/2*shift_frac+(1-shift_frac);
presample = datasample(Tlist{ij(1)},n,'Weights',wt);
intsample = datasample(dT(keepidx>0,1),n);
ShuffPost = sort([ShuffPost; presample+intsample]);
Tshuff=Tlist;
Tshuff{ij(2)}=ShuffPost;