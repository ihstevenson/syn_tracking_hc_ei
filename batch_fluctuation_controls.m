
%% Generate optimal window size estimate for control/shuffled data

glist = dir([output_dir '*_detection.mat']);
for rec_num=1:length(glist)
    prefix = strsplit(glist(rec_num).name,'_');
    prefix = prefix{1};
    load([output_dir prefix '_detection'])
    load([output_dir prefix '_opt_window'])

    clear opt_res_shuff smooth_estimate_shuff
    if ~isempty(detection.ijlist)
        for cnx = 1:size(detection.ijlist,1)
            ij = detection.ijlist(cnx,:);
            fsyn = fitout(detection.fidx(cnx));
            sparams=params;
            syn_v = fsyn.atoms.syn1(exp(fsyn.feat_params))>0.01;
            syn_t = params.t(syn_v);
            sparams.interval = [min(syn_t) max(syn_t)];

            % Shuffle randomly (negative control)
            % Tshuff = shuffle_preserve_efficacy(Tlist,ij,sparams);

            % Shuffle weighted (positive control)
            %get baseline fraction
            pfrac = 1-sum(fsyn.lamsmoo(syn_v))/sum(fsyn.yhat(syn_v));
            Tshuff = shuffle_manip_efficacy(Tlist,ij,sparams,pfrac,1/(60*40));

            params.nknots=6;
            params.doPlot = false;

            figure(2)
            [opt_res_shuff(cnx),smooth_estimate_shuff(cnx)] = optimizeWindowSize(ijlist(cnx,:),Tshuff,params,fsyn,20,1);
        end
        save([output_dir_pos prefix '_opt_window_shuff.mat'],'opt_res_shuff','smooth_estimate_shuff')
    end
end


%% Collate optimal window results for all pairs...

fluc.optimal_window_size=[];
fluc.cv=[];

output_dir_neg = 'E:\fluctuation_control_neg\';

cc=1;
glist = dir([output_dir '*_detection.mat']);
for rec_num=1:length(glist)
    prefix = strsplit(glist(rec_num).name,'_');
    prefix = prefix{1};
    load([output_dir prefix '_detection'])

    if ~isempty(detection.ijlist)
        load([output_dir prefix '_opt_window'])
        load([output_dir_neg prefix '_opt_window_shuff'])
        for cnx = 1:size(detection.ijlist,1)
            fsyn = detection.fitout(detection.fidx(cnx));
            fluc.optimal_window_size(cc,1) = smooth_estimate(cnx).window;
            fluc.cv(cc,1) = nanstd(smooth_estimate(cnx).eff)/fsyn.eff;

            fluc.optimal_window_size(cc,2) = smooth_estimate_shuff(cnx).window;
            fluc.cv(cc,2) = nanstd(smooth_estimate_shuff(cnx).eff)/fsyn.eff;
            cc=cc+1;
        end
    end
end
save('results\optimal_fluctuation_windows','fluc')

%%

