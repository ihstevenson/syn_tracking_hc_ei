
%% Run synapse detection on all recordings

dat_dir = 'D:\allen_neuropixel\';
output_dir = 'Z:\data\fluctuation_results\';
flist = dir([dat_dir '*_units.csv']);

for rec_num = [1:length(flist)]
    prefix = strsplit(flist(rec_num).name,'_');
    prefix = prefix{1};

    [data,detection] = syn_detection(dat_dir,prefix);
    save([output_dir prefix '_detection'],'data','detection')

    % optimize window size for detected pairs
    for cnx = 1:size(detection.ijlist,1)
        fsyn = detection.fitout(detection.fidx(cnx));
        detection.params.nknots=2;
        detection.params.doPlot = false;
        [opt_res(cnx),smooth_estimate(cnx)] = optimizeWindowSize(detection.ijlist(cnx,:),data.Tlist,detection.params,fsyn,32,0);
    end
    save([output_dir prefix '_opt_window'],'data','detection')
end

%% Generate index of detected pairs

flist = dir([output_dir '*_detection.mat']);
syn_table = table({'0'},0,0,{'0'},{'0'},0,0,0,0,0,0,'VariableNames',{'ecephys_session','pre_unit_id','post_unit_id','pre_anat','post_anat',...
    'pre_waveform_halfwidth','post_waveform_halfwidth','pre_waveform_duration','post_waveform_duration','pre_firing_rate','post_firing_rate'});
syn_table=table();
wave_all=[];
for rec_num = [1:length(flist)]
    prefix = strsplit(flist(rec_num).name,'_');
    prefix = prefix{1};
    load([output_dir prefix '_detection'])
    %load([output_dir prefix '_opt_window'])

    row_ids = [(height(syn_table)+1):(height(syn_table)+size(detection.ijlist,1))];
    if ~isempty(detection.ijlist)
        syn_table.ecephys_session(row_ids) = {prefix};
        syn_table.pre_unit_id(row_ids) = data.tab_sub.unit_id(detection.ijlist(:,1));
        syn_table.post_unit_id(row_ids) = data.tab_sub.unit_id(detection.ijlist(:,2));
        syn_table.pre_anat(row_ids) = data.tab_sub.ecephys_structure_acronym(detection.ijlist(:,1));
        syn_table.post_anat(row_ids) = data.tab_sub.ecephys_structure_acronym(detection.ijlist(:,2));

        syn_table.pre_waveform_halfwidth(row_ids) = data.tab_sub.waveform_halfwidth(detection.ijlist(:,1));
        syn_table.post_waveform_halfwidth(row_ids) = data.tab_sub.waveform_halfwidth(detection.ijlist(:,2));
        syn_table.pre_waveform_duration(row_ids) = data.tab_sub.waveform_duration(detection.ijlist(:,1));
        syn_table.post_waveform_duration(row_ids) = data.tab_sub.waveform_duration(detection.ijlist(:,2));
        syn_table.pre_firing_rate(row_ids) = data.tab_sub.firing_rate(detection.ijlist(:,1));
        syn_table.post_firing_rate(row_ids) = data.tab_sub.firing_rate(detection.ijlist(:,2));

        syn_table.pre_probeid(row_ids) = data.tab_sub.probe_id(detection.ijlist(:,1));
        syn_table.post_probeid(row_ids) = data.tab_sub.probe_id(detection.ijlist(:,2));

        syn_table.syn_latency(row_ids) = detection.fit_res.lat(detection.fit_res.kidx);
        syn_table.syn_time_const(row_ids) = detection.fit_res.dur(detection.fit_res.kidx);
        syn_table.syn_efficacy(row_ids) = detection.fit_res.eff(detection.fit_res.kidx);
        syn_table.syn_contrib(row_ids) = cellfun(@length,data.Tlist(detection.ijlist(:,1))).*detection.fit_res.eff(detection.fit_res.kidx)'./cellfun(@length,data.Tlist(detection.ijlist(:,2)));

        syn_table.rec_num(row_ids) = rec_num;
        syn_table.local_id(row_ids) = 1:length(row_ids);
    end
    w=[];
    w(:,1) = data.tab_sub.waveform_halfwidth;
    w(:,2) = data.tab_sub.waveform_duration;
    w(:,3) = data.tab_sub.firing_rate;
    w(:,4:5) = 0;
    if ~isempty(detection.ijlist)
        w(unique(detection.ijlist(:,1)),4)=1;
        w(unique(detection.ijlist(:,2)),5)=1;
    end
    w(:,6) = data.tab_sub.ecephys_structure_id;
    w(:,7) = str2num(prefix);
    wave_all=[wave_all; w];
end


% save('results/syn_index','syn_table')
save('results/wave_all','wave_all')

%% Generate paired spike times for all putative synapses

flist = dir([output_dir '*_detection.mat']);
cc=1;
Tlist_pre = cell(0);
Tlist_post = cell(0);
Tlist_inv = [];
T=[];
for rec_num = [1:length(flist)]
    prefix = strsplit(flist(rec_num).name,'_');
    prefix = prefix{1};
    load([output_dir prefix '_detection'])
    for cnx=1:size(detection.ijlist,1)
        Tlist_pre{cc}=data.Tlist{detection.ijlist(cnx,1)};
        Tlist_post{cc}=data.Tlist{detection.ijlist(cnx,2)};
        T(cc) = max(cellfun(@max,data.Tlist));
        Tlist_inv
        cc=cc+1;
    end
end

save([output_dir 'prepost_spike_times'],'Tlist_pre','Tlist_post','T')

%% Generate subset of pairs for demos and figures

% load([output_dir 'prepost_spike_times'])

% Fig 2D (rec_num-ij_idx)
% 1-13, 35-3
% 51-8, 42-10
% 29-26, 25-25

% Fig 4A [14-14]
% Fig 4B [25-9]
% [35-5, 25-13]
% Fig 4C [1-6,29-3,25-28]
% Fig 4D [51-10,25-1,40-10]
% Fig 4E [2-27,32-9,34-9]

% Fig 5A 100 % 13-7
% Fig 5B 174 % 19-14


subset_index = [13; % 1-13
    443; % 35-3
    658; % 51-8
    550; % 42-10
    310; % 29-26
    239; % 25-25

    124; % 14-14
    223; % 25-9
    6; % 1-6
    287; % 29-3
    242; % 25-28
    660; % 51-10
    215; % 25-1
    514; % 40-10
    42; % 2-27
    386; % 32-9
    433; % 34-9

    100; % 13-7
    174; % 19-14
    ];
% syn_table{full_index,2:3}
syn_table{subset_index,[end-1:end]}

Tlist_pre=Tlist_pre(subset_index);
Tlist_post=Tlist_post(subset_index);
T=T(subset_index);
save(['results/subset_prepost_spike_times'],'Tlist_pre','Tlist_post','T','subset_index')