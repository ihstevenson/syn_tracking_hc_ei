


%% load table with all detected pairs and metadata

load('results/syn_index.mat')
load('results/wave_all.mat')

% counts
height(syn_table) % total pairs
sum(strcmp(syn_table.pre_anat,syn_table.post_anat)) % intra-area pairs
sum(strcmp(syn_table.pre_anat,syn_table.post_anat) & strcmp(syn_table.pre_anat,'CA1')) % intra-CA1

% putative synapse statistics
median(wave_all(wave_all(:,4)>0,2)) % pre waveform durations
median(wave_all(wave_all(:,5)>0,2)) % post waveform durations

prctile(wave_all(wave_all(:,4)>0,3),[25 50 75]) % pre rates
prctile(wave_all(wave_all(:,5)>0,3),[25 50 75]) % post rates

prctile(syn_table.syn_latency,[25 50 75])
prctile(syn_table.syn_time_const,[25 50 75])
prctile(syn_table.syn_efficacy,[25 50 75])
prctile(syn_table.syn_contrib,[25 50 75])


%% fig 2a anatomy for putative pairs
area_list = {'CA1','CA3','DG','SUB','ProS'};
area_id = [382,463,726,502,484682470];
area_count = [];
syn_count=[];
for i=1:length(area_list)
    area_count(i)=sum(wave_all(:,6)==area_id(i));
    for j=1:length(area_list)
        syn_count(i,j)=sum(strcmp(syn_table.pre_anat,area_list{i}) & strcmp(syn_table.post_anat,area_list{j}));
    end
end
syn_count

urec = unique(wave_all(:,7));
pair_count=zeros(length(area_list));
total_pairs=0;
tmp=zeros(length(area_list),1);
for r=1:length(urec)
    total_pairs = total_pairs+sum(wave_all(:,7)==urec(r))*(sum(wave_all(:,7)==urec(r))-1);
    for i=1:length(area_list)
        tmp(i)=sum(wave_all(:,6)==area_id(i) & wave_all(:,7)==urec(r));
    end
    pair_count=pair_count+tmp*tmp'-diag(tmp);
end

figure(21)
subplot(5,1,1)
imagesc(log(area_count))
text(1:5,ones(1,5)-.4,area_list,'HorizontalAlignment','center')
text(1:5,ones(1,5),arrayfun(@num2str,area_count,'UniformOutput',false),'HorizontalAlignment','center')
box off; axis off

subplot(5,1,2:5)
imagesc(log(pair_count))
for i=1:length(area_list)
    text(1:5,ones(1,5)+i-1,arrayfun(@num2str,syn_count(i,:),'UniformOutput',false),'HorizontalAlignment','center')
end
axis equal
box off; axis off
colormap(cbrew(256,'blues'))

%% fig 2b waveform features

figure(22)
clf
[n1,x1e,y1e] = histcounts2(log10(wave_all(:,1)),log10(wave_all(:,2)),256);
[n2,x2e,y2e] = histcounts2(log10(wave_all(:,3)),log10(wave_all(:,2)),256);
x1 = x1e(1:end-1)+diff(x1e)/2;
y1 = y1e(1:end-1)+diff(y1e)/2;
x2 = x2e(1:end-1)+diff(x2e)/2;
y2 = y2e(1:end-1)+diff(y2e)/2;
gwin = gausswin(64);
gwin = gwin*gwin';
s1 = conv2(n1,gwin,'same');
s2 = conv2(n2,gwin,'same');

subplot(3,5,[6 7 11 12])
imagesc(x1,y1,log(s1'+50))
set(gca,'YDir','normal')
set(gca,'TickDir','out')
y1r = unique(round((10.^y1)*10)/10);
set(gca,'YTick',log10(y1r))
set(gca,'YTickLabels',y1r)
x1r = unique(round((10.^x1)*10)/10);
set(gca,'XTick',log10(x1r))
set(gca,'XTickLabels',x1r)
box off
xlabel('Waveform Half-width [ms]')
ylabel('Waveform Duration [ms]')
xlim([min(x1) log10(0.6)])
ylim([min(y1) log10(1.4)])


cmap = lines(2);
% avoid overplottting
jit_wave = wave_all(:,1:2);
jit_wave(:,1) = jit_wave(:,1)+(rand(size(jit_wave(:,1)))-.5)*0.0137;
jit_wave(:,2) = jit_wave(:,2)+(rand(size(jit_wave(:,2)))-.5)*0.0137;

hold on
scatter(log10(jit_wave(wave_all(:,4)>0,1)),log10(jit_wave(wave_all(:,4)>0,2)),10,'MarkerFaceColor',cmap(1,:),'MarkerEdgeColor','none')
scatter(log10(jit_wave(wave_all(:,5)>0,1)),log10(jit_wave(wave_all(:,5)>0,2)),10,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','none')
hold off

subplot(3,5,[8 9 13 14])
imagesc(x2,y2,log(s2'+50))
set(gca,'YDir','normal')
set(gca,'TickDir','out')
set(gca,'YDir','normal')
set(gca,'TickDir','out')
y2r = unique(round((10.^y2)*10)/10);
set(gca,'YTick',log10(y2r))
set(gca,'YTickLabels',y2r)
x2r = 2.^[-3:6];
set(gca,'XTick',log10(x2r))
set(gca,'XTickLabels',x2r)
box off
colormap(1-gray(256))
xlabel('Firing Rate [Hz]')
ylabel('Waveform Duration [ms]')
ylim([min(y1) log10(1.4)])
xl=xlim();

hold on
scatter(log10(wave_all(wave_all(:,4)>0,3)),log10(jit_wave(wave_all(:,4)>0,2)),10,'MarkerFaceColor',cmap(1,:),'MarkerEdgeColor','none')
scatter(log10(wave_all(wave_all(:,5)>0,3)),log10(jit_wave(wave_all(:,5)>0,2)),10,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','none')
hold off

subplot(3,5,1:2)
edges = linspace(log10(0.06),log10(0.65),256);
[fi,x]=ksdensity(log10(jit_wave(wave_all(:,4)>0,1)),edges,'bandwidth',range(edges)/15);
[fe,x]=ksdensity(log10(jit_wave(wave_all(:,5)>0,1)),edges,'bandwidth',range(edges)/15);
shadedErrorBar(x,fi*0,[fi; fi*0]);
hold on
shadedErrorBar(x,fe*0,[fe; fe*0]);
hold off
box off; set(gca,'TickDir','out')
xlim([min(x1) log10(0.6)])

subplot(3,5,3:4)
edges = linspace(log10(0.1),log10(100),256);
[fi,x]=ksdensity(log10(wave_all(wave_all(:,4)>0,3)),edges,'bandwidth',range(edges)/15);
[fe,x]=ksdensity(log10(wave_all(wave_all(:,5)>0,3)),edges,'bandwidth',range(edges)/15);
shadedErrorBar(x,fi*0,[fi; fi*0]);
hold on
shadedErrorBar(x,fe*0,[fe; fe*0]);
hold off
box off; set(gca,'TickDir','out')
xlim(xl)

subplot(3,5,[10 15])
edges = linspace(log10(0.08),log10(1.5),256);
[fi,x]=ksdensity(log10(wave_all(wave_all(:,4)>0,2)),edges,'bandwidth',range(edges)/15);
[fe,x]=ksdensity(log10(wave_all(wave_all(:,5)>0,2)),edges,'bandwidth',range(edges)/15);
shadedErrorBar(x,fi*0,[fi; fi*0]);
hold on
shadedErrorBar(x,fe*0,[fe; fe*0]);
hold off
box off; set(gca,'TickDir','out')
camroll(90)
set(gca,'XDir','reverse')
xlim([min(y1) log10(1.4)])

%% fig 2c - putative synapse statistics

figure(23)
subplot(3,1,1)
histogram(syn_table.syn_latency,25,'EdgeColor','none')
xlabel('Latency [ms]')
box off; set(gca,'TickDir','out')
subplot(3,1,2)
histogram(syn_table.syn_time_const,25,'EdgeColor','none')
xlabel('Time Constant [ms]')
box off; set(gca,'TickDir','out')
subplot(3,1,3)
tmp = syn_table.syn_efficacy*100;
tmp(tmp>23)=23;
histogram(tmp,25,'EdgeColor','none')
xlabel('Efficacy [%]')
box off; set(gca,'TickDir','out')
xlim([0 23])
