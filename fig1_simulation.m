
%% Simulate BCM rule with naturalistic fluctuating pre- and post-rates

clear all; clc; clf;

c1 = 3.5; c2 = -10; % coefficients for hebbian example BCM
% c1 = -2.5; c2 = 8; % coefficients for anti-hebbian example

% simulating average firing rate over a 2 hour recording
rl = 7200; % recording length [s]
tstep = 1; % timestep [s]

rng(10);
max_rate = 4.25;
tau_BCM  = 3e3;
tau_fast = 90;

samp1 = rand(200,1)*3;
firing_pre = exp(interp1(linspace(0,rl,length(samp1)),samp1,linspace(1,rl/tstep,rl/tstep),'spline'));
samp2 = rand(200,1)*3;
firing_post = exp(interp1(linspace(0,rl,length(samp2)),samp2,linspace(1,rl/tstep,rl/tstep),'spline'));

firing_pre = firing_pre/max(firing_pre)*max_rate;
firing_post = firing_post/max(firing_post)*max_rate;

dw(:,1) = firing_post.^2.*firing_pre;
dw(:,2) = firing_post.*firing_pre;
dw = dw-mean(dw);
w1 = zeros(size(dw));
w2 = zeros(size(dw));

for i = 2:length(w1)
    w1(i,:) = w1(i-1,:) - 1/tau_BCM*w1(i-1,:) + dw(i,:);
    w2(i,:) = w2(i-1,:) - 1/tau_fast*w2(i-1,:) + dw(i,:);
end

w = w1 - w2;
w_val = c1*w(:,1) + c2*w(:,2);

% Plot rate and weight timeseries
plot_tick = 5; %to make the plot separations between the pre and post firing rate series
figure(1)
ax1 = subplot(4,4,1:3);
plot(firing_pre, 'blue',Linewidth=2);
hold on;
plot(firing_post+plot_tick, 'red',LineWidth=2);
hold off;
box off; set(gca,'TickDir','out')
ax=gca;
ax2 = subplot(4,4,5:7);
plot(w_val, 'k', LineWidth=2);
hold off;
linkaxes([ax1,ax2],'x')
axis tight
box off; set(gca,'TickDir','out')

% Getting BCM graph
subplot(4,4,4)
pts = linspace(0,max(firing_post),256);
plot(pts,(c1*pts.^2+c2*pts),'k-','LineWidth',2);
yline(0,'--')
box off; set(gca,'TickDir','out')

subplot(4,4,8)
[pos,pre]=meshgrid(pts,pts);
deltaw = c1*pos.^2.*pre+c2*pos.*pre;
imagesc(pts,pts,deltaw)
set(gca,'YDir','normal')
set(gca,'CLim',[-1 1].*max(abs(deltaw(:))))
box off; set(gca,'TickDir','out')


%% Simulate BCM rule with controlled stimulation

firing_pre = zeros(1,rl/tstep);
firing_post = zeros(1,rl/tstep);
ind_fr = 1:rl/20:rl;
firing_pre(ind_fr) = ones(1,length(ind_fr))';%fixed pre
firing_post(ind_fr) = rand(length(ind_fr),1)'*max_rate;

clear weight;
dw(:,1) = firing_post.^2.*firing_pre;
dw(:,2) = firing_post.*firing_pre;
dw = dw-mean(dw);
w1 = zeros(size(dw));
w2 = zeros(size(dw));
for i = 2:length(w1)
    w1(i,:) = w1(i-1,:) - 1/tau_BCM*w1(i-1,:) + dw(i,:);
    w2(i,:) = w2(i-1,:) - 1/tau_fast*w2(i-1,:) + dw(i,:);
end

w = w1 - w2;
w_val = c1*w(:,1) + c2*w(:,2);

% plot rate and weight timeseries

plot_tick = 4;%to make the plot separations between the pre and post firing rate series
figure(1)
ax1 = subplot(4,4,9:11);
plot(firing_pre, 'blue',Linewidth=2);
hold on;
plot(firing_post+plot_tick, 'red',LineWidth=2);
hold off;
box off; set(gca,'TickDir','out')
ax=gca;
ax2 = subplot(4,4,13:15);
plot(w_val, 'k', LineWidth=2);
hold off;
linkaxes([ax1,ax2],'x')
axis tight
box off; set(gca,'TickDir','out')