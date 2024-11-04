% function [varargout] = chancompare(dat, idx_sub, good_chan, bad_chan)

% function chancompare is used to generate figure for a comparison of good 
% channels and bad channels

% Parameters:
%       data: <nirs.core.Data>, dOD data array
%       idx_sub: <int>, index of the subject in dOD data array
%       good_chan: <int>, good channel index, for 760nm
%       bad_chan:  <int>, bad cahnnel index, for 760nm
%       time_idx:  a certain interval wanted to be plotted
%       varargin:  other possible inputs
%           TODO: varargin

% Output:
%       

% XW last edited, 15.10.2024
%% Variables define and parameters settings

% % TODO: whether automating?
% % Using these lines, determine which subject is used to plot
% job = nirs.modules.QT;
% job.qThreshold = 0.65;
% job.sciThreshold = 0.6;
% job.pspThreshold = 0.1;
% SQ_all = job.run(dOD);
% 
% all_bad_windows = arrayfun(@(x) x.qMats.bad_windows, SQ_selected, 'UniformOutput', false);
% all_bad_channels = arrayfun(@(x) x.qMats.bad_links, SQ_selected, 'UniformOutput', false);
% 
% Nr_sub = 1; % with largest amount of bad channels
% 
% bad_win = all_bad_windows{Nr_sub};
% bad_win_selected = 48;



%%
idx_sub = 11;
bad_win_selected = 48;
bad_chan_idx = 1;
good_chan_idx = 2;

dOD_select = dOD(idx_sub);

job = nirs.modules.QT;
job.qThreshold = 0.65;
job.sciThreshold = 0.6;
job.pspThreshold = 0.1;
SQ_selected = job.run(dOD_select);

% number of windows calculation in 'qtnirs.m' line 692 693
% nirsplot_param.window in second, see in 'QT.m'
% window_samples_ = nirsplot_param.window*nirsplot_param.fs;
% n_windows_ = (length(raw.t)-overlap_samples_)/(window_samples_-overlap_samples_);
win_length = SQ_selected.qMats.sampPerWindow;
% win_length = SQ_selected(idx_sub).qMats.sampPerWindow;
time_idx = (bad_win_selected-20) * win_length : (bad_win_selected + 20) * win_length;
dOD_select = dOD(idx_sub);



bad_chan_list = SQ_selected.qMats.bad_links;
bad_chan = bad_chan_list(bad_chan_idx);
good_chan_list = find(SQ_selected.qMats.good_combo_link(:,3) == ...
    max(SQ_selected.qMats.good_combo_link(:,3)));
good_chan = good_chan_list(good_chan_idx);

good_chan = 10;

zoom_in_strat = bad_win_selected * win_length;
zoom_in = zoom_in_strat: zoom_in_strat + ceil(20 * dOD_select.Fs); % for 20 second

%% Plot A
% 1. Plot: bad channel plotting in time domain
figure;
fig_size = [0, 0, .70, 0.80];
set(gcf,'unit', 'normalized','Position', fig_size)

subplot(4,4,[1,2,5,6])
dat = dOD_select.data;
time_vec = dOD_select.time;
% draw plots
plot(time_vec(time_idx), dat(time_idx, bad_chan), 'b');
hold on
plot(time_vec(time_idx), dat(time_idx, bad_chan + 24), 'r');
yline(0, '--', 'Color', [0.5 0.5 0.5],'LineWidth',1); 

% plot settings
y1 = ylabel('ΔOD');
% xlabel('Time (s)');
ylim([-0.08,0.08])
title('Poor Channel Quality')

% Move ylabel towards plot
set(y1, 'Units', 'normalized');  
pos = get(y1, 'Position');       
pos(1) = pos(1) + 0.05;          
set(y1, 'Position', pos); 

ticks_x = [time_vec(time_idx(1)), time_vec(time_idx(end))];
ticks_y = [-0.06, 0, 0.06];
set(gca, 'XTick', ticks_x, ...
         'XTickLabel', {sprintf('%.1f', ticks_x(1)), sprintf('%.1f', ticks_x(2))}, ...
         'YTick', ticks_y, ...
         'YTickLabel', {'-0.06', '0', '0.06'}, ...
         'FontSize', 12);

ax = gca;
ax.XLim = [time_vec(time_idx(1)), time_vec(time_idx(end))];
ax.XAxisLocation = 'bottom';
ax.YAxis.TickLength = [0.02, 0.02];

% serial number
annotation('textbox', [0.10, 0.93, 0.05, 0.05], 'String', 'A', 'FontSize', 14, ...
    'FontWeight', 'bold', 'EdgeColor', 'none');

% index add into plot
SCI = mean(SQ_selected.qMats.sci_array(bad_chan,:));
PSP = mean(SQ_selected.qMats.power_array(bad_chan,:));

text(time_vec(time_idx(1)) + 5, 0.07, ['SCI = ', sprintf('%.2f', SCI)], 'FontSize', 9);
text(time_vec(time_idx(1)) + 5, 0.060,['PSP = ', sprintf('%.2f', PSP)], 'FontSize', 9);

% draw a zoom-in area

x = time_vec(zoom_in_strat);
y = -0.03;
w = 20; % in second
h = 0.08;

rectangle('Position', [x, y, w, h], ...
    'EdgeColor', 'k', 'LineWidth', 1);

annotation('line', [0.35, 0.54], [0.80, 0.78], 'Color', 'k');

%% Plot B
subplot(4,4,[3,7])

plot(time_vec(zoom_in), dat(zoom_in, bad_chan), 'b');
hold on
plot(time_vec(zoom_in), dat(zoom_in, bad_chan + 24), 'r');
yline(0, '--', 'Color', [0.5 0.5 0.5],'LineWidth',1); 

y2 = ylabel('ΔOD');
% xlabel('Time (s)');
ylim([-0.08,0.08])

% Move ylabel towards plot
set(y2, 'Units', 'normalized');  
pos = get(y2, 'Position');       
pos(1) = pos(1) + 0.05;          
set(y2, 'Position', pos); 

annotation('textbox', [0.50, 0.93, 0.05, 0.05], 'String', 'B', 'FontSize', 14, ...
    'FontWeight', 'bold', 'EdgeColor', 'none');

% Axis setting
ticks_x = [time_vec(zoom_in(1)), time_vec(zoom_in(end))];
ticks_y = [-0.06, 0, 0.06];
set(gca, 'XTick', ticks_x, ...
         'XTickLabel', {sprintf('%.1f', ticks_x(1)), sprintf('%.1f', ticks_x(2))}, ...
         'YTick', ticks_y, ...
         'YTickLabel', {'-0.06', '0', '0.06'}, ...
         'FontSize', 12);
ax = gca;
ax.XLim = [time_vec(zoom_in(1)), time_vec(zoom_in(end))];
ax.XAxisLocation = 'bottom';
ax.YAxis.TickLength = [0.02, 0.02];

%% Plot C
subplot(4,4,4)

Fs = dOD_select.Fs;
L = length(dat(zoom_in, bad_chan));
f  = Fs*(0:(floor(L/2)))/L;

dat_fft = fft(dat(zoom_in, bad_chan));
dat_fft = abs(dat_fft / L);
dat_fft = dat_fft(1 : floor(L/2) + 1);
dat_fft(2:end-1) = 2*dat_fft(2:end-1);

plot(f, dat_fft,'b');
xlim([0, 2])
ylabel('|P(data)|')

subplot(4,4,8)
dat_fft = fft(dat(zoom_in, bad_chan + 24));
dat_fft = abs(dat_fft / L);
dat_fft = dat_fft(1 : floor(L/2) + 1);
dat_fft(2:end-1) = 2*dat_fft(2:end-1);

plot(f, dat_fft, 'r');
xlim([0, 2])
ylabel('|P(data)|')


% [p1, f1] = pspectrum(dat(zoom_in, bad_chan));
% ylim([-100,-40])
% plot(f1,pow2db(p1),'b');
% ylabel('Power(dB)');
% % xlabel('Frequency(Hz)');

% subplot(4,4,8)
% [p2, f2] = pspectrum(dat(zoom_in, bad_chan + 24));
% ylim([-100,-40])
% plot(f2,pow2db(p2),'r');
% ylabel('Power(dB)');
% % xlabel('Frequency(Hz)');

annotation('textbox', [0.72, 0.93, 0.05, 0.05], 'String', 'C', 'FontSize', 14, ...
    'FontWeight', 'bold', 'EdgeColor', 'none');

%% Plot D
subplot(4,4,[9,10,13,14])

% draw plots
plot(time_vec(time_idx), dat(time_idx, good_chan), 'b');
hold on
plot(time_vec(time_idx), dat(time_idx, good_chan + 24), 'r');
yline(0, '--', 'Color', [0.5 0.5 0.5],'LineWidth',1); 

% plot settings
y3 = ylabel('ΔOD');
xlabel('Time (s)');
ylim([-0.08,0.08])
title('Good Channel Quality')

% Move ylabel towards plot
set(y3, 'Units', 'normalized');  
pos = get(y3, 'Position');       
pos(1) = pos(1) + 0.05;          
set(y3, 'Position', pos); 

ticks_x = [time_vec(time_idx(1)), time_vec(time_idx(end))];
ticks_y = [-0.06, 0, 0.06];
set(gca, 'XTick', ticks_x, ...
         'XTickLabel', {sprintf('%.1f', ticks_x(1)), sprintf('%.1f', ticks_x(2))}, ...
         'YTick', ticks_y, ...
         'YTickLabel', {'-0.06', '0', '0.06'}, ...
         'FontSize', 12);

ax = gca;
ax.XLim = [time_vec(time_idx(1)), time_vec(time_idx(end))];
ax.XAxisLocation = 'bottom';
ax.YAxis.TickLength = [0.02, 0.02];

% serial number
annotation('textbox', [0.10, 0.47, 0.05, 0.05], 'String', 'D', 'FontSize', 14, ...
    'FontWeight', 'bold', 'EdgeColor', 'none');

% index add into plot
SCI = mean(SQ_selected.qMats.sci_array(good_chan,:));
PSP = mean(SQ_selected.qMats.power_array(good_chan,:));

text(time_vec(time_idx(1)) + 5, 0.07, ['SCI = ', sprintf('%.2f', SCI)], 'FontSize', 9);
text(time_vec(time_idx(1)) + 5, 0.060,['PSP = ', sprintf('%.2f', PSP)], 'FontSize', 9);

% draw a zoom-in area

x = time_vec(zoom_in_strat);
y = -0.03;
w = 20; % in second
h = 0.08;

rectangle('Position', [x, y, w, h], ...
    'EdgeColor', 'k', 'LineWidth', 1);

annotation('line', [0.35, 0.54], [0.39, 0.35], 'Color', 'k');



%%
subplot(4,4,[11,15])

plot(time_vec(zoom_in), dat(zoom_in, good_chan), 'b');
hold on
plot(time_vec(zoom_in), dat(zoom_in, good_chan + 24), 'r');
yline(0, '--', 'Color', [0.5 0.5 0.5],'LineWidth',1); 

y4 = ylabel('ΔOD');
xlabel('Time (s)');
ylim([-0.08,0.08])

% Move ylabel towards plot
set(y4, 'Units', 'normalized');  
pos = get(y4, 'Position');       
pos(1) = pos(1) + 0.05;          
set(y4, 'Position', pos); 

annotation('textbox', [0.50, 0.47, 0.05, 0.05], 'String', 'E', 'FontSize', 14, ...
    'FontWeight', 'bold', 'EdgeColor', 'none');

% Axis setting
ticks_x = [time_vec(zoom_in(1)), time_vec(zoom_in(end))];
ticks_y = [-0.06, 0, 0.06];
set(gca, 'XTick', ticks_x, ...
         'XTickLabel', {sprintf('%.1f', ticks_x(1)), sprintf('%.1f', ticks_x(2))}, ...
         'YTick', ticks_y, ...
         'YTickLabel', {'-0.06', '0', '0.06'}, ...
         'FontSize', 12);
ax = gca;
ax.XLim = [time_vec(zoom_in(1)), time_vec(zoom_in(end))];
ax.XAxisLocation = 'bottom';
ax.YAxis.TickLength = [0.02, 0.02];

%%
% subplot(4,4,12)
% [p1, f1] = pspectrum(dat(zoom_in, good_chan));
% line1 = plot(f1,pow2db(p1),'b');
% ylim([-100,-40])
% ylabel('Power(dB)');
% % xlabel('Frequency(Hz)');
% 
% subplot(4,4,16)
% [p2, f2] = pspectrum(dat(zoom_in, good_chan + 24));
% line2 = plot(f2,pow2db(p2),'r');
% ylim([-100,-40])
% ylabel('Power(dB)');
% xlabel('Frequency(Hz)');
subplot(4,4,12)
dat_fft = fft(dat(zoom_in, good_chan));
dat_fft = abs(dat_fft / L);
dat_fft = dat_fft(1 : floor(L/2) + 1);
dat_fft(2:end-1) = 2*dat_fft(2:end-1);

line1 = plot(f, dat_fft,'b');
xlim([0, 2])
ylabel('|P(data)|')

subplot(4,4,16)
dat_fft = fft(dat(zoom_in, good_chan + 24));
dat_fft = abs(dat_fft / L);
dat_fft = dat_fft(1 : floor(L/2) + 1);
dat_fft(2:end-1) = 2*dat_fft(2:end-1);

line2 = plot(f, dat_fft, 'r');
xlim([0, 2])
ylabel('|P(data)|')
xlabel('Frequency(Hz)');

annotation('textbox', [0.72, 0.47, 0.05, 0.05], 'String', 'F', 'FontSize', 14, ...
    'FontWeight', 'bold', 'EdgeColor', 'none');
%%
fig_handles = [line1, line2];
lgd = legend(fig_handles, {'HbR (760nm)', 'HbO (850nm)'}, ...
            'Location', 'southoutside', 'Orientation', 'horizontal');
lgd.Position = [0.40, 0.005, 0.3, 0.05];
lgd.Box = 'off';