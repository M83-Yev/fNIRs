%function [varargout] = IIR_plot(dHbX, Nr_Sub, Nr_chan, IIR_order, FIR_order, lcut, hcut, type)

% IIR_PLOT Plots filtered and unfiltered fNIRS data.
%
% Usage:
%   IIR_plot(data, Nr_Sub, Nr_chan, order, lcut, hcut)
%
% Inputs:
%   dHbX     - <nirs.core.Data> object, dHbX data
%   Nr_Sub   - (integer) Subject number
%   Nr_chan  - (integer) Channel number
%   IIR_order- (integer) Filter order
%   lcut     - (double) Low-pass filter cutoff frequency in Hz
%   hcut     - (double) High-pass filter cutoff frequency in Hz
%   type     - (char) 'IIR', 'IIR_FIR'
%
% The function generates four plots:
%   - Unfiltered data
%   - Low-pass filtered data
%   - High-pass filtered data
%   - Band-pass filtered data
%
% Each plot shows the change in hemoglobin concentration for HbR and HbO.

% XW last edited, 11.10.2024, MATLAB R2024b

%% Testing
% values to define
Nr_Sub = 6;
Nr_chan = 1;
data = dHbX(Nr_Sub);
IIR_order = 2;
FIR_order = 1000;
lcut = 0.09; % Hz
hcut = 0.01; % Hz
type = 'IIR_FIR';   % 'IIR_FIR': IIR vs. FIR, 'IIR': all filters

% values can be calculated or extracted
fs = data.Fs;   % sampling rate
nqfreq = fs/2;  % nyquist frequency
onsets = cellfun(@(x) x.onset, data.stimulus.values, 'UniformOutput', false);   % extract data stimuli onsets
stim_dur = cellfun(@(x) x.dur, data.stimulus.values, 'UniformOutput', false);
stim_name = cellfun(@(x) x.name, data.stimulus.values, 'UniformOutput', false);
% generate a time index array, extract data and time points 100 samples
% before the first stimulus, and 300 samples after the last stimulus
onset_min = min(cellfun(@(x) min(x), onsets));
onset_max = max(cellfun(@(x) max(x), onsets));
time_idx = (find(data.time == onset_min))-100 : (find(data.time == onset_max))+300;
time_vec = data.time(time_idx);     % generate the time vector with time_idx

% data preparation (Nr_chan + 24: corresponding 850nm electrode for HbO)
data2plot = [data.data(time_idx, Nr_chan), data.data(time_idx, Nr_chan + 24)];

% IIR filtered data
[b_l, a_l] = butter(IIR_order, lcut/nqfreq,"low");  % low-pass filter design
data2plot_lf = filtfilt(b_l, a_l, data2plot);   % low-pass filtering
[b_h, a_h] = butter(IIR_order, hcut/nqfreq,"high"); % high-pass filter design
data2plot_hf = filtfilt(b_h, a_h, data2plot);   % high-pass filtering
data2plot_bf = filtfilt(b_h, a_h, data2plot_lf);% band-pass filtering

IIR_array = {data2plot, data2plot_lf, data2plot_hf, data2plot_bf};

% FIR filtered data
FIR_bpFilter = designfilt('bandpassfir', 'FilterOrder', FIR_order, ...
                    'CutoffFrequency1', hcut, 'CutoffFrequency2', lcut, ...
                    'SampleRate', fs);

FIR_array = filtfilt(FIR_bpFilter, data2plot);


% % following lines for online-data
% FIR_array = filter(IIR_bpFilter, data2plot);
% N = length(data2plot);
% grpdelay(IIR_bpFilter, N,fs);
% delay = mean(grpdelay(IIR_bpFilter));
% tt = time_vec(1:end-delay);
% sn = data2plot(1:end-delay,:);
% sf = FIR_array;
% sf(1:delay) = [];
% plot(tt,sn), hold on, plot(time_vec,xn(:,1))


%%
% define ylim for plot
ymax = max(data2plot, [], 'all');
ymin = min(data2plot, [], 'all');

patch_handles = [];
line_handles = [];
line_colors = {'#0072BD','#D95319'};
patch_colors = {'#AAAAAA', '#333333'};

switch type
    case 'IIR'
        figure;

        fig_size = [0, 0, .60, .80];    % size of the whole figure
        set(gcf,'unit', 'normalized','Position', fig_size)
        
        patch_handles = [];

        % plot settings
        pos = { [0.10, 0.93, 0.05, 0.05];
            [0.55, 0.93, 0.05, 0.05];
            [0.10, 0.48, 0.05, 0.05];
            [0.55, 0.48, 0.05, 0.05]};
        pos_nr = {'A','B','C','D'};
        plot_title = {  'Unfiltered';
            'low-pass filtered [0.09 Hz]';
            'high-pass filtered [0.01 Hz]';
            'band-pass filtered [0.01, 0.09] Hz'};

        for i = 1:4
            subplot(2,2,i)
            % line1 = plot(time_vec(:), IIR_array{i}(:,1),'Color',line_colors{1});
            % hold on
            % line2 = plot(time_vec(:), IIR_array{i}(:,2),'Color',line_colors{2});
            line1 = plot(time_vec(ceil(length(time_vec)/2):end), IIR_array{i}(ceil(length(time_vec)/2):end,1),'Color',line_colors{1});
            hold on
            line2 = plot(time_vec(ceil(length(time_vec)/2):end), IIR_array{i}(ceil(length(time_vec)/2):end,2),'Color',line_colors{2});

            ylabel('Δ[HbX] (\muM)')
            title(plot_title{i})
            ylim([ymin - 10, ymax + 10])
            if i == 3 || i == 4
                xlabel('Time (s)')
            end
            yline(0, '--', 'Color', [0.5 0.5 0.5],'LineWidth',1);


            ticks_x = [time_vec(ceil(length(time_vec)/2)), time_vec(end)];
            ticks_y = [ceil(ymin / 10) * 10, 0, floor(ymax / 10) * 10];
            set(gca, 'XTick', ticks_x, ...
                'XTickLabel', {sprintf('%.1f', ticks_x(1)), sprintf('%.1f', ticks_x(2))}, ...
                'YTick', ticks_y, ...
                'YTickLabel', {sprintf('%.1f', ticks_y(1)), sprintf('%.1f', ticks_y(2)), sprintf('%.1f', ticks_y(3))}, ...
                'FontSize', 12);

            ax = gca;
            % ax.XLim = [time_vec(1), time_vec(end)];
            ax.XLim = [time_vec(ceil(length(time_vec)/2)), time_vec(end)];
            ax.XAxisLocation = 'bottom';
            ax.YAxis.TickLength = [0.01, 0.01];
            
            % for k = 1:2
            for k = 2
                for j = 1:length(onsets{k})
                    x = [onsets{k}(j), onsets{k}(j) + stim_dur{k}(j), onsets{k}(j) + stim_dur{k}(j), onsets{k}(j)];
                    y = ylim;
                    y = [y(1), y(1), y(2), y(2)];
                    p1 = patch(x, y, [0.8 0.8 0.8], 'FaceColor',patch_colors{k},'EdgeColor', 'none','FaceAlpha', 0.3);
                    if i== 1 && j == 1
                        patch_handles = [patch_handles, p1];
                    end
                end
            end

            % serial number
            annotation('textbox', pos{i}, 'String', pos_nr{i}, 'FontSize', 14, ...
                'FontWeight', 'bold', 'EdgeColor', 'none');

            if i == 1
                line_handles = [line1, line2];
            end

        end

        fig_handles = [line_handles, patch_handles];
        % lgd = legend(fig_handles, {'HbR (760nm)', 'HbO (850nm)', 'ME-LEFT', 'ME-RIGHT'}, ...
        %     'Location', 'southoutside', 'Orientation', 'horizontal');
        lgd = legend(fig_handles, {'HbR (760nm)', 'HbO (850nm)', 'ME-RIGHT'}, ...
            'Location', 'southoutside', 'Orientation', 'horizontal');

        lgd.Position = [0.35, 0.0, 0.3, 0.05];
        lgd.Box = 'off';
        set(gca, 'LooseInset', get(gca, 'TightInset') + [0, 0.35, 0, 0.55]);

    %%
    case 'IIR_FIR'
        figure;

        fig_size = [0, 0, .60, .80];    % size of the whole figure
        set(gcf,'unit', 'normalized','Position', fig_size)
        fig_handles = [];
        line_handles = [];
        patch_handles = [];

        % numbering settings
        pos = { 
            [0.10, 0.93, 0.05, 0.05];
            [0.10, 0.48, 0.05, 0.05]
            };
        pos_nr = {'A','B'};
        plot_title = {'FIR and IIR Band-pass Filter (0.01~0.09 Hz) on HbR (760nm)';
            'FIR and IIR Band-pass Filter (0.01~0.09 Hz) on HbO (850nm)'};

        % subplots positions: [left, bottom, width, height]
        subplot_pos = {
            [0.1, 0.55, 0.7, 0.35];  
            [0.1, 0.1, 0.7, 0.35] 
            };

        for i = 1:2
            fig = subplot(2,1,i);
            set(fig, 'Position', subplot_pos{i});
            % handle_tem = plot(time_vec,data2plot(:,i),'Color',line_colors{i});
            handle_tem = plot(time_vec,data2plot(:,i),'Color',line_colors{i});
            line_handles = [line_handles, handle_tem];

            hold on
            % FIR_handle = plot(time_vec,FIR_array(:,i),'-','Color','#4F726C','linewidth',1.5);
            % IIR_handle = plot(time_vec, data2plot_bf(:,i), '-','Color','#0B1013','LineWidth',1.5);
            FIR_handle = plot(time_vec(ceil(length(time_vec)/2):end),FIR_array(ceil(length(time_vec)/2):end,i),'-','Color','#4F726C','linewidth',1.5);
            IIR_handle = plot(time_vec(ceil(length(time_vec)/2):end), data2plot_bf(ceil(length(time_vec)/2):end,i), '-','Color','#0B1013','LineWidth',1.5);

            if i == 2
                xlabel('Time (s)')
            end
            ylabel('Δ[HbX] (\muM)')
            ylim([ymin - 10, ymax + 10])
            
            title(plot_title{i})

            % ticks_x = [time_vec(1), time_vec(end)];
            ticks_x = [time_vec(ceil(length(time_vec)/2)), time_vec(end)];
            ticks_y = [ceil(ymin / 10) * 10, 0, floor(ymax / 10) * 10];
            set(gca, 'XTick', ticks_x, ...
                'XTickLabel', {sprintf('%.1f', ticks_x(1)), sprintf('%.1f', ticks_x(2))}, ...
                'YTick', ticks_y, ...
                'YTickLabel', {sprintf('%.1f', ticks_y(1)), sprintf('%.1f', ticks_y(2)), sprintf('%.1f', ticks_y(3))}, ...
                'FontSize', 12);

            ax = gca;
            % ax.XLim = [time_vec(1), time_vec(end)];
            ax.XLim = [time_vec(ceil(length(time_vec)/2)), time_vec(end)];
            ax.XAxisLocation = 'bottom';
            ax.YAxis.TickLength = [0.01, 0.01];

            % for k = 1:2
            for k = 2
                for j = 1:length(onsets{k})
                    x = [onsets{k}(j), onsets{k}(j) + stim_dur{k}(j), onsets{k}(j) + stim_dur{k}(j), onsets{k}(j)];
                    y = ylim;
                    y = [y(1), y(1), y(2), y(2)];
                    p1 = patch(x, y, [0.8 0.8 0.8], 'FaceColor',patch_colors{k},'EdgeColor', 'none','FaceAlpha', 0.3);
                    if i== 1 && j == 1
                        patch_handles = [patch_handles, p1];
                    end
                end
            end

        end
        fig_handles = [line_handles, FIR_handle, IIR_handle, patch_handles];
        % lgd = legend(fig_handles, {'HbR (760nm)', 'HbO (850nm)', ...
        %     'Filtered Signal(FIR)','Filtered Signal(IIR)','ME-LEFT', 'ME-RIGHT'}, ...
        %     'Location', 'eastoutside', 'Orientation', 'vertical');
        lgd = legend(fig_handles, {'HbR (760nm)', 'HbO (850nm)', ...
            'FIR Filtered','IIR Filtered','ME-RIGHT'}, ...
            'Location', 'eastoutside', 'Orientation', 'vertical');

        lgd.Position = [0.90, 0.50, 0.05, 0.05];
        lgd.Box = 'off';
        set(gca, 'LooseInset', get(gca, 'TightInset') + [0.3, 0.35, 0.3, 0.55]);


end
