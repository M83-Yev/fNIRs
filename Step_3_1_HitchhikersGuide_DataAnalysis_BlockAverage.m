% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%%            A Hitchhiker's Guide to fNIRS Data Analysis                %%
%                 - DATA ANALYSIS: BLOCK AVERAGING (Manuscript Section 3.3.1) -                      %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Created by Franziska Klein
% Commented by Mojtaba Soltanlou
% Using Matlab 2021b
% Last Updated 24/08/2023

% What you need:

% NIRS Brain AnalyzIR Toolbox: https://github.com/huppertt/nirs-toolbox
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

close all; clear all; clc;
% Change below to your project folder path.
% Here are the directories that you need to set, which will be used in the next steps. Note that the input 
% file is the output file of the previous script (preprocessing).
MAINPATH = ['your_path_here’ , ‘Hitchhikers_Guide’]; 
PATHIN      = [MAINPATH, 'DATA_OUT\preprocessed\'];
PATHOUT  = [MAINPATH, 'DATA_OUT\block_average_ana\'];

% List of the variables. The epoch limit is the time window for the analysis (i.e., averaging), which is 
% from 5 s before the block onset (i.e., baseline) to 25 s after that. This could be changed based on 
% your design or needs. For instance, you can select a longer or shorter baseline, but be careful that 
% your baseline doesn’t include the previous block. A practical way to find the most appropriate 
% baseline and even the length of your block or long event (especially if it is a self-paced design with 
% various lengths of the events; e.g., Soltanlou et al., 2022), is to plot your grand average (average of 
% all channels and all participants). You can then decide where the signal is on the baseline (for the 
% length of the baseline) and when it comes back to the baseline after an oscillation (for the length of 
% the block or long event).
epoch_lim = [-5, 25];
% Name of the conditions in the current data set. You need to change it based on your conditions.
stim = {'ME_LEFT', 'ME_RIGHT'};
% Plot positions according to the channel layout.
idx = [5, 13, 14, 19, 21, 29, 15, 23, 25, 33, 27, 31, 39, 41, 35, 43];
% Source-Detector pairs for ROI ME LEFT.
opt_M1_LEFT = [3 3; 3 4; 3 6; 6 4; 6 6];
% Source-Detector pairs for ROI ME RIGHT.
opt_M1_RIGHT = [4 5; 4 7; 5 5; 7 5; 7 7];
% This is in seconds.
start_val = 6;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~ Load Preprocessed Data ~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

load([PATHIN, 'data_preprocessed.mat'], 'Hb_SDCcor');
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~ Normalize Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ % 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Add description to the data. 
for sub = 1:size(Hb_SDCcor, 1)
    Hb_SDCcor(sub).data = zscore(Hb_SDCcor(sub).data);
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~ Run Block Average ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

job = nirs.modules.BlockAverage;
% Define the pre-stimulus duration (i.e., baseline) in seconds. This refers to the first value in the above 
% defined epoch duration.
job.pre = abs(epoch_lim(1));
% Define the post-stimulus duration (i.e., length of the block or event) in seconds. This refers to the 
% second value in the above defined epoch duration.
job.post = epoch_lim(2);
job.stim_names = stim;

% You may decide to apply no baseline correction. This needs a very good justification as the 
% uncorrected block contaminates stimulus unrelated changes that it inherited from the pre-stimulus or 
% baseline duration.
% job.baseline = 'none';
HbX_BA = job.run(Hb_SDCcor);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~ VISUAL INSPECTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% To visually inspect data before taking the next steps, select HbO and HbR of LDC. This is not a part 
% of the analysis process, but rather a step to double check the processing signal.
ch_HbO = find(strcmp(Hb_SDCcor(1).probe.link.type, 'hbo') & Hb_SDCcor(1).probe.link.ShortSeperation == 0);
ch_HbR = find(strcmp(Hb_SDCcor(1).probe.link.type, 'hbr') & Hb_SDCcor(1).probe.link.ShortSeperation == 0);

t = linspace(epoch_lim(1), epoch_lim(2), size(HbX_BA.stim_ME_LEFT_allSubj_BA, 1));

% ME LEFT
figure;
for ch = 1:length(ch_HbO)
    subplot(5, 9, idx(ch))
    hold on
    plot(t, HbX_BA.stim_ME_LEFT_allSubj_BA(:, ch_HbO(ch)), 'Color', [248, 105, 35]./255, 'LineWidth', 4);
    plot(t, HbX_BA.stim_ME_LEFT_allSubj_BA(:, ch_HbR(ch)), 'Color', [102, 134, 176]./255, 'LineWidth', 4);
    xlim([epoch_lim(1), epoch_lim(2)]);
    ylim([min(min(HbX_BA.stim_ME_LEFT_allSubj_BA(:, [ch_HbO, ch_HbR]))), ...
        max(max(HbX_BA.stim_ME_LEFT_allSubj_BA(:, [ch_HbO, ch_HbR])))]);
    hline(0, 'k');
    % Stimulus onset:
    vline(0, 'k');
    % Stimulus offset:
    vline(15, 'k');
    if idx(ch) == 5
        title('BA ME LEFT');
    end
end

% ME RIGHT
figure;
for ch = 1:length(ch_HbO)
    subplot(5, 9, idx(ch))
    hold on
    plot(t, HbX_BA.stim_ME_RIGHT_allSubj_BA(:, ch_HbO(ch)), 'Color', [248, 105, 35]./255, 'LineWidth', 4);
    plot(t, HbX_BA.stim_ME_RIGHT_allSubj_BA(:, ch_HbR(ch)), 'Color', [102, 134, 176]./255, 'LineWidth', 4);
    xlim([epoch_lim(1), epoch_lim(2)]);
    ylim([min(min(HbX_BA.stim_ME_RIGHT_allSubj_BA(:, [ch_HbO, ch_HbR]))), ...
        max(max(HbX_BA.stim_ME_RIGHT_allSubj_BA(:, [ch_HbO, ch_HbR])))]);
    hline(0, 'k');
    % Stimulus onset:
    vline(0, 'k');
    % Stimulus offset:
    vline(15, 'k');
    if idx(ch) == 5
        title('BA ME RIGHT');
    end
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~ Channel Selection ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% We don’t always use all the recorded channels, but rather some channels that make our region of 
% interest (ROI). In the current data, select the channels for the ROI M1 LEFT.
hbo_M1_LEFT = find(ismember(Hb_SDCcor(1).probe.link.source, opt_M1_LEFT(:, 1)) & ...
                   ismember(Hb_SDCcor(1).probe.link.detector, opt_M1_LEFT(:, 2)) & ...
                   strcmp(Hb_SDCcor(1).probe.link.type, 'hbo'));
hbr_M1_LEFT = find(ismember(Hb_SDCcor(1).probe.link.source, opt_M1_LEFT(:, 1)) & ...
                   ismember(Hb_SDCcor(1).probe.link.detector, opt_M1_LEFT(:, 2)) & ...
                   strcmp(Hb_SDCcor(1).probe.link.type, 'hbr'));

% Now select the channels for the ROI M1 RIGHT.
hbo_M1_RIGHT = find(ismember(Hb_SDCcor(1).probe.link.source, opt_M1_RIGHT(:, 1)) & ...
                   ismember(Hb_SDCcor(1).probe.link.detector, opt_M1_RIGHT(:, 2)) & ...
                   strcmp(Hb_SDCcor(1).probe.link.type, 'hbo'));
hbr_M1_RIGHT = find(ismember(Hb_SDCcor(1).probe.link.source, opt_M1_RIGHT(:, 1)) & ...
                   ismember(Hb_SDCcor(1).probe.link.detector, opt_M1_RIGHT(:, 2)) & ...
                   strcmp(Hb_SDCcor(1).probe.link.type, 'hbr'));

% You could be more selective and instead of a group of channels (i.e., ROI), select individual best 
% channel based on max/min peak in each recording region.
% To find the best channel for each chromophore and stimulus type, follow the below steps.

for sub = 1:size(HbX_BA.stim_ME_LEFT_indiv_BA, 1);
   

    % HbO ME LEFT:
    [~, tmp] = max(max(squeeze(HbX_BA.stim_ME_LEFT_indiv_BA(sub, :, hbo_M1_RIGHT))));
    bestChan_HbO.ME_LEFT.M1_RIGHT(sub) = hbo_M1_RIGHT(tmp);

    [~, tmp] = max(max(squeeze(HbX_BA.stim_ME_LEFT_indiv_BA(sub, :, hbo_M1_LEFT))));
    bestChan_HbO.ME_LEFT.M1_LEFT(sub) = hbo_M1_LEFT(tmp);

    % HbO ME RIGHT:
    [~, tmp] = max(max(HbX_BA.stim_ME_RIGHT_indiv_BA(sub, :, hbo_M1_LEFT)));
    bestChan_HbO.ME_RIGHT.M1_LEFT(sub) = hbo_M1_LEFT(tmp);

    [~, tmp] = max(max(HbX_BA.stim_ME_RIGHT_indiv_BA(sub, :, hbo_M1_RIGHT)));
    bestChan_HbO.ME_RIGHT.M1_RIGHT(sub) = hbo_M1_RIGHT(tmp);

    % HbR ME LEFT:
    [~, tmp] = min(min(squeeze(HbX_BA.stim_ME_LEFT_indiv_BA(sub, :, hbr_M1_RIGHT))));
    bestChan_HbR.ME_LEFT.M1_RIGHT(sub) = hbr_M1_RIGHT(tmp);

    [~, tmp] = min(min(squeeze(HbX_BA.stim_ME_LEFT_indiv_BA(sub, :, hbr_M1_LEFT))));
    bestChan_HbR.ME_LEFT.M1_LEFT(sub) = hbr_M1_LEFT(tmp);

    % HbR ME RIGHT:
    [~, tmp] = min(min(HbX_BA.stim_ME_RIGHT_indiv_BA(sub, :, hbr_M1_LEFT)));
    bestChan_HbR.ME_RIGHT.M1_LEFT(sub) = hbr_M1_LEFT(tmp);

    [~, tmp] = min(min(HbX_BA.stim_ME_RIGHT_indiv_BA(sub, :, hbr_M1_RIGHT)));
    bestChan_HbR.ME_RIGHT.M1_RIGHT(sub) = hbr_M1_RIGHT(tmp);
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~ DATA EXTRACTION - Best Channel ~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% If you have decided to select the best channel in each region, you can extract those channels here.
start = round(start_val*Hb_SDCcor(1).Fs);

for sub = 1:length(bestChan_HbO.ME_LEFT.M1_LEFT)

    % To extract individual average across individual BA for ME LEFT:
    ME_LEFT.HbO.M1_LEFT(sub)  = mean(squeeze(HbX_BA.stim_ME_LEFT_indiv_BA(sub, start:end, bestChan_HbO.ME_LEFT.M1_LEFT(sub))));
    ME_LEFT.HbO.M1_RIGHT(sub) = mean(squeeze(HbX_BA.stim_ME_LEFT_indiv_BA(sub, start:end, bestChan_HbO.ME_LEFT.M1_RIGHT(sub))));

    ME_LEFT.HbR.M1_LEFT(sub)  = mean(squeeze(HbX_BA.stim_ME_LEFT_indiv_BA(sub, start:end, bestChan_HbR.ME_LEFT.M1_LEFT(sub))));
    ME_LEFT.HbR.M1_RIGHT(sub) = mean(squeeze(HbX_BA.stim_ME_LEFT_indiv_BA(sub, start:end, bestChan_HbR.ME_LEFT.M1_RIGHT(sub))));

    % To extract individual average across individual BA for ME RIGHT:
    ME_RIGHT.HbO.M1_LEFT(sub)  = mean(squeeze(HbX_BA.stim_ME_RIGHT_indiv_BA(sub, start:end, bestChan_HbO.ME_RIGHT.M1_LEFT(sub))));
    ME_RIGHT.HbO.M1_RIGHT(sub) = mean(squeeze(HbX_BA.stim_ME_RIGHT_indiv_BA(sub, start:end, bestChan_HbO.ME_RIGHT.M1_RIGHT(sub))));

    ME_RIGHT.HbR.M1_LEFT(sub)  = mean(squeeze(HbX_BA.stim_ME_RIGHT_indiv_BA(sub, start:end, bestChan_HbR.ME_RIGHT.M1_LEFT(sub))));
    ME_RIGHT.HbR.M1_RIGHT(sub) = mean(squeeze(HbX_BA.stim_ME_RIGHT_indiv_BA(sub, start:end, bestChan_HbR.ME_RIGHT.M1_RIGHT(sub))));
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Save Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Check if the path where data is saved exists.
if(~exist(PATHOUT, ‘dir’))
	mkdir(PATHOUT);
end

save([PATHOUT, 'BA_values.mat'], 'ME_LEFT', 'ME_RIGHT', 'bestChan_HbO', 'bestChan_HbR');
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~ STATISTICS  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~ Descriptive Statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% This section will save the descriptive statistics that are usually needed to be reported.

tab_hbo = table([length(ME_LEFT.HbO.M1_LEFT); mean(ME_LEFT.HbO.M1_LEFT); std(ME_LEFT.HbO.M1_LEFT); ...
            std(ME_LEFT.HbO.M1_LEFT)\sqrt(length(ME_LEFT.HbO.M1_LEFT)); min(ME_LEFT.HbO.M1_LEFT); max(ME_LEFT.HbO.M1_LEFT)], ...
            [length(ME_LEFT.HbO.M1_RIGHT); mean(ME_LEFT.HbO.M1_RIGHT); std(ME_LEFT.HbO.M1_RIGHT); ...
            std(ME_LEFT.HbO.M1_RIGHT)\sqrt(length(ME_LEFT.HbO.M1_RIGHT)); min(ME_LEFT.HbO.M1_RIGHT); max(ME_LEFT.HbO.M1_RIGHT)], ...
            [length(ME_RIGHT.HbO.M1_LEFT); mean(ME_RIGHT.HbO.M1_LEFT); std(ME_RIGHT.HbO.M1_LEFT); ...
            std(ME_RIGHT.HbO.M1_LEFT)\sqrt(length(ME_RIGHT.HbO.M1_LEFT)); min(ME_RIGHT.HbO.M1_LEFT); max(ME_RIGHT.HbO.M1_LEFT)], ...
            [length(ME_RIGHT.HbO.M1_RIGHT); mean(ME_RIGHT.HbO.M1_RIGHT); std(ME_RIGHT.HbO.M1_RIGHT); ...
            std(ME_RIGHT.HbO.M1_RIGHT)\sqrt(length(ME_RIGHT.HbO.M1_RIGHT)); min(ME_RIGHT.HbO.M1_RIGHT); max(ME_RIGHT.HbO.M1_RIGHT)], ...
            'VariableNames', {'ME_LEFT_M1_LEFT','ME_LEFT_M1_RIGHT', 'ME_RIGHT_M1_LEFT','ME_RIGHT_M1_RIGHT'}, ...
            'RowNames', {'N', 'MEAN', 'STD', 'STDERR', 'MIN', 'MAX'});

tab_hbr = table([length(ME_LEFT.HbR.M1_LEFT); mean(ME_LEFT.HbR.M1_LEFT); std(ME_LEFT.HbR.M1_LEFT); ... 
                std(ME_LEFT.HbR.M1_LEFT)\sqrt(length(ME_LEFT.HbR.M1_LEFT)); min(ME_LEFT.HbR.M1_LEFT); max(ME_LEFT.HbR.M1_LEFT)], ...
            [length(ME_LEFT.HbR.M1_RIGHT); mean(ME_LEFT.HbR.M1_RIGHT); std(ME_LEFT.HbR.M1_RIGHT); ...
            std(ME_LEFT.HbR.M1_RIGHT)\sqrt(length(ME_LEFT.HbR.M1_RIGHT)); min(ME_LEFT.HbR.M1_RIGHT); max(ME_LEFT.HbR.M1_RIGHT)], ...
            [length(ME_RIGHT.HbR.M1_LEFT); mean(ME_RIGHT.HbR.M1_LEFT); std(ME_RIGHT.HbR.M1_LEFT); ...
            std(ME_RIGHT.HbR.M1_LEFT)\sqrt(length(ME_RIGHT.HbR.M1_LEFT)); min(ME_RIGHT.HbR.M1_LEFT); max(ME_RIGHT.HbR.M1_LEFT)], ...
            [length(ME_RIGHT.HbR.M1_RIGHT); mean(ME_RIGHT.HbR.M1_RIGHT); std(ME_RIGHT.HbR.M1_RIGHT); ...
            std(ME_RIGHT.HbR.M1_RIGHT)\sqrt(length(ME_RIGHT.HbR.M1_RIGHT)); min(ME_RIGHT.HbR.M1_RIGHT); max(ME_RIGHT.HbR.M1_RIGHT)], ...
            'VariableNames', {'ME_LEFT_M1_LEFT','ME_LEFT_M1_RIGHT', 'ME_RIGHT_M1_LEFT','ME_RIGHT_M1_RIGHT'}, ...
            'RowNames', {'N', 'MEAN', 'STD', 'STDERR', 'MIN', 'MAX'});
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Inferential Statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Following the current design, we run repeated measure analysis of variance (rmANOVA). You need 
% to adjust the statistical analysis based on your hypothesis testing and/or your study design.
% The results of the HbO analysis.
HEMI = [repmat({'M1_LEFT'}, length(ME_LEFT.HbO.M1_LEFT), 1); repmat({'M1_RIGHT'}, length(ME_LEFT.HbO.M1_LEFT), 1)];
TASK(:, 1) = [ME_LEFT.HbO.M1_LEFT'; ME_LEFT.HbO.M1_RIGHT'];
TASK(:, 2) = [ME_RIGHT.HbO.M1_LEFT'; ME_RIGHT.HbO.M1_RIGHT'];

t = table(HEMI, TASK(:, 1), TASK(:, 2),...
'VariableNames', {'HEMI', 'ME_LEFT', 'ME_RIGHT'});
task = table([1 2]', 'VariableNames', {'TASK'});

rmANAOVA_hbo = fitrm(t, 'ME_LEFT-ME_RIGHT~HEMI', 'WithinDesign', task);

rmANAOVA_hbo = ranova(rmANAOVA_hbo);

% The results of HbR analysis.
HEMI = [repmat({'M1_LEFT'}, length(ME_LEFT.HbR.M1_LEFT), 1); repmat({'M1_RIGHT'}, length(ME_LEFT.HbR.M1_LEFT), 1)];
TASK(:, 1) = [ME_LEFT.HbR.M1_LEFT'; ME_LEFT.HbR.M1_RIGHT'];
TASK(:, 2) = [ME_RIGHT.HbR.M1_LEFT'; ME_RIGHT.HbR.M1_RIGHT'];

t = table(HEMI, TASK(:, 1), TASK(:, 2),...
'VariableNames', {'HEMI', 'ME_LEFT', 'ME_RIGHT'});
task = table([1 2]', 'VariableNames', {'TASK'});

rmANAOVA_hbr = fitrm(t, 'ME_LEFT-ME_RIGHT~HEMI', 'WithinDesign', task);

rmANAOVA_hbr = ranova(rmANAOVA_hbr);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% You can use Gramm to visualize your results. It is a complete data visualization toolbox for 
% MATLAB that provides an easy-to-use and high-level interface to produce publication-quality 
% plots of complex data with varied statistical visualizations. 
% https://github.com/piermorel/gramm



