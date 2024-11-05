% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%%            A Hitchhiker's Guide to fNIRS Data Analysis                %%
%                        - DATA ANALYSIS: GLM (Manuscript Section 3.3.2) -                           %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Created by Franziska Klein
% Commented by Mojtaba Soltanlou
% Using Matlab 2021b
% Last Updated 25/08/2023

% What you need:

% NIRS Brain AnalyzIR Toolbox: https://github.com/huppertt/nirs-toolbox
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

close all; clear all; clc;

% Change below to your project folder path.
% Here are the directories that you need to set, which will be used in the next steps. Note that the input 
% file is the output file of the previous script (preprocessing).
MAINPATH = ['your_path_here’, ‘Hitchhikers_Guide’];
PATHIN   = [MAINPATH, 'DATA_OUT\preprocessed\'];
PATHOUT  = [MAINPATH, 'DATA_OUT\glm_ana\'];
if ~exist(PATHOUT, 'dir')
    mkdir(PATHOUT)
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~ Load Preprocessed Data ~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

load([PATHIN, 'data_preprocessed.mat'], 'Hb_SDCcor');
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~ Data Normalization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Add description to the data.
for sub = 1:size(Hb_SDCcor, 1)
    Hb_SDCcor(sub).data = zscore(Hb_SDCcor(sub).data);
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~ Run General Linear Model (GLM) ~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% There are a few different ways of running the General Linear Model (GLM) on the fNIRS data. We will 
% present two commonly used models here.

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% GLM with Canonical HRF (Manuscript Section 3.3.3)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% In case of a priori assumptions about the shape of the hemodynamic response function (HRF),
% we chose a canonical HRF with a peak at 6 s as the basis function.
basis = nirs.design.basis.Canonical;

% Here are further possible options and their default values that you may want to adjust based on 
% your needs.
% The HRF peaks at 6 s after stimulus onset.
basis.peakTime = 6;
 % To add the first and the second derivatives for each regressor to account for time and dispersion.
basis.incDeriv = true;

job = nirs.modules.GLM;
% Chose the GLM type. The default is AR-IRLS (Barker et al., 2016). The other options are: 'OLS', 
% 'NIRS-SPM','MV-GLM', and 'Nonlinear'.
job.type = 'AR-IRLS';
job.basis('default') = basis;
% Add the constant for trend to the model. The other options are legendre and dctmtx.
job.trend_func = @(t) nirs.design.trend.constant(t);

SubjStats_Canonical = job.run(Hb_SDCcor);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% GLM with FIR (Manuscript Section 3.3.4)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% This is a flexible model with no assumptions for the HRF; basis function is therefore a finite impulse 
% response (FIR) function. 

% In the current data, we have 15 s task, followed by 10 s rest. With a sampling frequency of ~10 Hz 
% this would result in (15+10) s x 10 Hz = 250 bins or data points in the output file. When the model is 
% estimated, there would be 250 corresponding coefficients β that model the amplitude of the response 
% over this entire time period from stimulus onset until 10 s after the end of the task. % Because of the 
% size of the data, running this GLM model would take forever and therefore we increased the bin width 
% to 5 by grouping 5 data points together, which results in 250 / 5 = 50 bins and accordingly 50 
% corresponding coefficients β. This means one β describes ~0.5 s of the signal.

basis = nirs.design.basis.FIR;
% Each bin is 5 samples wide.
basis.binwidth = 10;
basis.nbins = 25;

job = nirs.modules.GLM;

% Chose the GLM type. The default is AR-ILS (Barker et al., 2016). The other options are: 'OLS', 
% 'NIRS-SPM','MV-GLM', and 'Nonlinear'.
job.type = 'AR-IRLS';
job.basis('default') = basis;

% Add the constant for trend to the model. The other options are legendre and dctmtx.
job.trend_func = @(t) nirs.design.trend.constant(t);

SubjStats_FIR = job.run(Hb_SDCcor);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~ Remove Outlier Subjects~~~ ~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%This optional step allows you to remove outlier subjects. To run this step, you need to  make yourself familiar with setting appropriate outlier criteria.

% job = nirs.modules.RemoveOutlierSubjects;
% 
% ~~~~~~~~~~~ %
% Canonical model
% ~~~~~~~~~~~ %
% SubjStats_Canonical_Pruned = job.run(SubjStats_Canonical);
% tbl_Canonical = nirs.util.grouplevelleveragestats(SubjStats_Canonical);


% ~~~~~~ %
% FIR model
% ~~~~~~ %

% SubjStats_FIR_Pruned = job.run(SubjStats_FIR);
% tbl_FIR = nirs.util.grouplevelleveragestats(SubjStats_FIR);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~ Save Channel Statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% This step is optional due to the large size of the output files and the time it takes.

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Canonical model
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
 save([PATHOUT, 'GLM_Channel_Subject_Stats_Canonical.mat'], 'SubjStats_Canonical');

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% FIR model
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
 save([PATHOUT, 'GLM_Channel_Subject_Stats_FIR.mat'], 'SubjStats_FIR', '-v7.3','-nocompression');
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STATISTICS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Variables:
COND      = {'ME_LEFT:01', 'ME_RIGHT:01'};
COND_NAME = {'ME_LEFT', 'ME_RIGHT'};
ROIs      = {'M1_LEFT', 'M1_RIGHT'};

% plotting positions according to channel layout
% add PLOTTING Code

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Group-level Analysis I ~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% We extract the betas to find the largest HbO and smallest HbR betas for each condition and each 
% ROI. You can use this data in MATLAB or export it to your favorit statistics software (e.g., R, JASP).

for roi = 1:size(ROI_opt, 1)
    % To find the indices for ROI channels.
    idx_roi = find(ismember(optodes(:, 1), ROI_opt(roi, :, 1)) & ...
                       ismember(optodes(:, 2), ROI_opt(roi, :, 2)));

    for cond = 1:length(COND)
        % To find the largest betas for HbO.
        max_HbO_betas.(COND_NAME{cond}).(ROIs{roi}) = max(betas_HbO.(COND_NAME{cond})(:, idx_roi), [], 2);

        % To find the smallest betas for HbR.
        min_HbR_betas.(COND_NAME{cond}).(ROIs{roi}) = min(betas_HbR.(COND_NAME{cond})(:, idx_roi), [], 2);
    end
end


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Inferential Statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Following the current design, we run repeated measure analysis of variance (rmANOVA). You need 
% to adjust the statistical analysis based on your hypothesis testing and/or your study design.

% The results of the HbO analysis.
HEMI = [repmat({'M1_LEFT'}, length(max_HbO_betas.ME_LEFT.M1_LEFT), 1); 
repmat({'M1_RIGHT'}, length(max_HbO_betas.ME_LEFT.M1_RIGHT), 1)];
TASK(:, 1) = [max_HbO_betas.ME_LEFT.M1_LEFT; max_HbO_betas.ME_LEFT.M1_RIGHT];
TASK(:, 2) = [max_HbO_betas.ME_RIGHT.M1_LEFT; max_HbO_betas.ME_RIGHT.M1_RIGHT];

t = table(HEMI, TASK(:, 1), TASK(:, 2),...
'VariableNames', {'HEMI', 'ME_LEFT', 'ME_RIGHT'});
task = table([1 2]', 'VariableNames', {'TASK'});

rmANAOVA_hbo = fitrm(t, 'ME_LEFT-ME_RIGHT~HEMI', 'WithinDesign', task);

rmANAOVA_hbo = ranova(rmANAOVA_hbo);

% The results of the HbR analysis.
HEMI = [repmat({'M1_LEFT'}, length(min_HbR_betas.ME_LEFT.M1_LEFT), 1); repmat({'M1_RIGHT'}, length(min_HbR_betas.ME_LEFT.M1_RIGHT), 1)];
TASK(:, 1) = [min_HbR_betas.ME_LEFT.M1_LEFT; min_HbR_betas.ME_LEFT.M1_RIGHT];
TASK(:, 2) = [min_HbR_betas.ME_RIGHT.M1_LEFT; min_HbR_betas.ME_RIGHT.M1_RIGHT];

t = table(HEMI, TASK(:, 1), TASK(:, 2),...
'VariableNames', {'HEMI', 'ME_LEFT', 'ME_RIGHT'});
task = table([1 2]', 'VariableNames', {'TASK'});

rmANAOVA_hbr = fitrm(t, 'ME_LEFT-ME_RIGHT~HEMI', 'WithinDesign', task);

rmANAOVA_hbr = ranova(rmANAOVA_hbr);


% TO DO: Add same stuff for FIR model
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~ Group-level Analysis II (Manuscript Section 3.3.5) ~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% The following runs a second-level GLM analysis with NIRS toolbox functions using mixed 
% model analysis.

job = nirs.modules.MixedEffects;  

% This model follows the Wikinson-Roger's notation. Any demographic information (e.g., age) can be 
% incorporated into the model.
% MixedEffects with properties:
%                 formula: 'beta ~ -1 + cond + (1|subject)'   
%             dummyCoding: 'full'
%              centerVars: 1
%     include_diagnostics: 0
%                  robust: 0
%                weighted: 1
%                 verbose: 1
%                    name: 'Mixed Effects Model'
%                 prevJob: []


% For simplicity, let's just do this one. Note that you can’t add these covariates to the rmANOVA that 
% was explained above.
job.formula = 'beta ~ -1 + cond + (1|sub)';

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Canonical model
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% To find the results of the canonical model.
GroupStatsME_Canonical = job.run(SubjStats_Canonical);
disp(GroupStatsME_Canonical.conditions)
ContrastStats = GroupStatsME_Canonical.ttest({'ME_LEFT:01'; ...
                                              'ME_RIGHT:01'});       
ContrastStats.draw('tstat', [-5 5], 'q<.05')

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% FIR model
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% To find the results of the FIR model.
GroupStatsME_FIR = job.run(SubjStats_FIR);  
disp(GroupStatsME_FIR.conditions)
ContrastStats = GroupStatsME_FIR.ttest({'ME_LEFT:01'; ...
                                              'ME_RIGHT:01'});       
ContrastStats.draw('tstat', [-5 5], 'q<.05')


% We define the ROIs based on source-detector pairs.
M1_left  = [3 4; ...
            3 6; ...
            6 4; ...
            6 6];

M1_right = [4 5; ...
            4 7; ...
            7 5; ...
            7 7];

ROIs{1} = table(M1_left(:, 1),  M1_left(:, 2),  'VariableNames', {'source', 'detector'});
ROIs{2} = table(M1_right(:, 1), M1_right(:, 2), 'VariableNames', {'source', 'detector'});

ROItable = nirs.util.roiAverage(ContrastStats, ROIs{1}, {'M1 LEFT'});
disp(ROItable)

ROItable = nirs.util.roiAverage(ContrastStats, ROIs{2}, {'M1 RIGHT'});
disp(ROItable)

% TO DO: Add same stuff for FIR model


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~ BRAIN VISUALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Follow these steps to plot the data on a brain.
% To use an actual MRI brain, register and download the mesh
% http://mcx.sourceforge.net/cgi-bin/index.cgi?MMC/Colin27AtlasMesh
% Add the file "MMC_Collins_Atlas_Mesh_Version_2L.mat" to the 
% MATLAB path:
% addpath('yourpath\MMC_Collins');
% In addition, you need to download and install the iso2mesh package and add it to the MATLAB path.
% http://iso2mesh.sourceforge.net

% 'C:\users\yevge\Downloads\Compressed\MCXStudio-win64-v2024.2\MCXStudio\MATLAB\mcxlab\examples\colin27_v3.mat'

addpath(genpath('yourpath\iso2mesh-1.9.6'));

for i = 1:length(data)
lambda = unique(data(i).probe.link.type);
fwdBEM = nirs.registration.Colin27.BEM(lambda);
end

% By default, the mesh will draw the fiducial points when plotting. This is controlled in the 
% mesh(1).fiducials field. To turn off all of the 10-20 labels stop the function.
fwdBEM.mesh(1).fiducials.Draw(:) = false;
fwdBEM.mesh(1).transparency=0;
fwdBEM.mesh(2).transparency=0;
fwdBEM.mesh(3).transparency=0;
fwdBEM.draw;

% This command will get the head shape from the mesh.
headshape = nirs.registration.getheadshape(fwdBEM.mesh(1));

% Likewise, this will register a mesh onto your probe.
for i = 1:length(data)
data(i).probe = data(i).probe.register_mesh2probe(fwdBEM.mesh);
end

% You can also draw depth maps of the probe or 10-20 space. 
% This command will display a list of all the available labels.
disp(nirs.util.depthmap);

% This will plot the probe and depth of the nearest cortical point in the
% ROI to the surface of the head. Note that a depth greater than 30 mm is inaccessible to NIRS.
nirs.util.depthmap('BA-4',rawEditc(1).probe);

% You can now plot your data.
GroupLevelStats.probe.defaultdrawfcn='3D mesh (superior)';
GroupLevelStats.draw('tstat', [-5 5], 'q < 0.05');
resultsdir = 'yourpath\RESULTS\GroupLevelStats\';
folder = [resultsdir filesep 'figures'];
GroupLevelStats.printAll('tstat', [-5 5], 'q < 0.05', folder,'png')
