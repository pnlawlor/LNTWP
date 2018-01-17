%% Load scripts

addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Scripts\Pat'))

%% Load data

fpath = 'C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Source_data\Miller_lab\Processed\Brian_uncertainty\';
load([fpath 'Brian_uncertainty_correct_only.mat'])


%% User inputs

%% Initialize

%% Do it

for block_num = length(Data)
    for trial_num = length(Data(block_num).kinematics);
        dir_reach = get_reach_dir(kinematics{trial_num});
        wonky_reach = is_reach_wonky(kinematics{trial_num});
    end
end