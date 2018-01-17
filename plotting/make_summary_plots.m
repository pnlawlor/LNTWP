%% Load scripts

addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Scripts\Pat'))
addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Projects\FEF\Scripts\Pavan_GLM_code_v2\FEF\utilities'))
addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Scripts_TW'))

%% Simulations: Load data

% Summary stats for Warp Condition x Num Neurons
sim_results_fpath = 'C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Simulation_data\Results\SummaryStats\';
sim_results_fname = 'SummaryStats_warptype_x_numnrn_all_priors.mat';

load([sim_results_fpath sim_results_fname])
% THIS DOESN'T WORK BECAUSE SAVING THE SUMMARYSTATS FILE DIDN'T WORK. 
% So just regenerate the file using 'organize_results.m'

% Actually, you need to run "organize_results.m" to generate the variable
% SummaryStats, which is needed for this script. The file was too large to
% save.

%% Simulations: Plot warp matrix distance

rng('default')

figure_num = 1;
prior_color = {[1 0 0],[0 0 1],[0 1 0]}; % Red = 111, Blue = 121, Green = 131
num_reps = length(SummaryStats(1,1,1).Results);
shift_shuffle = .2;
include_shuffle = 0;

% X points for plotting
% x = [1 2 3 4]';
ind = 1;
num_col = size(SummaryStats,3);
num_warp_cond = size(SummaryStats,2);
num_prior = size(SummaryStats,1);

ax = [0 (num_col + 1) 0 .1];


% For each row (warp condition)
for j = 1:num_warp_cond
    ind = 1;
    x = []; y = []; y_EB = []; colors = [];
    subplot(3,1,j)
    
    % For each prior
    for i = 1:num_prior
        
        % For each column (each num of neurons used)
        for k = 1:num_col
            
            % For each repetition of the same conditions
            for rep = 1:num_reps
                % For the unshuffled data
                x(ind,1) = k + .01*rep;
                y(ind,1) = SummaryStats(i,j,k).Results(rep).param_recovery.tau.median_dist_tau;
                y_EB(ind,1:2) = SummaryStats(i,j,k).Results(rep).param_recovery.tau.median_dist_tau_CI';
                colors = [colors; prior_color{i}];
                ind = ind + 1;
                
                % For the shuffled data
                if include_shuffle
                    x(ind,1) = k + shift_shuffle + .01*rep;
                    y(ind,1) = SummaryStats(i,j,k).Results(rep).param_recovery.tau.median_dist_tau_shuffle;
                    y_EB(ind,1:2) = SummaryStats(i,j,k).Results(rep).param_recovery.tau.median_dist_tau_CI';
                    colors = [colors; prior_color{i}];
                    ind = ind + 1;
                end
            end
        end

    end
    
    % Plot
    scatter_better(x,y,y_EB,ax,colors)
end

fig = gcf;
fig.Units = 'normalized';
fig.Position = [.1 0 .4 1];

%% Simulations: Plot effect sizes

rng('default')

num_reps = length(SummaryStats(1,1,1).Results);
prior_color = {[1 0 0],[0 0 1],[0 1 0]}; % Red = 111, Blue = 121, Green = 131

% X points for plotting
% x = [1 2 3 4]';

ind = 1;
num_col = size(SummaryStats,3);
num_warp_cond = size(SummaryStats,2);
num_prior = size(SummaryStats,1);

ax = [.5 (num_col + .5) -.1 .1];


% For each row (warp condition)
for j = 1:num_warp_cond
    ind = 1;
    x = []; y = []; y_EB = []; colors = [];
    subplot(3,1,j)
    
    % For each prior condition
    for i = 1:num_prior
        
        % For each column (each num of neurons used)
        for k = 1:num_col
            
            % For each repetition of the same conditions
            for rep = 1:num_reps
                x(ind,1) = k + .01*rep;
                mean_effect = median(SummaryStats(i,j,k).Results(rep).RPR2(:,1));
%                 mean_effect = SummaryStats(i,j,k).Results(rep).mean_effect;
                sem = SummaryStats(i,j,k).Results(rep).mean_effect_sem;
                y(ind,1) = mean_effect;
                y_EB(ind,1:2) = mean_effect + [sem, -sem];
                colors = [colors; prior_color{i}];
                ind = ind + 1;
            end
        end
    end
    
    % Plot
    scatter_better(x,y,y_EB,ax,colors)
end

fig = gcf;
fig.Units = 'normalized';
fig.Position = [.1 0 .4 1];

%% Simulations: Plot beta (LNP) correlations

rng('default')

figure_num = 1;
prior_color = {[1 0 0],[0 0 1],[0 1 0]}; % Red = 111, Blue = 121, Green = 131
num_reps = length(SummaryStats(1,1,1).Results);
shift_shuffle = .2;
include_shuffle = 1;
include_nodtw = 1;

% X points for plotting
% x = [1 2 3 4]';
ind = 1;
num_col = size(SummaryStats,3);
num_warp_cond = size(SummaryStats,2);
num_prior = size(SummaryStats,1);

ax = [0 (num_col + 1) 0 1];


% For each row (warp condition)
for j = 1:num_warp_cond
    ind = 1;
    x = []; y = []; y_EB = []; colors = [];
    subplot(3,1,j)
    
    % For each prior
    for i = 1:num_prior
        
        % For each column (each num of neurons used)
        for k = 1:num_col
            
            % For each repetition of the same conditions
            for rep = 1:num_reps
                % For the unshuffled data
                x(ind,1) = k + .01*rep;
                y(ind,1) = SummaryStats(i,j,k).Results(rep).param_recovery.beta.median_corr_beta;
                y_EB(ind,1:2) = SummaryStats(i,j,k).Results(rep).param_recovery.beta.median_corr_beta_CI';
                colors = [colors; prior_color{i}];
                ind = ind + 1;
                
                % For the noDTW data
                if include_nodtw
                    x(ind,1) = k + shift_shuffle + .01*rep;
                    y(ind,1) = SummaryStats(i,j,k).Results(rep).param_recovery.beta.median_corr_beta_noDTW;
                    y_EB(ind,1:2) = SummaryStats(i,j,k).Results(rep).param_recovery.beta.median_corr_beta_noDTW_CI';
                    colors = [colors; prior_color{i}];
                    ind = ind + 1;
                end
                
                % For the shuffled data
                if include_shuffle
                    x(ind,1) = k + 2*shift_shuffle + .01*rep;
                    y(ind,1) = SummaryStats(i,j,k).Results(rep).param_recovery.beta.median_corr_beta_shuffle;
                    y_EB(ind,1:2) = SummaryStats(i,j,k).Results(rep).param_recovery.beta.median_corr_beta_shuffle_CI';
                    colors = [colors; prior_color{i}];
                    ind = ind + 1;
                end
            end
        end

    end
    
    % Plot
    scatter_better(x,y,y_EB,ax,colors)
end

fig = gcf;
fig.Units = 'normalized';
fig.Position = [.1 0 .4 1];

%% Real data: Plot effect sizes for all effects (main text)

% Inputs
eb_size = 0; % In SEM units
% eb_size = 1.96; % In SEM units
axes = [0 110 -.02 .06];
nrn_color = [1 0 0];
pop_color = [0 1 0];

% Load data
% Monkey 1
file_loc = 'C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Results\6-19-17\Conditions\M1\';
fname = 'EffectSizes_M1_conditions.mat';
temp = load([file_loc fname]);
EffectSizes = temp.EffectSizes;

% Monkey 2
file_loc = 'C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Results\6-19-17\Conditions\M2\';
fname = 'EffectSizes_M2_conditions.mat';
temp = load([file_loc fname]);
EffectSizes2 = temp.EffectSizes; clear temp

% Plot locations for mean/median effects
num_nrn = size(EffectSizes(1).Results.RPR2,1);
xm1 = num_nrn + 8; xm2 = xm1 + 3; xm3 = xm2 + 8; xm4 = xm3 + 3; xm5 = xm4 + 8; xm6 = xm5 + 3; xm7 = xm6 + 8; xm8 = xm7 + 3; xm9 = xm8 + 8; xm10 = xm9 + 3;

% Initialize
x = []; y = []; eb = []; colors = [];


% First scatterplot RPR2 for the individual neurons in the basic case
num_nrn = size(EffectSizes(1).Results.RPR2,1);
num_nrn2 = size(EffectSizes2(1).Results.RPR2,1);
RPR2_indiv_mean = EffectSizes(1).Results.RPR2(:,1);
RPR2_indiv_EB = EffectSizes(1).Results.RPR2(:,2:3);
jitt_base1 = linspace(-.3,.3,num_nrn)';
jitt_base2 = linspace(-.3,.3,num_nrn2)';


x_nrn = 1:num_nrn;
x = [x_nrn'];
y = [RPR2_indiv_mean];
eb = [RPR2_indiv_EB];
colors = repmat(nrn_color,[num_nrn,1]);

% ========== Next, do individual population effects
% Basic TW effect: (LNP+TW) - (LNP)
% Monkey 1
idx = 1;
pop_effect_TW_basic_mean = EffectSizes(idx).Results.mean_effect;
pop_effect_TW_basic_median = median(EffectSizes(idx).Results.RPR2(:,1));
pop_effect_TW_basic_EB = EffectSizes(idx).Results.mean_effect + EffectSizes(idx).Results.mean_effect_sem*eb_size*[-1, 1];
x = [x; xm1];
% y = [y; pop_effect_TW_basic_mean];
y = [y; pop_effect_TW_basic_median];
eb = [eb; pop_effect_TW_basic_EB];
colors = [colors; pop_color];
% Basic TW effect: jittered data points
x_jitt1 = jitt_base1 + xm1;
y_jitt1 = EffectSizes(idx).Results.RPR2(:,1);
eb_jitt1 = repmat(EffectSizes(idx).Results.RPR2(:,1),1,2); % dummy EBs
jitt_color = [.1 .1 .1];
x = [x; x_jitt1];
y = [y; y_jitt1];
eb = [eb; eb_jitt1];
colors = [colors; repmat(jitt_color,num_nrn,1)];
ns_tw1 = sum(EffectSizes(idx).Results.RPR2(:,2)>0); % num significant neurons
nv_tw1 = sum((EffectSizes(idx).Results.RPR2(:,2)>axes(3)).*(EffectSizes(idx).Results.RPR2(:,3)<axes(4))); % num visible neurons


% Monkey 2
idx = 1;
pop_effect_TW_basic_mean = EffectSizes2(idx).Results.mean_effect;
pop_effect_TW_basic_median = median(EffectSizes2(idx).Results.RPR2(:,1));
pop_effect_TW_basic_EB = EffectSizes2(idx).Results.mean_effect + EffectSizes2(idx).Results.mean_effect_sem*eb_size*[-1, 1];
x = [x; xm2];
% y = [y; pop_effect_TW_basic_mean];
y = [y; pop_effect_TW_basic_median];
eb = [eb; pop_effect_TW_basic_EB];
colors = [colors; pop_color];
% Basic TW effect: jittered data points
x_jitt1 = jitt_base2 + xm2;
y_jitt1 = EffectSizes2(idx).Results.RPR2(:,1);
eb_jitt1 = repmat(EffectSizes2(idx).Results.RPR2(:,1),1,2); % dummy EBs
jitt_color = [.1 .1 .1];
x = [x; x_jitt1];
y = [y; y_jitt1];
eb = [eb; eb_jitt1];
colors = [colors; repmat(jitt_color,num_nrn2,1)];
ns_tw2 = sum(EffectSizes2(idx).Results.RPR2(:,2)>0); % num significant neurons
nv_tw2 = sum((EffectSizes2(idx).Results.RPR2(:,2)>axes(3)).*(EffectSizes2(idx).Results.RPR2(:,3)<axes(4))); % num visible neurons


% ======== Shuffled TW: (LNP+TW_shuffled) - (LNP)
% Monkey 1
idx = 4;
pop_effect_shuffle_mean = EffectSizes(idx).Results.mean_effect;
pop_effect_shuffle_median = median(EffectSizes(idx).Results.RPR2(:,1));
pop_effect_shuffle_EB = EffectSizes(idx).Results.mean_effect + EffectSizes(idx).Results.mean_effect_sem*eb_size*[-1, 1];
x = [x; xm3];
% y = [y; pop_effect_shuffle_mean];
y = [y; pop_effect_shuffle_median];
eb = [eb; pop_effect_shuffle_EB];
colors = [colors; pop_color];
% Shuffled TW: jittered data points
x_jitt1 = jitt_base1 + xm3;
y_jitt1 = EffectSizes(idx).Results.RPR2(:,1);
eb_jitt1 = repmat(EffectSizes(idx).Results.RPR2(:,1),1,2); % dummy EBs
jitt_color = [.1 .1 .1];
x = [x; x_jitt1];
y = [y; y_jitt1];
eb = [eb; eb_jitt1];
colors = [colors; repmat(jitt_color,num_nrn,1)];
ns_tws1 = sum(EffectSizes(idx).Results.RPR2(:,2)>0); % num significant neurons
nv_tws1 = sum((EffectSizes(idx).Results.RPR2(:,2)>axes(3)).*(EffectSizes(idx).Results.RPR2(:,3)<axes(4))); % num visible neurons

% Monkey 2
idx = 4;
pop_effect_shuffle_mean = EffectSizes2(idx).Results.mean_effect;
pop_effect_shuffle_median = median(EffectSizes2(idx).Results.RPR2(:,1));
pop_effect_shuffle_EB = EffectSizes2(idx).Results.mean_effect + EffectSizes2(idx).Results.mean_effect_sem*eb_size*[-1, 1];
x = [x; xm4];
% y = [y; pop_effect_shuffle_mean];
y = [y; pop_effect_shuffle_median];
eb = [eb; pop_effect_shuffle_EB];
colors = [colors; pop_color];
% Shuffled TW: jittered data points
x_jitt1 = jitt_base2 + xm4;
y_jitt1 = EffectSizes2(idx).Results.RPR2(:,1);
eb_jitt1 = repmat(EffectSizes2(idx).Results.RPR2(:,1),1,2); % dummy EBs
jitt_color = [.1 .1 .1];
x = [x; x_jitt1];
y = [y; y_jitt1];
eb = [eb; eb_jitt1];
colors = [colors; repmat(jitt_color,num_nrn2,1)];
ns_tws2 = sum(EffectSizes2(idx).Results.RPR2(:,2)>0); % num significant neurons
nv_tws2 = sum((EffectSizes2(idx).Results.RPR2(:,2)>axes(3)).*(EffectSizes2(idx).Results.RPR2(:,3)<axes(4))); % num visible neurons


% ========= TW after spike history and kinematics: (LNP+TW+SH+Kin) - (LNP+SH+Kin)
% Monkey 1
idx = 3;
pop_effect_TW_mean = EffectSizes(idx).Results.mean_effect;
pop_effect_TW_median = median(EffectSizes(idx).Results.RPR2(:,1));
pop_effect_TW_EB = EffectSizes(idx).Results.mean_effect + EffectSizes(idx).Results.mean_effect_sem*eb_size*[-1, 1];
x = [x; xm5];
% y = [y; pop_effect_TW_mean];
y = [y; pop_effect_TW_median];
eb = [eb; pop_effect_TW_EB];
colors = [colors; pop_color];
% TW: jittered data points
x_jitt1 = jitt_base1 + xm5;
y_jitt1 = EffectSizes(idx).Results.RPR2(:,1);
eb_jitt1 = repmat(EffectSizes(idx).Results.RPR2(:,1),1,2); % dummy EBs
jitt_color = [.1 .1 .1];
x = [x; x_jitt1];
y = [y; y_jitt1];
eb = [eb; eb_jitt1];
colors = [colors; repmat(jitt_color,num_nrn,1)];
ns_twf1 = sum(EffectSizes(idx).Results.RPR2(:,2)>0); % num significant neurons
nv_twf1 = sum((EffectSizes(idx).Results.RPR2(:,2)>axes(3)).*(EffectSizes(idx).Results.RPR2(:,3)<axes(4))); % num visible neurons

% Monkey 2
idx = 3;
pop_effect_TW_mean = EffectSizes2(idx).Results.mean_effect;
pop_effect_TW_median = median(EffectSizes2(idx).Results.RPR2(:,1));
pop_effect_TW_EB = EffectSizes2(idx).Results.mean_effect + EffectSizes2(idx).Results.mean_effect_sem*eb_size*[-1, 1];
x = [x; xm6];
% y = [y; pop_effect_TW_mean];
y = [y; pop_effect_TW_median];
eb = [eb; pop_effect_TW_EB];
colors = [colors; pop_color];
% TW: jittered data points
x_jitt1 = jitt_base2 + xm6;
y_jitt1 = EffectSizes2(idx).Results.RPR2(:,1);
eb_jitt1 = repmat(EffectSizes2(idx).Results.RPR2(:,1),1,2); % dummy EBs
jitt_color = [.1 .1 .1];
x = [x; x_jitt1];
y = [y; y_jitt1];
eb = [eb; eb_jitt1];
colors = [colors; repmat(jitt_color,num_nrn2,1)];
ns_twf2 = sum(EffectSizes2(idx).Results.RPR2(:,2)>0); % num significant neurons
nv_twf2 = sum((EffectSizes2(idx).Results.RPR2(:,2)>axes(3)).*(EffectSizes2(idx).Results.RPR2(:,3)<axes(4))); % num visible neurons

% === Dividing line/space === 
% === Other effects for calibration === 

% ====== Spike history effect: (LNP+SH) - (LNP)
% Monkey 1
idx = 2;
pop_effect_SH_mean = EffectSizes(idx).Results.mean_effect;
pop_effect_SH_median = median(EffectSizes(idx).Results.RPR2(:,1));
pop_effect_SH_EB = EffectSizes(idx).Results.mean_effect + EffectSizes(idx).Results.mean_effect_sem*eb_size*[-1, 1];
x = [x; xm7];
% y = [y; pop_effect_SH_mean];
y = [y; pop_effect_SH_median];
eb = [eb; pop_effect_SH_EB];
colors = [colors; pop_color];
% SH: jittered data points
x_jitt1 = jitt_base1 + xm7;
y_jitt1 = EffectSizes(idx).Results.RPR2(:,1);
eb_jitt1 = repmat(EffectSizes(idx).Results.RPR2(:,1),1,2); % dummy EBs
jitt_color = [.1 .1 .1];
x = [x; x_jitt1];
y = [y; y_jitt1];
eb = [eb; eb_jitt1];
colors = [colors; repmat(jitt_color,num_nrn,1)];
ns_sh1 = sum(EffectSizes(idx).Results.RPR2(:,2)>0); % num significant neurons
nv_sh1 = sum((EffectSizes(idx).Results.RPR2(:,2)>axes(3)).*(EffectSizes(idx).Results.RPR2(:,3)<axes(4))); % num visible neurons

% Monkey 2
idx = 2;
pop_effect_SH_mean = EffectSizes2(idx).Results.mean_effect;
pop_effect_SH_median = median(EffectSizes2(idx).Results.RPR2(:,1));
pop_effect_SH_EB = EffectSizes2(idx).Results.mean_effect + EffectSizes2(idx).Results.mean_effect_sem*eb_size*[-1, 1];
x = [x; xm8];
% y = [y; pop_effect_SH_mean];
y = [y; pop_effect_SH_median];
eb = [eb; pop_effect_SH_EB];
colors = [colors; pop_color];
% SH: jittered data points
x_jitt1 = jitt_base2 + xm8;
y_jitt1 = EffectSizes2(idx).Results.RPR2(:,1);
eb_jitt1 = repmat(EffectSizes2(idx).Results.RPR2(:,1),1,2); % dummy EBs
jitt_color = [.1 .1 .1];
x = [x; x_jitt1];
y = [y; y_jitt1];
eb = [eb; eb_jitt1];
colors = [colors; repmat(jitt_color,num_nrn2,1)];
ns_sh2 = sum(EffectSizes2(idx).Results.RPR2(:,2)>0); % num significant neurons
nv_sh2 = sum((EffectSizes2(idx).Results.RPR2(:,2)>axes(3)).*(EffectSizes2(idx).Results.RPR2(:,3)<axes(4))); % num visible neurons

% ====== Tuning: (LNP_full) - (LNP_untuned)
% Monkey 1
idx = 5;
pop_effect_tuning_mean = EffectSizes(idx).Results.mean_effect;
pop_effect_tuning_median = median(EffectSizes(idx).Results.RPR2(:,1));
pop_effect_tuning_EB = EffectSizes(idx).Results.mean_effect + EffectSizes(idx).Results.mean_effect_sem*eb_size*[-1, 1];
x = [x; xm9];
% y = [y; pop_effect_tuning_mean];
y = [y; pop_effect_tuning_median];
eb = [eb; pop_effect_tuning_EB];
colors = [colors; pop_color];
% Tuning: jittered data points
x_jitt1 = jitt_base1 + xm9;
y_jitt1 = EffectSizes(idx).Results.RPR2(:,1);
eb_jitt1 = repmat(EffectSizes(idx).Results.RPR2(:,1),1,2); % dummy EBs
jitt_color = [.1 .1 .1];
x = [x; x_jitt1];
y = [y; y_jitt1];
eb = [eb; eb_jitt1];
colors = [colors; repmat(jitt_color,num_nrn,1)];
ns_t1 = sum(EffectSizes(idx).Results.RPR2(:,2)>0); % num significant neurons
nv_t1 = sum((EffectSizes(idx).Results.RPR2(:,2)>axes(3)).*(EffectSizes(idx).Results.RPR2(:,3)<axes(4))); % num visible neurons

% Monkey 2
idx = 5;
pop_effect_tuning_mean = EffectSizes2(idx).Results.mean_effect;
pop_effect_tuning_median = median(EffectSizes2(idx).Results.RPR2(:,1));
pop_effect_tuning_EB = EffectSizes2(idx).Results.mean_effect + EffectSizes2(idx).Results.mean_effect_sem*eb_size*[-1, 1];
x = [x; xm10];
% y = [y; pop_effect_tuning_mean];
y = [y; pop_effect_tuning_median];
eb = [eb; pop_effect_tuning_EB];
colors = [colors; pop_color];
% Tuning: jittered data points
x_jitt1 = jitt_base2 + xm10;
y_jitt1 = EffectSizes2(idx).Results.RPR2(:,1);
eb_jitt1 = repmat(EffectSizes2(idx).Results.RPR2(:,1),1,2); % dummy EBs
jitt_color = [.1 .1 .1];
x = [x; x_jitt1];
y = [y; y_jitt1];
eb = [eb; eb_jitt1];
colors = [colors; repmat(jitt_color,num_nrn2,1)];
ns_t2 = sum(EffectSizes2(idx).Results.RPR2(:,2)>0); % num significant neurons
nv_t2 = sum((EffectSizes2(idx).Results.RPR2(:,2)>axes(3)).*(EffectSizes2(idx).Results.RPR2(:,3)<axes(4))); % num visible neurons


% ========= Plot
scatter_better(x,y,eb,axes,colors)

% ========= Write text
% Number of significant neurons
txt_ht = .06;
text(xm1,txt_ht,[num2str(ns_tw1) '/' num2str(num_nrn)])
text(xm3,txt_ht,[num2str(ns_tws1) '/' num2str(num_nrn)])
text(xm5,txt_ht,[num2str(ns_twf1) '/' num2str(num_nrn)])
text(xm7,txt_ht,[num2str(ns_sh1) '/' num2str(num_nrn)])
text(xm9,txt_ht,[num2str(ns_t1) '/' num2str(num_nrn)])
text(xm2,txt_ht,[num2str(ns_tw2) '/' num2str(num_nrn2)])
text(xm4,txt_ht,[num2str(ns_tws2) '/' num2str(num_nrn2)])
text(xm6,txt_ht,[num2str(ns_twf2) '/' num2str(num_nrn2)])
text(xm8,txt_ht,[num2str(ns_sh2) '/' num2str(num_nrn2)])
text(xm10,txt_ht,[num2str(ns_t2) '/' num2str(num_nrn2)])

% Number of neurons not in field of view
txt_ht2 = .055;
text(xm1,txt_ht2,[num2str(nv_tw1) '/' num2str(num_nrn)])
text(xm3,txt_ht2,[num2str(nv_tws1) '/' num2str(num_nrn)])
text(xm5,txt_ht2,[num2str(nv_twf1) '/' num2str(num_nrn)])
text(xm7,txt_ht2,[num2str(nv_sh1) '/' num2str(num_nrn)])
text(xm9,txt_ht2,[num2str(nv_t1) '/' num2str(num_nrn)])
text(xm2,txt_ht2,[num2str(nv_tw2) '/' num2str(num_nrn2)])
text(xm4,txt_ht2,[num2str(nv_tws2) '/' num2str(num_nrn2)])
text(xm6,txt_ht2,[num2str(nv_twf2) '/' num2str(num_nrn2)])
text(xm8,txt_ht2,[num2str(nv_sh2) '/' num2str(num_nrn2)])
text(xm10,txt_ht2,[num2str(nv_t2) '/' num2str(num_nrn2)])

%% Real data: Plot effect sizes for different speed thresholds MEANS (supplement)

% Inputs
eb_size = 1.96; % In SEM units
axes = [0 5 -.01 .03];
nrn_color = [1 0 0];
pop_color = [0 1 0];

% Initialize
x = []; y = []; eb = []; colors = [];

% Load data
file_loc = 'C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Results\8-27-15\';
fname = 'EffectSizes_thresholds';
load([file_loc fname])


% Do each threshold
% Threshold = 6 cm/s
idx = 3;
effect_mean = EffectSizes(idx).Results.mean_effect;
effect_EB = EffectSizes(idx).Results.mean_effect + EffectSizes(idx).Results.mean_effect_sem*eb_size*[-1, 1];
x = [x; 1];
y = [y; effect_mean];
eb = [eb; effect_EB];
colors = [colors; pop_color];

% Threshold = 8 cm/s
idx = 4;
effect_mean = EffectSizes(idx).Results.mean_effect;
effect_EB = EffectSizes(idx).Results.mean_effect + EffectSizes(idx).Results.mean_effect_sem*eb_size*[-1, 1];
x = [x; 2];
y = [y; effect_mean];
eb = [eb; effect_EB];
colors = [colors; pop_color];

% Threshold = 10 cm/s
idx = 1;
effect_mean = EffectSizes(idx).Results.mean_effect;
effect_EB = EffectSizes(idx).Results.mean_effect + EffectSizes(idx).Results.mean_effect_sem*eb_size*[-1, 1];
x = [x; 3];
y = [y; effect_mean];
eb = [eb; effect_EB];
colors = [colors; pop_color];

% Threshold = 12 cm/s
idx = 2;
effect_mean = EffectSizes(idx).Results.mean_effect;
effect_EB = EffectSizes(idx).Results.mean_effect + EffectSizes(idx).Results.mean_effect_sem*eb_size*[-1, 1];
x = [x; 4];
y = [y; effect_mean];
eb = [eb; effect_EB];
colors = [colors; pop_color];

scatter_better(x,y,eb,axes,colors)

%% Real data: Plot effect sizes for different speed thresholds MEDIANS(supplement)

% Inputs
eb_size = 0; % In SEM units
axes = [0 5 -.01 .05];
nrn_color = [1 0 0];
pop_color = [0 1 0];
jitter_color = [0 0 0];

% Initialize
x = []; y = []; eb = []; colors = [];

% Load data
file_loc = 'C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Data\6-19-17\Results_files\Thresh\';
fname = 'EffectSizes_thresholds';
load([file_loc fname])

% Initialize
% num_nrn = size(EffectSizes(3).Results.RPR2,1);


% Do each threshold
% Threshold = 6 cm/s
idx = 3;
num_nrn = size(EffectSizes(idx).Results.RPR2,1);
jitter_base = linspace(-.05,.05,num_nrn)';
effect_mean = median(EffectSizes(idx).Results.RPR2(:,1));
effect_EB = EffectSizes(idx).Results.mean_effect + EffectSizes(idx).Results.mean_effect_sem*eb_size*[-1, 1];
x = [x; 1];
y = [y; effect_mean];
eb = [eb; effect_EB];
colors = [colors; pop_color];
% Jitter
x = [x; 1+jitter_base];
y = [y; EffectSizes(idx).Results.RPR2(:,1)];
eb = [eb; repmat(effect_EB,num_nrn,1)];
colors = [colors; repmat(jitter_color,num_nrn,1)];
ns1 = sum(EffectSizes(idx).Results.RPR2(:,2)>0);


% Threshold = 8 cm/s
idx = 4;
num_nrn = size(EffectSizes(idx).Results.RPR2,1);
jitter_base = linspace(-.05,.05,num_nrn)';
effect_mean = median(EffectSizes(idx).Results.RPR2(:,1));
effect_EB = EffectSizes(idx).Results.mean_effect + EffectSizes(idx).Results.mean_effect_sem*eb_size*[-1, 1];
x = [x; 2];
y = [y; effect_mean];
eb = [eb; effect_EB];
colors = [colors; pop_color];
% Jitter
x = [x; 2+jitter_base];
y = [y; EffectSizes(idx).Results.RPR2(:,1)];
eb = [eb; repmat(effect_EB,num_nrn,1)];
colors = [colors; repmat(jitter_color,num_nrn,1)];
ns2 = sum(EffectSizes(idx).Results.RPR2(:,2)>0);

% Threshold = 10 cm/s
idx = 1;
num_nrn = size(EffectSizes(idx).Results.RPR2,1);
jitter_base = linspace(-.05,.05,num_nrn)';
effect_mean = median(EffectSizes(idx).Results.RPR2(:,1));
effect_EB = EffectSizes(idx).Results.mean_effect + EffectSizes(idx).Results.mean_effect_sem*eb_size*[-1, 1];
x = [x; 3];
y = [y; effect_mean];
eb = [eb; effect_EB];
colors = [colors; pop_color];
% Jitter
x = [x; 3+jitter_base];
y = [y; EffectSizes(idx).Results.RPR2(:,1)];
eb = [eb; repmat(effect_EB,num_nrn,1)];
colors = [colors; repmat(jitter_color,num_nrn,1)];
ns3 = sum(EffectSizes(idx).Results.RPR2(:,2)>0);

% Threshold = 12 cm/s
idx = 2;
num_nrn = size(EffectSizes(idx).Results.RPR2,1);
jitter_base = linspace(-.05,.05,num_nrn)';
effect_mean = median(EffectSizes(idx).Results.RPR2(:,1));
effect_EB = EffectSizes(idx).Results.mean_effect + EffectSizes(idx).Results.mean_effect_sem*eb_size*[-1, 1];
x = [x; 4];
y = [y; effect_mean];
eb = [eb; effect_EB];
colors = [colors; pop_color];
% Jitter
x = [x; 4+jitter_base];
y = [y; EffectSizes(idx).Results.RPR2(:,1)];
eb = [eb; repmat(effect_EB,num_nrn,1)];
colors = [colors; repmat(jitter_color,num_nrn,1)];
ns4 = sum(EffectSizes(idx).Results.RPR2(:,2)>0);

scatter_better(x,y,eb,axes,colors)

% Add text
text(1,.03,[num2str(ns1) '/' num2str(num_nrn)])
text(2,.03,[num2str(ns2) '/' num2str(num_nrn)])
text(3,.03,[num2str(ns3) '/' num2str(num_nrn)])
text(4,.03,[num2str(ns4) '/' num2str(num_nrn)])

%% Real data: Plot effect sizes for different time bin sizes (supplement)

% Inputs
eb_size = 0; % In SEM units
axes = [0 5 -.01 .05];
nrn_color = [1 0 0];
pop_color = [0 1 0];
jitter_color = [0 0 0];

% Initialize
x = []; y = []; eb = []; colors = [];

% Load data
file_loc = 'C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Data\6-19-17\Results_files\Bin_sizes\';
fname = 'EffectSizes_bins';
load([file_loc fname])


% Do each time bin
% Bin size = 5 ms
idx = 4;
effect_mean = EffectSizes(idx).Results.median_effect;
effect_MAD = mad(EffectSizes(idx).Results.RPR2(:,1),1);
effect_EB = EffectSizes(idx).Results.mean_effect + EffectSizes(idx).Results.mean_effect_sem*eb_size*[-1, 1];
num_nrn = size(EffectSizes(idx).Results.RPR2,1);
jitter_base = linspace(-.05,.05,num_nrn)';
x = [x; 1];
y = [y; effect_mean];
eb = [eb; effect_EB];
colors = [colors; pop_color];
% Jitter
x = [x; 1+jitter_base];
y = [y; EffectSizes(idx).Results.RPR2(:,1)];
eb = [eb; repmat(effect_EB,num_nrn,1)];
colors = [colors; repmat(jitter_color,num_nrn,1)];
ns1 = sum(EffectSizes(idx).Results.RPR2(:,2)>0);


% Bin size = 10 ms
idx = 1;
effect_mean = EffectSizes(idx).Results.median_effect;
effect_MAD = mad(EffectSizes(idx).Results.RPR2(:,1),1);
effect_EB = EffectSizes(idx).Results.mean_effect + EffectSizes(idx).Results.mean_effect_sem*eb_size*[-1, 1];
num_nrn = size(EffectSizes(idx).Results.RPR2,1);
jitter_base = linspace(-.05,.05,num_nrn)';
x = [x; 2];
y = [y; effect_mean];
eb = [eb; effect_EB];
colors = [colors; pop_color];
% Jitter
x = [x; 2+jitter_base];
y = [y; EffectSizes(idx).Results.RPR2(:,1)];
eb = [eb; repmat(effect_EB,num_nrn,1)];
colors = [colors; repmat(jitter_color,num_nrn,1)];
ns1 = sum(EffectSizes(idx).Results.RPR2(:,2)>0);

% Bin size = 20 ms
idx = 2;
effect_mean = EffectSizes(idx).Results.median_effect;
effect_MAD = mad(EffectSizes(idx).Results.RPR2(:,1),1);
effect_EB = EffectSizes(idx).Results.mean_effect + EffectSizes(idx).Results.mean_effect_sem*eb_size*[-1, 1];
num_nrn = size(EffectSizes(idx).Results.RPR2,1);
jitter_base = linspace(-.05,.05,num_nrn)';
x = [x; 3];
y = [y; effect_mean];
eb = [eb; effect_EB];
colors = [colors; pop_color];
% Jitter
x = [x; 3+jitter_base];
y = [y; EffectSizes(idx).Results.RPR2(:,1)];
eb = [eb; repmat(effect_EB,num_nrn,1)];
colors = [colors; repmat(jitter_color,num_nrn,1)];
ns1 = sum(EffectSizes(idx).Results.RPR2(:,2)>0);

% Bin size = 40 ms
idx = 3;
effect_mean = EffectSizes(idx).Results.median_effect;
effect_MAD = mad(EffectSizes(idx).Results.RPR2(:,1),1);
effect_EB = EffectSizes(idx).Results.mean_effect + EffectSizes(idx).Results.mean_effect_sem*eb_size*[-1, 1];
num_nrn = size(EffectSizes(idx).Results.RPR2,1);
jitter_base = linspace(-.05,.05,num_nrn)';
x = [x; 4];
y = [y; effect_mean];
eb = [eb; effect_EB];
colors = [colors; pop_color];
% Jitter
x = [x; 4+jitter_base];
y = [y; EffectSizes(idx).Results.RPR2(:,1)];
eb = [eb; repmat(effect_EB,num_nrn,1)];
colors = [colors; repmat(jitter_color,num_nrn,1)];
ns1 = sum(EffectSizes(idx).Results.RPR2(:,2)>0);

scatter_better(x,y,eb,axes,colors)

% Add text
text(1,.03,[num2str(ns1) '/' num2str(num_nrn)])
text(2,.03,[num2str(ns2) '/' num2str(num_nrn)])
text(3,.03,[num2str(ns3) '/' num2str(num_nrn)])
text(4,.03,[num2str(ns4) '/' num2str(num_nrn)])


