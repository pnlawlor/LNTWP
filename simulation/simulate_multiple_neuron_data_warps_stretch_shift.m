% Simulate single neuron data for time warping project

%% Load scripts
rng('default')

addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Scripts\Pat'))
addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Projects\FEF\Scripts\Pavan_GLM_code_v2\FEF\utilities'))
addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Scripts_TW'))
addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Scripts\Pat'))

%% User inputs

num_nrn = 1;
num_trials = 200;
bin_size = 10; % is ms
bins_per_trial = 61;
dt = bin_size/1000;
FR_baseline = 2;

baseline_term_mean = log(FR_baseline*dt); 

tuning_strength_mean = log(1 + 3);
untuned_strength_mean = log(1 + 1);

% tuning_strength_mean = log(1 + 5);
% untuned_strength_mean = log(1 + 3);

measured_onset = 41; % When the movement started (bin#)

transition_prob = [ 0 1; 0 1; 0 1 ];
run_len = 1;

transition_prob = transition_prob / sum(sum(transition_prob)); % Normalize

visualize_data = 0;

fname_base = 'SimData';
fpath_save = 'C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Simulation_data\';

rng('default')

%% Choose model to simulate from

model_num = 8;
filt_struct = TW_define_filters(0,model_num,dt); 
num_filt = length(filt_struct);

% For model num 8, there are 14 filters, all generating untuned, cos(), sin()

%% Initialize

firing_rates = repmat({zeros(bins_per_trial,1)},num_trials,num_nrn);
spikes = firing_rates;
movement_data = repmat({zeros(bins_per_trial,2)},num_trials,1);
X_cell = cell(num_trials,1);
X_cell_warped = cell(num_trials,1);
warp_matrix = cell(num_trials,1);
FR_all = nan(num_trials*bins_per_trial,num_nrn);
FR_all_prewarp = nan(num_trials*bins_per_trial,num_nrn);
FR_cell = cell(num_trials,num_nrn);
FR_prewarp_cell = cell(num_trials,num_nrn);
spikes_all = nan(num_trials*bins_per_trial,num_nrn);
spikes = cell(num_trials,num_nrn);
spikes_all_prewarp = nan(num_trials*bins_per_trial,num_nrn);
spikes_prewarp = cell(num_trials,num_nrn);

predictions = cell(num_nrn,1);
fit_parameters = cell(num_nrn,1);
fit_info = cell(num_nrn,1);
pseudo_R2 = cell(num_nrn,1);

%% Randomize parameters of simulated cells

% Choose basleine term for each cell
baseline_term = normrnd(baseline_term_mean,.1*abs(baseline_term_mean), num_nrn,1);

% Choose preferred direction of each cell
dir_pref = 2*pi*rand(num_nrn,1);
% dir_pref = pi*ones(num_nrn,1);

% Choose tuning strength of each cell
tuning_strength = normrnd(tuning_strength_mean,.2*abs(tuning_strength_mean), num_nrn,1);
untuned_strength = normrnd(untuned_strength_mean,.2*abs(untuned_strength_mean), num_nrn,1);

% Choose temporal parameters for each cell
untuned_temporal_param_temp = normrnd(0,.15,[num_nrn, num_filt]);
tuned_temporal_param_temp = normrnd(0,.15,[num_nrn, num_filt]);

% Multiply them all together
untuned_param_all = ...
    repmat(untuned_strength,1,num_filt) .* ...
    untuned_temporal_param_temp;
    
cos_param_all = ...
    repmat(tuning_strength,1,num_filt) .* ...
    repmat(cos(dir_pref),1,num_filt) .* ...
    tuned_temporal_param_temp;

sin_param_all = ...
    repmat(tuning_strength,1,num_filt) .* ...
    repmat(sin(dir_pref),1,num_filt) .* ...
    tuned_temporal_param_temp;

% Form weight vector
beta = [baseline_term, untuned_param_all, cos_param_all, sin_param_all];

% Simple model
% beta(:,2) = 1;
% beta(:,3:end) = 0;

%% Randomize parameters of simulated behavioral data

% Choose movement direction
movement_angle = 2*pi*rand(num_trials,1);

% movement_angle = pi*ones(num_trials,1);
% movement_angle(2:2:end) = zeros(num_trials/2,1);

%% Generate un-warped movement data

% Generate data for stereotypical trial
temp = zeros(bins_per_trial,3);
temp(measured_onset,1) = 1;

movement_data_prewarp = repmat({temp},num_trials,1);

for tr = 1:num_trials
    % Add cos() and sin() components for each trial
    movement_data_prewarp{tr}(measured_onset,2) = cos(movement_angle(tr));
    movement_data_prewarp{tr}(measured_onset,3) = sin(movement_angle(tr));
    
    % Filter (old way)
%     X_cell{tr} = filter_and_insert(movement_data_prewarp{tr},filt_struct);
end

%% Filter movement data to create covariates

% Make duplicate with 0s for padding
temp = [movement_data_prewarp cellfun(@(x)x*0,movement_data_prewarp,'UniformOutput',0)];
% Reshape so that every other trial is a blank
temp2 = cell2mat(reshape(temp',[],1));
% Filter
X_temp = full(filter_and_insert(temp2,filt_struct));
% Remove padding
bins_per_trial_rep = repmat(bins_per_trial,2*num_trials,1);
X_cell_temp = mat2cell(X_temp,bins_per_trial_rep,size(X_temp,2));
X_cell(:,1,1) = X_cell_temp(1:2:end);
% Clear temp variable
clear X_cell_temp

%% Warp each trial

% Generate warp
for tr = 1:num_trials
    warp_matrix{tr,1} = generate_warp([bins_per_trial, bins_per_trial], transition_prob, run_len);
end

% Apply warp to covariates
X_cell_warped = apply_warp_cov_trialwise_WV( ...
                        warp_matrix, X_cell, ...
                        1:size(X_cell{1},2) ...
                            );
                        
%% Warp each trial - choose stretch/shift
% Strategy: simulate warps; reject those that don't meet criteria. E.g., if
% you want stretches only, reject warps that are shifts and continue
% simulating until you get the number of desired trials.
% NB: This is only for one cell at a time
% Details: User specifies proportion of trials desired to be "shifts"

% User inputs
prop_shift = 1; % Proportion of trials desired to be shifts
crit_low = 1; % Threshold below which relative FR change is denoted a shift
crit_high = 1; % Threshold above which relative FR change is denoted a stretch
warp_min = 0; % Minimum warp distance required to denote a shift

% Initialize
num_sim_trials = 0;
tr = 1;
num_shift = 0;
num_stretch = 0;
nrn = 1;
FR_rel_array = nan;
X_prewarp = stdize(cell2mat(X_cell));

% One trial at a time
while tr <= num_trials
    
    % Generate a warp
    warp_matrix{tr,1} = generate_warp([bins_per_trial, bins_per_trial], transition_prob, run_len);
    
    % Apply warp
    X_cell_warped = apply_warp_cov_trialwise_WV( ...
                        warp_matrix(tr,1), X_cell(tr,1), ...
                        1:size(X_cell{1},2) ...
                            );
                        
    % Get FR for this trial
    X_all = stdize(cell2mat(X_cell_warped));
    FR = glmval(beta(nrn,:)',X_all,'log');
    FR_prewarp = glmval(beta(nrn,:)',X_prewarp,'log');
    FR_rel = abs(mean(FR) - mean(FR_prewarp)) / mean(FR_prewarp); % relative change in firing rate due to warp
    FR_rel_array = [FR_rel_array; FR_rel];
    
    disp(['Relative FR: ' num2str(FR_rel)])
    
    % Check criterion
    is_shift = (FR_rel < crit_low);
    is_stretch = (FR_rel > crit_high);
    
    % If meets criterion, increment and move on to next trial
    if is_shift
        % Check to make sure there aren't too many shifts
        prop_shift_new = (num_shift + 1)/num_trials;
        
        % If not too many, keep this warp
        if prop_shift_new <= prop_shift
            tr = tr + 1;
            num_shift = num_shift + 1;
            disp(['Shift accepted for trial ' num2str(tr-1)])
        end
        
    elseif is_stretch
        % Check to make sure there aren't too many stretches (1-shifts)
        prop_stretch_new = (num_stretch + 1)/num_trials;
        
        % If not too many, keep this warp
        if prop_stretch_new <= (1 - prop_shift)
            tr = tr + 1;
            num_stretch = num_stretch + 1;
            disp(['Stretch accepted for trial ' num2str(tr-1)])
        end
    end
        
    
end

% Apply all warps to all trials
X_cell_warped = apply_warp_cov_trialwise_WV( ...
                        warp_matrix, X_cell, ...
                        1:size(X_cell{1},2) ...
                            );

%% Simulate neural data based upon warped covariates

% Standardize covariates
X_all = stdize(cell2mat(X_cell_warped));
X_prewarp = stdize(cell2mat(X_cell));

for nrn = 1:num_nrn
    FR_all(:,nrn) = glmval(beta(nrn,:)',X_all,'log');
    FR_cell(:,nrn) = mat2cell(FR_all(:,nrn),repmat(bins_per_trial,num_trials,1),1);
    
    FR_all_prewarp(:,nrn) = glmval(beta(nrn,:)',X_prewarp,'log');
    FR_prewarp_cell(:,nrn) = mat2cell(FR_all_prewarp(:,nrn),repmat(bins_per_trial,num_trials,1),1);
    
    spikes_all(:,nrn) = poissrnd(FR_all(:,nrn));
    spikes(:,nrn) = mat2cell(spikes_all(:,nrn),repmat(bins_per_trial,num_trials,1),1);
    
    spikes_all_prewarp(:,nrn) = poissrnd(FR_all_prewarp(:,nrn));
    spikes_prewarp(:,nrn) = mat2cell(spikes_all_prewarp(:,nrn),repmat(bins_per_trial,num_trials,1),1);
end

%% Plot gain vs TW strength

TW_strength = nan(num_trials,1);
FR_rel_plot = nan(num_trials,1);

for tr = 1:num_trials
    TW_strength(tr) = compare_warp_matrices(warp_matrix{tr,1},flipud(eye(size(warp_matrix{tr,1},1))));
    
%     FR1 = mean(FR_cell{tr,nrn});
    FR1 = mean(spikes{tr,nrn});
    FR0 = mean(FR_prewarp_cell{tr,nrn});
    FR_rel_plot(tr) = abs(FR1 - FR0)/FR0;
end    

figure,
scatter(TW_strength,FR_rel_plot)
xlim([0 .75])
ylim([0 .75])

%% Plot spikes and FR

wind = [1:100];
nrns = [];

for nrn_idx = 1:length(nrns)
    nrn = nrns(nrn_idx);
    figure
    plot(spikes_all(wind,nrn)); hold on
    plot(FR_all(wind,nrn),'r'); hold on
    plot(FR_all_prewarp(wind,nrn),'g'); hold off
    
    title([ ...
        'Nrn: ' num2str(nrn) ...
        ', MeanFR: ' num2str(mean(spikes_all(:,nrn))) ...
        ', PD: ' num2str(dir_pref(nrn)) ...
        ', Tuning str: ' num2str(tuning_strength(nrn)) ...
        ', Untuned str: ' num2str(untuned_strength(nrn))])
    
    ylim([0 1.2])
end

%% Plot rasters

trs = [];

for tr_idx = 1:length(trs)
    tr = trs(tr_idx);
    
    spikes_to_plot = cell2mat(spikes(tr,:))';
    FR_to_plot = cell2mat(FR_cell(tr,:))';
    FR_prewarp_to_plot = cell2mat(FR_prewarp_cell(tr,:))';
    
    figure
    subplot(3,1,1)
    imagesc(spikes_to_plot,[0 2])
    subplot(3,1,2)
    imagesc(FR_to_plot, [0, 1])
    subplot(3,1,3)
    imagesc(FR_prewarp_to_plot, [0, 1])
    
end

%% Fit each neuron to make sure it works

nrn_to_fit = 1;
num_CV = 2;
dt = 1;
lambda = [.005];
alpha = .01 ;
fit_method = 'lassoglm';
bins_per_trial = 61;
num_BS_GLM = 1;

cov_X = cov(X_all);
inv_cov_X = inv(cov_X);
% inv_cov_X = eye(60);

for nrn_idx = 1:length(nrn_to_fit)

    nrn = nrn_to_fit(nrn_idx);
    
    [predictions{nrn,1}, ...
                fit_parameters{nrn,1}, ...
                fit_info{nrn,1}, ...
                pseudo_R2{nrn,1}] = ...
                                        fit_poiss_GLM( X_all, spikes_all(:,nrn), ...
                                                num_CV, ...
                                                dt, ...
                                                lambda, ... % lambda
                                                alpha, ... % alpha
                                                fit_method, ...
                                                bins_per_trial, ...
                                                num_BS_GLM);
                                            
    % Calculate agreement between simulated and inferred parameters
    ang = calc_vector_distance(beta(nrn,2:end)',fit_parameters{nrn}{1}(2:end));
    disp(['Neuron: ' num2str(nrn) '. Angle between simulated and inferred parameters: ' num2str(ang*180/pi)])
    
    rand_ind = randperm(length(fit_parameters{nrn}{1}(2:end))) + 1;
    ang1b = calc_vector_distance(beta(nrn,2:end)',fit_parameters{nrn}{1}(rand_ind));
    disp(['Neuron: ' num2str(nrn) '. Angle between simulated and shuffled parameters: ' num2str(ang1b*180/pi)])
    
    % Calculate performance metrics
    ang2 = calc_mahalanobis_dist(beta(nrn,2:end)',fit_parameters{nrn}{1}(2:end),inv_cov_X);
    disp(['Neuron: ' num2str(nrn) '. Mahalanobis distance between simulated and inferred parameters: ' num2str(ang2)])
    ang3 = calc_mahalanobis_dist(fit_parameters{nrn}{2}(2:end),fit_parameters{nrn}{1}(2:end),inv_cov_X);
    disp(['Neuron: ' num2str(nrn) '. Mahalanobis distance between cross validations: ' num2str(ang3)])
    ang4 = calc_mahalanobis_dist(beta(nrn,2:end)',fit_parameters{nrn}{1}(rand_ind),inv_cov_X);
    disp(['Neuron: ' num2str(nrn) '. Mahalanobis distance between simulated and shuffled parameters: ' num2str(ang4)])
    
    % Plot parameters
    param_plot = nan*repmat(fit_parameters{nrn}{2}(2:end),1,(num_CV + 1));
    param_plot_shuffle = nan*repmat(fit_parameters{nrn}{2}(2:end),1,(num_CV + 1));
    param_plot(:,1) = beta(nrn,2:end)';
    param_plot_shuffle(:,1) = beta(nrn,2:end)';
    
    for cv = 1:num_CV
        param_plot(:,cv+1) = fit_parameters{nrn}{cv}(2:end);
        param_plot_shuffle(:,cv+1) = fit_parameters{nrn}{cv}(rand_ind);
    end
    
%     figure('units','normalized','outerposition',[.25 .25 .5 .5])
%     subplot(1,2,1)
%     imagesc(param_plot)
%     title(['Correlation: ' num2str(corr(param_plot(:,1),param_plot(:,2)))])
%     subplot(1,2,2)
%     imagesc(param_plot_shuffle)
%     title(['Correlation: ' num2str(corr(param_plot_shuffle(:,1),param_plot_shuffle(:,2)))])

    % Show rasters
    [~,idx_sort] = sort(movement_angle);

    figure('units','normalized','outerposition',[.25 .25 .5 .5])
    subplot(1,2,1)
    imagesc(reshape(cell2mat(spikes(idx_sort,nrn)),61,[])')
    title(['FR: ' num2str(mean(spikes_all(:,nrn))) ', PR2: ' num2str(pseudo_R2{nrn,1}(1))])
    subplot(1,2,2)
    imagesc(reshape(cell2mat(FR_cell(idx_sort,nrn)),61,[])')
    title(['Corr: ' num2str(corr(param_plot(:,1),param_plot(:,2))) ', ShuffCorr: ' num2str(corr(param_plot_shuffle(:,1),param_plot_shuffle(:,2)))])
    
end

%% Prepare data for saving

% Meta-information
description.num_nrn = num_nrn;
description.num_trials = num_trials;
description.bin_size = bin_size;
description.bins_per_trial = bins_per_trial;
description.FR_baseline = FR_baseline;
description.baseline_term_mean = baseline_term_mean;
description.tuning_strength_mean = tuning_strength_mean;
description.tuning_strength = tuning_strength;
description.untuned_strength_mean = untuned_strength_mean;
description.untuned_strength = untuned_strength;
description.transition_prob = transition_prob;

% GLM parameters
parameters.beta = beta;

% Write to struct
Data_win.description = description;
Data_win.parameters = parameters;
Data_win.warp_matrices = warp_matrix;
Data_win.spikes_PMd = spikes;
Data_win.FR_postwarp = FR_cell;
Data_win.FR_prewarp = FR_prewarp_cell;
Data_win.X_cell_prewarp = X_cell;
Data_win.X_cell_postwarp = X_cell_warped;
Data_win.movement_cov = movement_data_prewarp;
Data_win.movement_angle = movement_angle;
Data_win.target_dir = movement_angle;
Data_win.dir_pref = dir_pref;

%% Save data

datetimestr = datestr(now,'mm-dd-yy--HH-MM');
fname = [fname_base '_' num2str(num_nrn) 'nrn_' num2str(num_trials) 'tr-' datetimestr];

save([fpath_save fname],'Data_win')

disp(['Data saved: ' fname])
