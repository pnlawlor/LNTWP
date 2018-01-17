% Simulate single neuron data for time warping project

%% Load scripts
rng('default')

addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Scripts\Pat'))
addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Projects\FEF\Scripts\Pavan_GLM_code_v2\FEF\utilities'))
addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Scripts_TW'))
addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Scripts\Pat'))

%% User inputs

num_nrn = 250;
num_trials = 400;
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
run_len = 3;

transition_prob = transition_prob / sum(sum(transition_prob)); % Normalize

visualize_data = 0;

fname_base = 'SimData';
fpath_save = 'C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Simulation_data\';

rng('default')

%% Choose model to simulate from

model_num = 8;
filt_struct = TW_define_filters(0,model_num); 
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

%% Filter movement data to create covariates (new way)

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

%% Plot spikes and FR

wind = [1:1000];
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
