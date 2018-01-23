% Open correct parallel pool
if quest
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool
    end
else
    poolobj = gcp('nocreate');
    if isempty(poolobj)
%         parpool % Uncomment this to use parallel pool on local machine.
    end
end

% More user input
trim_times = [0.4, 0.2]; % The size of the time window to use, in seconds. Default [0.4, 0.2] which is 400ms before and 200ms after reach start.

% Check if save folder(s) exist. If not, create them
if ~(exist(Results_fpath,'dir') == 7)
    mkdir(Results_fpath);
end

% File name for saving
Results_fname = [Results_fname_base '_' datetimestr];

% Select only desired trials
movement_data = movement_data_all(trials);
spikes_all_preselect = spikes_all_preselect(trials,:);
if real_data
    kinematics_data = kinematics(trials);
end

% Get rid of low FR neurons
spikes_all_preselect_trimmed = trim_edges(spikes_all_preselect,trim_times,dt,2,movement_data); % This is just to get the neurons with the right FRs
mean_FR = mean(cell2mat(spikes_all_preselect_trimmed));
good_nrn = mean_FR > FR_thresh;
spikes_all = spikes_all_preselect(:,good_nrn); % NB: This uses spikes for the good neurons (non-low FR neurons), but keeps the full-length data. The full-length data is important for spike history terms. The spike data will be trimmed to the right length later.

% Choose the right neurons
if size(spikes_all,2) > num_nrn_desired % If you desire fewer neurons than provided
    nrn_to_use = randsample(sum(good_nrn),num_nrn_desired);

    % Keep track of original indices of neurons (mostly for simulation purposes)
    good_nrn_orig_idx = find(good_nrn);
    nrn_to_use_orig_idx = good_nrn_orig_idx(nrn_to_use);

    spikes_all = spikes_all(:,nrn_to_use);
    neurons = 1:size(spikes_all,2);

elseif size(spikes_all,2) == num_nrn_desired % If you've requested exactly the right number of neurons
    neurons = 1:size(spikes_all,2);

elseif size(spikes_all,2) < num_nrn_desired % If you've requested too many neurons
    neurons = 1:size(spikes_all,2);
    disp(['You have requested more neurons than you have available in this dataset. Using ' num2str(length(neurons)) ' instead.'])

end

num_nrn_max = size(spikes_all,2);

disp(['Using ' num2str(size(spikes_all,2)) '/' num2str(size(Data_win.spikes_PMd,2)) ' neurons, FR threshold: ' num2str(round(FR_thresh/dt)) ' Hz'])

% Add extra iterations for model with spike history terms. I just store
% this as an extra iteration of the alternation.
max_iter = max_iter + include_spk_history;

if include_spk_history
    disp('Including spike history covariates. This adds an extra iteration.')
end

spikes = spikes_all(:,neurons);

num_trials = length(trials);
num_nrn = length(neurons);

warps = cell(num_trials,max_iter);
spikes_warped = cell(num_trials,num_nrn,max_iter);
spikes_warped(:,:,1) = spikes;
X_cell = cell(num_trials,num_WV,max_iter);
acc_matrix = cell(num_trials,num_WV, max_iter);
cost_matrix = cell(num_trials,num_WV,max_iter);
LLHD_test = cell(num_trials,num_WV,max_iter);
LLHD_test_diag = cell(num_trials,num_WV,max_iter);
LLHD_total = nan(max_iter,1);
LLHD_mean = nan(max_iter,1);
warp_path = cell(num_trials,num_WV,max_iter);
warp_path_matrix = cell(num_trials,num_WV,max_iter);
WV_folds = cell(max_iter,1);
idx_spk_hist_cov = cell(num_nrn,1);

predictions_combined = cell(num_nrn,max_iter);
predictions = cell(num_trials,num_nrn,max_iter);
predictions_combined_SH = cell(num_nrn,1);
predictions_SH = cell(num_trials,num_nrn,1);
predictions_combined2 = cell(num_nrn,max_iter);
predictions2 = cell(num_trials,num_nrn,max_iter);
predictions_temp = cell(num_trials,num_nrn,1);
fit_parameters = cell(num_nrn,max_iter);
fit_parameters_SH = cell(num_nrn,1);
fit_parameters_temp = cell(num_nrn,1);
fit_info = cell(num_nrn,max_iter);
fit_info_SH = cell(num_nrn,1);
fit_info_temp = cell(num_nrn,1);
pseudo_R2 = cell(num_nrn,max_iter);
pseudo_R2_SH = cell(num_nrn,1);
pseudo_R2_temp = cell(num_nrn,1);
