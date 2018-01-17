% Generate the model covariates. I have custom scripts
% that form the model covariates from the experimental data. It generally
% involves filtering the experimental data with various types of smooth
% functions. My ecosystem allows the user to choose which experimental
% variables to filter, which functions to filter with, and other parameters. 
% The "define_filters" functions generate a structure with this information
% which is later used by the function that actually does the filtering,
% filter_and_insert.m

filt_struct = TW_define_filters(0,model_num, dt); % For external covariates
filt_struct2 = TW_define_filters2( num_nrn , num_spk_history_bf ,0, dt); % For spike history covariates

% Bookkeeping: Get number of bins per trial
bins_per_trial = cellfun(@size,movement_data,'UniformOutput',0);
bins_per_trial = cell2mat(bins_per_trial);
bins_per_trial = bins_per_trial(:,1); % for the non-padded cell array
bins_per_trial2 = reshape([bins_per_trial bins_per_trial]',[],1); % for the padded cell array

% Filter neural data for spike history terms
if include_spk_history
    % Add padding between trials to prevent spillover
    temp = cell2mat(spikes);
    spikes_to_filter = mat2cell(temp,bins_per_trial,num_nrn);
    spikes_to_filter = [spikes_to_filter cellfun(@(x)x*0,spikes_to_filter,'UniformOutput',0)]; % pad
    spikes_to_filter = cell2mat(reshape(spikes_to_filter',[],1)); % reshape
    
    % Filter
    spk_hist_cov = full(filter_and_insert( spikes_to_filter, filt_struct2 ));
    spk_hist_cov = mat2cell(spk_hist_cov,bins_per_trial2,num_nrn*num_spk_history_bf);
    spk_hist_cov = spk_hist_cov(1:2:end);
    
    % Trim
    spk_hist_cov = trim_edges(spk_hist_cov,trim_times,dt,2,movement_data);
    
    spk_hist_cov = cell2mat(spk_hist_cov);
    

end

% Filter external covariates
if model_num == 8 || model_num == 9
        % Make duplicate with 0s for padding (every other "trial" is
        % blank). This prevents filter spillover into adjacent trials.
        temp = [movement_data cellfun(@(x)x*0,movement_data,'UniformOutput',0)];
        % Reshape
        temp2 = cell2mat(reshape(temp',[],1));
        % Filter
        X_temp = full(filter_and_insert(temp2,filt_struct));
        % Remove padding
        X_cell_temp = mat2cell(X_temp,bins_per_trial2,size(X_temp,2));
        X_cell_temp2 = X_cell_temp(1:2:end);
        % Trim edges
        X_cell_temp2 = trim_edges(X_cell_temp2,trim_times,dt,2,movement_data);
        
        X_cell(:,1,1) = X_cell_temp2;
        
        clear X_cell_temp
    elseif model_num == 10
        % Make a cell for speed (not included in original kinematic data)
        kin_all = cell2mat(kinematics_data(trials));
        speed_temp = sqrt(kin_all(:,3).^2 + kin_all(:,4).^2);
        speed_cell = mat2cell(speed_temp,bins_per_trial,1);
        
        % Gather all covariates to be filtered into a single cell array
        cov_to_filter = [cell2mat(movement_data) cell2mat(kinematics_data(trials)) cell2mat(speed_cell(trials))];
        cov_to_filter(:,4:5) = []; % Get rid of x,y-position
        cov_to_filter2 = mat2cell(cov_to_filter,bins_per_trial,size(cov_to_filter,2));

        % Make duplicate with 0s for padding (every other "trial" is
        % blank). This prevents filter spillover into adjacent trials.
        temp = [cov_to_filter2 cellfun(@(x)x*0,cov_to_filter2,'UniformOutput',0)];
        % Reshape
        temp2 = cell2mat(reshape(temp',[],1));
        % Filter
        X_temp = full(filter_and_insert(temp2,filt_struct));
        % Remove padding
        X_cell_temp = mat2cell(X_temp,bins_per_trial2,size(X_temp,2));
        X_cell_temp2 = X_cell_temp(1:2:end);
        % Trim
        X_cell_temp2 = trim_edges(X_cell_temp2,trim_times,dt,2,movement_data);
        
        X_cell(:,1,1) = X_cell_temp2;
        
        clear X_cell_temp
    elseif model_num == 11
        % Make a cell for speed (not included in original kinematic data)
        kin_all = cell2mat(kinematics_data(trials));
        speed_temp = sqrt(kin_all(:,3).^2 + kin_all(:,4).^2);
        speed_cell = mat2cell(speed_temp,bins_per_trial,1);
        
        % Gather all covariates to be filtered into a single cell array
        cov_to_filter = [cell2mat(movement_data(trials)) cell2mat(kinematics_data(trials)) cell2mat(speed_cell(trials))];
        cov_to_filter(:,4:5) = []; % Get rid of x,y-position
        cov_to_filter2 = mat2cell(cov_to_filter,bins_per_trial,size(cov_to_filter,2));

        % Make duplicate with 0s for padding (every other "trial" is
        % blank). This prevents filter spillover into adjacent trials.
        temp = [cov_to_filter2 cellfun(@(x)x*0,cov_to_filter2,'UniformOutput',0)];
        % Reshape
        temp2 = cell2mat(reshape(temp',[],1));
        % Filter
        X_temp = full(filter_and_insert(temp2,filt_struct));
        % Remove padding
        X_cell_temp = mat2cell(X_temp,bins_per_trial2,size(X_temp,2));
        X_cell_temp2 = X_cell_temp(1:2:end);
        % Trim
        X_cell_temp2 = trim_edges(X_cell_temp2,trim_times,dt,2,movement_data);
        
        X_cell(:,1,1) = X_cell_temp2;        
        
        clear X_cell_temp
end


% Trim edges for spikes
spikes2 = spikes;
spikes = trim_edges(spikes,trim_times,dt,2,movement_data);
spikes_warped(:,:,1) = spikes;

% After trimming, get new number of bins per trial
bins_per_trial = cellfun(@size,X_cell,'UniformOutput',0);
bins_per_trial = cell2mat(bins_per_trial);
bins_per_trial = bins_per_trial(:,1);

% Other pre-processing
X_cell(:,2:num_WV,1) = repmat(X_cell(:,1,1),[1 (num_WV-1)]); % Duplicate starting X matrix across WV_folds

% Set up some spike history stuff 
num_exp_cov = size(cell2mat(X_cell(:,1,1)),2);
for nrn = 1:num_nrn
    idx_spk_hist_cov{nrn} = num_spk_history_bf*(nrn - 1) + [1:num_spk_history_bf] + num_exp_cov;
end


