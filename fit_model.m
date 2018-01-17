% Initialize termination condition (convergence or max iterations)
converged = 0;

% Prior to alternation, fit data with vanilla GLM
X_all = cell2mat(X_cell(:,1,1));

% Include spike history terms if desired
if include_spk_history
    X_all = [X_all spk_hist_cov];
end

% Standardize covariates
X_all = stdize(X_all);

% Convert y from cell to matrix
y_all = cell2mat(spikes);

% Which covariates to "warp" - i.e., the experimental ones, not SH
X_col_to_warp = 1:num_exp_cov;

num_BS_GLM = 1; % LEGACY: Number of times to bootstrap the Pseudo-R2 during fitting. Default is 1. 

%% Fit initial GLM

iter = 1; % Initialize at iter 1

parfor nrn_idx = 1:length(neurons)
    nrn_num = neurons(nrn_idx);
    
    % Define which covariates to include
    cov_to_include = [1:num_exp_cov];
    
    % Fit GLM without DTW to get baseline
    disp(['Now fitting neuron: ' num2str(nrn_num)])
    
    [predictions_combined{nrn_idx,1}, ...
        fit_parameters{nrn_idx,1}, ...
        fit_info{nrn_idx,1}, ...
        pseudo_R2{nrn_idx,1}] = ...
                                fit_poiss_GLM( X_all(:,cov_to_include), y_all(:,nrn_num), ...
                                        num_CV, ...
                                        dt, ...
                                        lambda, ... % lambda
                                        alpha, ... % alpha
                                        fit_method, ...
                                        bins_per_trial, ...
                                        num_BS_GLM);

    % Fit GLM with spike history but without DTW
    if include_spk_history
        disp(['Now fitting neuron: ' num2str(nrn_num) ' with spike history.'])
        
        cov_to_include = [1:num_exp_cov idx_spk_hist_cov{nrn_num}]; % Define which covariates to include
        
        [predictions_combined_SH{nrn_idx,1}, ...
            fit_parameters_SH{nrn_idx,1}, ...
            fit_info_SH{nrn_idx,1}, ...
            pseudo_R2_SH{nrn_idx,1}] = ...
                                    fit_poiss_GLM( X_all(:,cov_to_include), y_all(:,nrn_num), ...
                                            num_CV, ...
                                            dt, ...
                                            lambda, ... % lambda
                                            alpha, ... % alpha
                                            fit_method, ...
                                            bins_per_trial, ...
                                            num_BS_GLM);
    end

end

% Transfer some stuff
if ~include_spk_history
    predictions_combined_SH = predictions_combined;
end

for nrn_idx = 1:length(neurons)
    % Break up predictions for all trials into trialwise predictions
    predictions(:,nrn_idx,1) = mat2cell(nansum(cell2mat(predictions_combined{nrn_idx,1}),1)', ...
        bins_per_trial,1);
    predictions_SH(:,nrn_idx,1) = mat2cell(nansum(cell2mat(predictions_combined_SH{nrn_idx,1}),1)', ...
        bins_per_trial,1);
    predictions2(:,nrn_idx,1) = mat2cell(nansum(cell2mat(predictions_combined{nrn_idx,1}),1)', ...
        bins_per_trial,1);
end

% Calculate LLHD for initial GLM
opt.type = 'poisson';
[LLHD_total(1), LLHD_mean(1)] = calc_LLHD(y_all,predictions_combined(:,1),opt);

%% Now do alternation: 1) fit DTW, 2) fit GLM

while ~converged && (iter < max_iter)
    % Update iteration counter
    iter = iter + 1;
    
    disp(['Beginning iteration: ' num2str(iter)])
        

        weights = mean(y_all,1); % Avg firing rates across neurons
        
        % Fit DTW
        [acc_matrix(:,:,iter), cost_matrix(:,:,iter), ...
            warp_path(:,:,iter), warp_path_matrix(:,:,iter), ...
            WV_folds{iter}, ...
            LLHD_test(:,:,iter), LLHD_test_diag(:,:,iter)] = ...
                                                            align_signals_trialwise_WV( ...
                                                                    predictions2(:,:,iter-1), spikes_warped(:,:,1), transition_priors, weights, num_WV); % Fit new warp to data from last iteration

                                                                
        % Apply warp to original covariates
        X_cell(:,:,iter) = apply_warp_cov_trialwise_WV( ...
                                warp_path_matrix(:,:,iter), X_cell(:,:,1), X_col_to_warp);

    
    y_all = cell2mat(spikes);
    
    % Fit GLM to all neurons by WV fold
    for WV = 1:num_WV
        disp(['Iteration: ' num2str(iter) ', WV: ' num2str(WV)])
        wv_test_nrns = find(WV_folds{iter}(WV).test);
        
        X_all = cell2mat(X_cell(:,WV,iter));
        X2 = cell2mat(X_cell(:,1,1));
        
        % Standardize covariates for regularization
        X_all = stdize(X_all);
        X2 = stdize(X2);
        
        % Fit GLM (test set neurons for given warp validation)
        parfor nrn_idx = 1:length(wv_test_nrns)
            nrn_num = wv_test_nrns(nrn_idx);
            y = y_all(:,nrn_num);
            
            % Choose exp covariates only
            cov_to_use = [1:num_exp_cov];            

            % Fit GLM
            disp(['Now fitting neuron: ' num2str(wv_test_nrns(nrn_idx))])
            [predictions_combined_temp{nrn_idx,1}, ... % These indexes are fine; just have to use 1 instead of iter for parfor reasons
                fit_parameters_temp{nrn_idx,1}, ...
                fit_info_temp{nrn_idx,1}, ...
                pseudo_R2_temp{nrn_idx,1}] = ...
                                            fit_poiss_GLM( X_all(:,cov_to_use), y, ...
                                                    num_CV, ...
                                                    dt, ... 
                                                    lambda, ... % lambda
                                                    alpha, ... % alpha
                                                    fit_method, ...
                                                    bins_per_trial, ...
                                                    num_BS_GLM);
        end
        
        % Transfer some stuff
        for nrn_idx = 1:length(wv_test_nrns)
            nrn_num = wv_test_nrns(nrn_idx);
            cov_to_use = [1:num_exp_cov];
            
            % Make predictions using WEIGHTS from current GLM fit but NOT
            % covariates from current GLM fit
            predictions_combined2_temp = predictions_combined_temp;
            for cv = 1:num_CV
                idx = ~isnan(predictions_combined_temp{nrn_idx}{cv}); % Find data points in test set
                predictions_combined2_temp{nrn_idx,1}{cv}(idx) = glmval( ... % Predict on test data
                            fit_parameters_temp{nrn_idx,1}{cv}, ...
                            X2(idx,cov_to_use), ...
                            'log', ...
                            'constant', 'on', ...
                            'offset',repmat(log(dt),size(X2(idx,1))));
            end
            
            % Replace everything as needed
            predictions_combined{nrn_num,iter} = predictions_combined_temp{nrn_idx,1};
            predictions_combined2{nrn_num,iter} = predictions_combined2_temp{nrn_idx,1};
            fit_parameters{nrn_num,iter} = fit_parameters_temp{nrn_idx};
            fit_info{nrn_num,iter} = fit_info_temp{nrn_idx};
            pseudo_R2{nrn_num,iter} = pseudo_R2_temp{nrn_idx};
            
            % Break up predictions for all trials into trialwise predictions
            predictions(:,nrn_num,iter) = mat2cell(nansum(cell2mat(predictions_combined{nrn_num,iter}),1)', ...
                bins_per_trial,1);
            predictions2(:,nrn_num,iter) = mat2cell(nansum(cell2mat(predictions_combined2{nrn_num,iter}),1)', ...
                bins_per_trial,1);
        end
    end
    
    % Calculate LLHD
    [LLHD_total(iter), LLHD_mean(iter)] = calc_LLHD(y_all,predictions_combined(:,iter),opt);

    
    % Determine if fit has coverged
    % Choose tolerance calculation type
    LL0 = LLHD_mean(iter-1);
    LL1 = LLHD_mean(iter);
    
    switch tol_type
        case 'relTol'
            delta = abs(LL1 - LL0) / min(abs(LL0), abs(LL1));
            converged = delta < conv_thresh;
        case 'absTol'
            delta = abs(LL1 - LL0);
            converged = delta < conv_thresh;
        case 'maxIter'
            delta = 0;
            converged = (iter >= max_iter);
    end
    
    disp(['Model convergence: ' num2str(converged) '. Delta: ' num2str(delta) '. Using ' tol_type])
    
end

%% If including spike history

if include_spk_history
    
    % Increment iter so as not to overwrite last one
    iter = iter + 1;
    
    % Keep DTW info from last iter (iter-1). Not fitting DTW again.
    acc_matrix(:,:,iter) = acc_matrix(:,:,iter-1);
    cost_matrix(:,:,iter) = cost_matrix(:,:,iter-1);
    warp_path(:,:,iter) = warp_path(:,:,iter-1);
    warp_path_matrix(:,:,iter) = warp_path_matrix(:,:,iter-1);
    WV_folds{iter} = WV_folds{iter-1};
    LLHD_test(:,:,iter) = LLHD_test(:,:,iter-1);
    LLHD_test_diag(:,:,iter) = LLHD_test_diag(:,:,iter-1);
    X_cell(:,:,iter) = X_cell(:,:,iter-1);
    
    % Fit GLM to all neurons by WV fold
    for WV = 1:num_WV
        disp(['Iteration: ' num2str(iter) ', WV: ' num2str(WV)])
        wv_test_nrns = find(WV_folds{iter}(WV).test);
        
        % Form covariate matrix
        X_all = cell2mat(X_cell(:,WV,iter));
        
        % Include spike history terms
        X_all = [X_all spk_hist_cov];
        
        % Fit GLM (test set neurons for given warp validation)
        parfor nrn_idx = 1:length(wv_test_nrns)
            nrn_num = wv_test_nrns(nrn_idx);
            y = y_all(:,nrn_num);
            
            % Overhead for spike history terms
            cov_to_use = [1:num_exp_cov, idx_spk_hist_cov{nrn_num}];

            % Fit GLM
            disp(['Now fitting neuron: ' num2str(wv_test_nrns(nrn_idx))])
            [predictions_combined_temp{nrn_idx,1}, ... % These indexes are fine; just have to use 1 instead of iter for parfor reasons
                fit_parameters_temp{nrn_idx,1}, ...
                fit_info_temp{nrn_idx,1}, ...
                pseudo_R2_temp{nrn_idx,1}] = ...
                                            fit_poiss_GLM( X_all(:,cov_to_use), y, ...
                                                    num_CV, ...
                                                    dt, ... 
                                                    lambda, ... % lambda
                                                    alpha, ... % alpha
                                                    fit_method, ...
                                                    bins_per_trial, ...
                                                    num_BS_GLM);
        end
        
        % Transfer some stuff
        for nrn_idx = 1:length(wv_test_nrns)
            nrn_num = wv_test_nrns(nrn_idx);
            
            cov_to_use = [1:num_exp_cov, idx_spk_hist_cov{nrn_num}];
            
            predictions_combined2_temp = predictions_combined_temp;
            
            % Replace everything as needed
            predictions_combined{nrn_num,iter} = predictions_combined_temp{nrn_idx,1};
            predictions_combined2{nrn_num,iter} = predictions_combined2_temp{nrn_idx,1};
            fit_parameters{nrn_num,iter} = fit_parameters_temp{nrn_idx};
            fit_info{nrn_num,iter} = fit_info_temp{nrn_idx};
            pseudo_R2{nrn_num,iter} = pseudo_R2_temp{nrn_idx};
            
            % Break up predictions for all trials into trialwise predictions
            predictions(:,nrn_num,iter) = mat2cell(nansum(cell2mat(predictions_combined{nrn_num,iter}),1)', ...
                bins_per_trial,1);
            predictions2(:,nrn_num,iter) = mat2cell(nansum(cell2mat(predictions_combined2{nrn_num,iter}),1)', ...
                bins_per_trial,1);
        end
    end

% Calculate LLHD
[LLHD_total(iter), LLHD_mean(iter)] = calc_LLHD(y_all,predictions_combined(:,iter),opt);

end

%%  If doing shuffle control
% Shuffle trials for shuffle control. Will apply warp to "wrong trials"
% with the same reach direction
if do_shuffle
    
    % Increment iter so as not to overwrite last one
    iter = iter + 1;
    
    % Keep DTW info from last iter (iter-1). Not fitting DTW again.
    acc_matrix(:,:,iter) = acc_matrix(:,:,iter-1);
    cost_matrix(:,:,iter) = cost_matrix(:,:,iter-1);
    warp_path(:,:,iter) = warp_path(:,:,iter-1);
    warp_path_matrix(:,:,iter) = warp_path_matrix(:,:,iter-1);
    warp_path_matrix_shuffle(:,:,iter) = warp_path_matrix(:,:,iter-1);
    WV_folds{iter} = WV_folds{iter-1};
    LLHD_test(:,:,iter) = LLHD_test(:,:,iter-1);
    LLHD_test_diag(:,:,iter) = LLHD_test_diag(:,:,iter-1);
    X_cell(:,:,iter) = X_cell(:,:,iter-1);
    
    X_shuffle = X_cell;

    % Find all of the binned reach directions (e.g., 1:8, -8:8)
    possible_reach_dirs = unique(cell2mat(Data_win.target_dir_num));

    % For each of these, shuffle the covariate matrices
    for dir_reach_idx = 1:length(possible_reach_dirs);
        dir_reach = possible_reach_dirs(dir_reach_idx);

        trials_to_shuffle = find(cell2mat(Data_win.target_dir_num)==dir_reach); % Find the trials with this reach direction
        idx_shuffle = trials_to_shuffle(randperm(length(trials_to_shuffle))); % Mix up the trial numbers
        
        warp_path_matrix_shuffle(trials_to_shuffle,:,iter) = warp_path_matrix_shuffle(idx_shuffle,:,iter);
    end
    
    % Apply shuffled warp matrices to trials 
    X_shuffle(:,:,iter) = apply_warp_cov_trialwise_WV( ...
                                warp_path_matrix_shuffle(:,:,iter), X_cell(:,:,1), X_col_to_warp);
    
    
    % Fit GLM to all neurons by WV fold
    for WV = 1:num_WV
        disp(['Iteration: ' num2str(iter) ', WV: ' num2str(WV)])
        wv_test_nrns = find(WV_folds{iter}(WV).test);
        
        X_all = cell2mat(X_shuffle(:,WV,iter));
        
        % Fit GLM (test set neurons for given warp validation)
        parfor nrn_idx = 1:length(wv_test_nrns)
            nrn_num = wv_test_nrns(nrn_idx);
            y = y_all(:,nrn_num);
            
            % Overhead for spike history terms
            cov_to_use = [1:num_exp_cov];

            % Fit GLM
            disp(['Now fitting neuron: ' num2str(wv_test_nrns(nrn_idx))])
            [predictions_combined_temp{nrn_idx,1}, ... % These indexes are fine; just have to use 1 instead of iter for parfor reasons
                fit_parameters_temp{nrn_idx,1}, ...
                fit_info_temp{nrn_idx,1}, ...
                pseudo_R2_temp{nrn_idx,1}] = ...
                                            fit_poiss_GLM( X_all(:,cov_to_use), y, ...
                                                    num_CV, ...
                                                    dt, ... 
                                                    lambda, ... % lambda
                                                    alpha, ... % alpha
                                                    fit_method, ...
                                                    bins_per_trial, ...
                                                    num_BS_GLM);
        end
        
        % Transfer some stuff
        for nrn_idx = 1:length(wv_test_nrns)
            nrn_num = wv_test_nrns(nrn_idx);
            
            % If including spike history terms
            if include_spk_history
                cov_to_use = [1:num_exp_cov, idx_spk_hist_cov{nrn_num}];
            end
            
            predictions_combined2_temp = predictions_combined_temp;
            
            % Replace everything as needed
            predictions_combined{nrn_num,iter} = predictions_combined_temp{nrn_idx,1};
            predictions_combined2{nrn_num,iter} = predictions_combined2_temp{nrn_idx,1};
            fit_parameters{nrn_num,iter} = fit_parameters_temp{nrn_idx};
            fit_info{nrn_num,iter} = fit_info_temp{nrn_idx};
            pseudo_R2{nrn_num,iter} = pseudo_R2_temp{nrn_idx};
            
            % Break up predictions for all trials into trialwise predictions
            predictions(:,nrn_num,iter) = mat2cell(nansum(cell2mat(predictions_combined{nrn_num,iter}),1)', ...
                bins_per_trial,1);
            predictions2(:,nrn_num,iter) = mat2cell(nansum(cell2mat(predictions_combined2{nrn_num,iter}),1)', ...
                bins_per_trial,1);
        end
    end
end

% Calculate LLHD
[LLHD_total(iter), LLHD_mean(iter)] = calc_LLHD(y_all,predictions_combined(:,iter),opt);

pr2 = cell2mat(pseudo_R2(:,1));
[~,idx_sort] = sort(weights,'descend');

%% Bookkeeping for simulation repetitions

if ~real_data
    % Save simulated (ground truth) and inferred DTW parameters
    Results(1).warp_matrices_sim = cellfun(@sparse,Data_win.warp_matrices,'UniformOutput',false); % Save as sparse matrix
    Results(rep).warp_matrices_inf = cellfun(@sparse,warp_path_matrix,'UniformOutput',false); % Save as sparse matrix
    
    % Save simulated (ground truth) GLM parameters
    Results(1).beta_sim = Data_win.parameters.beta(nrn_to_use_orig_idx,:); % Get the sim parameters from the *actually chosen* simulated neurons
   
    % Collect inferred GLM parameters
    beta_inferred = nan*Data_win.parameters.beta(nrn_to_use_orig_idx,:); % initialize
    beta_inferred = repmat(beta_inferred,[1 1 num_CV]); % make a copy for each cv fold
    beta_inferred_noDTW = beta_inferred;
    
    for nrn_idx = 1:length(neurons)
        nrn_num = neurons(nrn_idx);
        
        for cv = 1:num_CV
            % Model with DTW
            fit_param_temp = fit_parameters{nrn_num,iter}{cv}; % stored as column (bad design, i know)
            fit_param_temp = fit_param_temp'; % convert to row
            beta_inferred(nrn_num,:,cv) = fit_param_temp;
            
            % Model without DTW
            fit_param_temp2 = fit_parameters{nrn_num,1}{cv}; % stored as column (bad design, i know)
            fit_param_temp2 = fit_param_temp2'; % convert to row
            beta_inferred_noDTW(nrn_num,:,cv) = fit_param_temp2;
        end
    end
    
    Results(rep).beta_inf = beta_inferred;
    Results(rep).beta_inf_noDTW = beta_inferred_noDTW;
    
end

% Save results for analyses
description.num_nrn = num_nrn;
description.num_trials = num_trials;
description.sim_file = data_fname;

parameters.transition_priors = transition_priors;
parameters.transition_prior_scale = transition_prior_scale;
parameters.fit_method = fit_method; 
parameters.num_CV = num_CV;
parameters.num_WV = num_WV; 
parameters.dt = dt;
parameters.num_iter = iter;
parameters.max_iter = max_iter;
parameters.tol_type = tol_type;
parameters.conv_thresh = conv_thresh; 
parameters.delta = delta;
parameters.LLHD_mean = LLHD_mean;
parameters.alpha = alpha;
parameters.lambda = lambda;
parameters.do_shuffle = do_shuffle; 
parameters.include_spk_history = include_spk_history;
parameters.num_spk_history_bf = num_spk_history_bf;

Results(1).description = description;
Results(1).parameters = parameters;

% Save results
if ~quest
    save([Results_fpath Results_fname],'Results','-v7.3')
else
    save(Results_fname,'Results','-v7.3')
end

disp(['Results file saved: ' Results_fname])