function [ X_cell_warped ] = apply_warp_cov_trialwise( warp_path_matrices, X_cell, col_to_warp )

    verbose = 0;

    % Initialize
    num_trials = size(X_cell,1);
    X_cell_warped = X_cell;
    
    % Fit warp for each trial
    for tr = 1:num_trials
        for col_idx = 1:length(col_to_warp)
            col = col_to_warp(col_idx);
                                        
            X_cell_warped{tr}(:,col) = ...
                                            apply_warp_cov(warp_path_matrices{tr}, X_cell{tr}(:,col));
            
        end
    end

    disp(['Covariates warped.'])

end

