function [ X_cell_warped ] = apply_warp_cov_trialwise_WV( warp_path_matrices, X_cell, col_to_warp )

    verbose = 0;

    % Initialize
    num_WV = size(X_cell,2);
    num_trials = size(X_cell,1);
    X_cell_warped = X_cell;
    
    % Apply warp for each trial
    for tr = 1:num_trials
        for WV = 1:num_WV
            for col_idx = 1:length(col_to_warp)
                col = col_to_warp(col_idx);

                X_cell_warped{tr,WV}(:,col) = ...
                                                apply_warp_cov(warp_path_matrices{tr,WV}, X_cell{tr,WV}(:,col));

            end
        end
    end

    disp(['Covariates warped.'])

end

