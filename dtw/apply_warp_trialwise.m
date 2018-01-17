function [ signals_rewarped ] = apply_warp_trialwise( warp_path_matrices, spikes, FR )

    verbose = 0;

    % Initialize
    num_trials = size(spikes,1);
    num_nrn = size(spikes,2);
    signals_rewarped = cell(num_trials,num_nrn);
    
    % Fit warp for each trial
    for tr = 1:num_trials
        for nrn_num = 1:num_nrn
            tic
            signals_rewarped{tr,nrn_num} =  ...
                                            apply_warp(warp_path_matrices{tr}, spikes{tr,nrn_num});
                                        
            signals_rewarped{tr,nrn_num} = ...
                                            apply_warp_approximate(warp_path_matrices{tr}, spikes{tr,nrn_num}, FR{tr,nrn_num});
            
            if verbose
                disp(['Neuron ' num2str(nrn_num) ',trial ' num2str(tr) ' un-warped in ' num2str(toc) ' seconds'])
            end
        end
    end


end

