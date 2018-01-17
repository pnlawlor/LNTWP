function [ acc_matrix, cost_matrix, warp_path, warp_path_matrix, WV_struct, LLHD_test, LLHD_test_diag ] = align_signals_trialwise_WV( signal_true, signal_warped, transition_priors, weights, num_WV )

    verbose = 1;

    if nargin < 5
        num_WV = 1;
    end
    
    
    % Initialize
    num_nrn = size(signal_true,2);
    num_trials = size(signal_true,1);
    acc_matrix = cell(num_trials,num_WV);
    cost_matrix = cell(num_trials,num_WV);
    warp_path = cell(num_trials,num_WV);
    warp_path_matrix = cell(num_trials,num_WV);
    LLHD_test = cell(num_trials,num_WV);
    LLHD_test_diag = cell(num_trials,num_WV);
    
    % Organize neurons for warp validation (cross validation) - sort by FR to split uniformly
    WV_struct = WV_partition(weights,num_WV);
    
    
    % Fit warp for each trial

    parfor tr = 1:num_trials
        tic
        
        for WV = 1:num_WV
        
        nrn_train = WV_struct(WV).train;
        nrn_test = WV_struct(WV).test;
            
        % For multiple neurons
        [acc_matrix{tr,WV}, cost_matrix{tr,WV}, ... 
            warp_path{tr,WV}, warp_path_matrix{tr,WV}] =  ...
                                                    align_signals_mult_nrn( ...
                                                                signal_true(tr,nrn_train), signal_warped(tr,nrn_train), transition_priors, weights);     % Grab right neurons
        
        
        % Calculate test set cost matrix (to get test set LLHD) - legacy;
        % not used anymore
        cost_matrix_test = score_function_vec2(signal_true(tr,nrn_test), signal_warped(tr,nrn_test));
        LLHD_test{tr,WV} = cost_matrix_test .* warp_path_matrix{tr,WV};
        LLHD_test_diag{tr,WV} = cost_matrix_test .* flipud(eye(size(cost_matrix_test,1)));
                                                            

        end
        
        if verbose
            disp(['Trial ' num2str(tr) ' fit in ' num2str(toc) ' seconds'])
        end
    end

end

function [ signal_out ] = downsample_conserve( signal_in, ds_factor )
    
    remainder = mod(length(signal_in),ds_factor);
    signal_in = [signal_in; nan(remainder,1)];
    num_bins_new = length(signal_in)/ds_factor;
    signal_in = reshape(signal_in,[ds_factor num_bins_new]);
    signal_out = nansum(signal_in,1)';    

end

function [ cell_out ] = downsample_conserve_cellwise( cell_in, ds_factor )

    cell_out = cellfun(@downsample_conserve,cell_in, repmat({ds_factor},size(cell_in)),'UniformOutput',false);

end

function [WV_struct] = WV_partition( weights, num_WV )

    if num_WV > 1 % If warp validation desired
        [~,idx_sort] = sort(weights,'descend');
        num_nrn = length(weights);

        for wv = 1:num_WV
            test_nrn_temp = wv:num_WV:num_nrn; % Set up train/test alternation
            test_nrn = idx_sort(test_nrn_temp); % Alternate within the sorted neurons
            WV_struct(wv).test = logical(zeros(num_nrn,1));
            WV_struct(wv).test(test_nrn) = 1;
            WV_struct(wv).train = ~(WV_struct(wv).test);
        end
    elseif num_WV == 1 % If no warp validation desired
        num_nrn = length(weights);
        WV_struct(1).test = logical(ones(num_nrn,1));
        WV_struct(1).train = logical(ones(num_nrn,1));
    end

end

function [ score ] = score_function_vec2 (FR, spikes)

    nn = size(FR,2);

    FR = cell2mat(FR);
    FR = flipud(FR);
    spikes = cell2mat(spikes);
    
    FR_rep = repmat(FR,[1 1 size(spikes,1)]);
    FR_rep = permute(FR_rep,[1 3 2]);

    spikes_rep = repmat(spikes,[1 1 size(FR,1)]);
    spikes_rep = permute(spikes_rep,[3 1 2]);

    score = spikes_rep.*log(FR_rep) - FR_rep - log(factorial(spikes_rep));
    
    score = sum(score,3);

end