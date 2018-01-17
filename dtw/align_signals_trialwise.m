function [ acc_matrix, cost_matrix, warp_path, warp_path_matrix ] = align_signals_trialwise( signal_true, signal_warped, transition_priors, weights, num_WV )

    verbose = 1;

    if nargin < 5
        num_WV = 1;
    end
    
    
    % Initialize
    num_nrn = size(signal_true,2);
    num_trials = size(signal_true,1);
    acc_matrix = cell(num_trials,1);
    cost_matrix = cell(num_trials,1);
    warp_path = cell(num_trials,1);
    warp_path_matrix = cell(num_trials,1);
    
    % Downsample signals
%     k = 10;
%     tic
%     signal_true_ds = downsample_conserve_cellwise( signal_true, k);
%     signal_warped_ds = downsample_conserve_cellwise( signal_warped, k);
%     disp(['Downsampling done in ' num2str(toc) ' seconds'])

    % Fit warp for each trial
    parfor tr = 1:num_trials
        tic
        
        % For single neuron
%         [acc_matrix{tr}, cost_matrix{tr}, warp_path{tr}, warp_path_matrix{tr}] =  ...
%                             align_signals(signal_true{tr,:}, signal_warped{tr,:}, transition_priors);
        
        % For multiple neurons
        [acc_matrix{tr}, cost_matrix{tr}, ... 
            warp_path{tr}, warp_path_matrix{tr}] =  ...
                                                    align_signals_mult_nrn( ...
                                                                signal_true(tr,:), signal_warped(tr,:), transition_priors, weights);      
        
        
%         [acc_matrix{tr}, cost_matrix{tr}, ... 
%             warp_path{tr}, warp_path_matrix{tr}] =  ...
%                                                     align_signals_mult_nrn( ...
%                                                                 signal_true_ds(tr,:), signal_warped_ds(tr,:), transition_priors, weights);
        
        
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