function [ neural_data_out, bad_nrn_out ] = remove_duplicate_neurons( neural_data_in, corr_threshold, FR_thresh )
    
    if nargin < 3
        FR_thresh = 0;
        verbose = 1;
    end
    
    if nargin < 4
        verbose = 1;
    end

    if iscell(neural_data_in)
        neural_data = cell2mat(neural_data_in);
        cellmode = 1;
    else
        neural_data = neural_data_in;
        cellmode = 0;
    end

    % Calculate correlation
    corr_structure = abs(corr(neural_data));
    mean_FR = mean(neural_data);
    
    % Find nrn below FR threshold
    bad_nrn_FR = find(mean_FR<FR_thresh);
    
    % Initialize
    bad_pairs = corr_structure > corr_threshold;
    bad_pairs(logical(eye(size(bad_pairs)))) = 0;
    bad_nrn = [];
    
    % For each bad pair, throw out the neuron with lower FR
    bad_pairs2 = find(bad_pairs);
    
    for bp = 1:length(bad_pairs2)
        bp_num = bad_pairs2(bp);
        
        [bp_sub1, bp_sub2] = ind2sub(size(bad_pairs),bp_num);
        
        if mean_FR(bp_sub1) > mean_FR(bp_sub2)
            bad_nrn = [bad_nrn bp_sub2];
        else
            bad_nrn = [bad_nrn bp_sub1];
        end
    end
    
    bad_nrn_all = unique([bad_nrn bad_nrn_FR]);
    
    % Another format (Boolean)
    bad_nrn_out = zeros(size(neural_data_in,2),1);
    bad_nrn_out(bad_nrn_all) = 1;
    
    neural_data_out = neural_data_in;
    neural_data_out(:,bad_nrn_all) = [];
    
    if verbose
        disp(['Using ' num2str(sum(~bad_nrn_out)) '/' num2str(length(bad_nrn_out)) ' neurons.'])
    end
    
end

