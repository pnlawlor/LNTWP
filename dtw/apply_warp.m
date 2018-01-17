function [ warped_signal ] = apply_warp( warp_path_matrix, signal )

    warped_signal = zeros(size(warp_path_matrix,1),1);
    
    for row = length(warped_signal):-1:1
        idx_warped_signal = length(warped_signal) - row + 1;
        
        idx_this_row = warp_path_matrix(row,:) > 0;
        
%         warped_signal(idx_warped_signal) = sum(signal(idx_this_row)); % Take the average value for the new signal; round b/c we're doing Poisson
        warped_signal(idx_warped_signal) = round(mean(signal(idx_this_row))); % Take the average value for the new signal; round b/c we're doing Poisson
%         warped_signal(idx_warped_signal) = mean(signal(idx_this_row)); % Take the average value for the new signal; no round for linear
    end

    
%% Approximate unwarping
    
end

