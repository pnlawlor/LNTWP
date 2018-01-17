function [ warped_signal ] = apply_warp_cov( warp_path_matrix, signal )

    warp_path_matrix(warp_path_matrix==0) = nan;

    signal_rep = repmat(signal,[1 size(warp_path_matrix,2)]);
    
    warped_signal = nanmean(signal_rep.*flipud(warp_path_matrix),1)';
    
end

