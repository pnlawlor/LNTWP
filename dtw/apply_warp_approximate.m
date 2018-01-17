function [ warped_signal ] = apply_warp_approximate( warp_path_matrix, spikes, FR )

    warped_signal = zeros(size(warp_path_matrix,1),1);
    
    num_rows = size(warp_path_matrix,1);
    num_cols = size(warp_path_matrix,2);
    
    for row = 1:num_rows
        if sum(warp_path_matrix(row,:)) > 1
            idx = warp_path_matrix(row,:) > 0
        end
    end
    
    for col = 1:num_cols
        if sum(warp_path_matrix(:,col)) > 1
            
        end
    end   
    
end

