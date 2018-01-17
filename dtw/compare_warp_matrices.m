function [ distance_norm, distance_matrix ] = compare_warp_matrices( matrix1, matrix2 )
% Calculate the distance between two warps.

distance_matrix = abs(cumsum(matrix1,1) - cumsum(matrix2,1));
distance_norm = sum(sum(distance_matrix)) / (size(matrix1,1)*size(matrix1,2));


end

