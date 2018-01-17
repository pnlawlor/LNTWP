function [ output_args ] = plot_signals_and_matrix( sig1, sig2, matrix, path_matrix, caxis )

    if nargin < 5
        caxis = nan;
        plot_path = 0;
        path_weight = 0;
        path_matrix = matrix*0;
    else
        plot_path = 1;
        path_weight = prctile(matrix(:),50);% + .1*max(matrix(:));
    end

    figure
    subplot(4,4,[1 5 9])
    plot(sig1,1:length(sig1))
    ylim([1 length(sig1)])
    
    subplot(4,4,[2 3 4 6 7 8 10 11 12])
    
    if ~isnan(caxis)
        imagesc(matrix + plot_path*path_weight*path_matrix,caxis)
    else
        imagesc(matrix + plot_path*path_weight*path_matrix) 
    end
    
    subplot(4,4,[14 15 16])
    plot(1:length(sig2),sig2)
    xlim([1 length(sig2)])

end

