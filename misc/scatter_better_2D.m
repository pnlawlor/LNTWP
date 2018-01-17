function [ ] = scatter_better_2D( x, EB_x, y, EB_y, ax, color )

% Initialization
num_dp = size(x,1);

% figure

% Draw EB lines first
for dp = 1:num_dp
    % Vertical line
    line([x(dp) x(dp)], [EB_y(dp,1) EB_y(dp,2)],'Color','Black'); hold on
    % Horizontal line
    line([EB_x(dp,1) EB_x(dp,2)], [y(dp) y(dp)],'Color','Black'); hold on
end

% Scatterplot the centers
% scatter(x,y,'filled','MarkerFaceColor',color); hold off
scatter(x,y,[],color,'o'); hold on

% Optional: add identity line
min_val = min(ax(1),ax(3));
max_val = max(ax(2),ax(4));
line([min_val max_val],[min_val max_val]); hold off

axis(ax)
axis square

end
