function [ X_cell_out ] = trim_edges( X_cell_in, trim_times, dt, mode, movement_data )
%UNTITLED Trims temporal window
%   To avoid edge effects of filtering, I filter a larger time window than
%   ultimately desired. After filtering, this function is used to get rid
%   of the extra time added to the time window. More specifically, this
%   function takes X_cell_in (a cell array of design matrices), and trims
%   trim_times(1) time from the beginning, and trim_times(2) from the end. It
%   then returns X_cell_out. 
%  
% Assumes that the data of interest is in the first column of X_cell_in

if nargin < 4
    mode = 1;
end

% Initialize
X_cell_out = cell(size(X_cell_in)); % initialize
num_trials = size(X_cell_in,1); % get num trials
num_col = size(X_cell_in,2); % get num columns

if mode == 1
    
if sum(trim_times) == 0
    disp(['Not trimming... continuing'])
    X_cell_out = X_cell_in;
else
    disp(['Cutting ' num2str(trim_times(1)) ' sec from beginning, ' num2str(trim_times(2)) ' sec from end.'])
end

% Mode 1: trim fixed length from each edge
num_bins_to_cut1 = round(trim_times(1)/dt); % get #bins from time
num_bins_to_cut2 = round(trim_times(2)/dt);

for tr_num = 1:num_trials
    for col_num = 1:num_col
        temp = X_cell_in{tr_num,col_num};
        temp(1:num_bins_to_cut1,:) = []; % cut beginning bins
        temp(end:-1:(end - num_bins_to_cut2 + 1),:) = []; % cut end bins
        X_cell_out{tr_num,col_num} = temp;
    end
end

elseif mode == 2
% Mode 2: trim to desired length

disp(['Trimming to ' num2str(trim_times(1)) ' s before, and ' num2str(trim_times(2)) ' s after.'])

for tr_num = 1:num_trials
    
    num_bins_trial = size(X_cell_in{tr_num},1);
    center_bin = find(movement_data{tr_num}(:,1));
    pre_bins = round(trim_times(1)/dt);
    post_bins = round(trim_times(2)/dt);
    num_bins_to_cut1 = center_bin - 1 - pre_bins;
    num_bins_to_cut2 = (num_bins_trial - center_bin) - post_bins;
    
    for col_num = 1:num_col
        temp = X_cell_in{tr_num,col_num};
        temp(1:num_bins_to_cut1,:) = []; % cut beginning bins
        temp(end:-1:(end - num_bins_to_cut2 + 1),:) = []; % cut end bins
        X_cell_out{tr_num,col_num} = temp;
    end
end

end

end

