function [ new_X, filt_cov_struct, new_cov_names ] = filter_and_insert( old_X, filt_struct, cov_struct )
%filter_and_insert filters columns in the X matrix and inserts the new,
%filtered version

verbose = 0;

if nargin < 3
    do_cov_names = 0;
else
    do_cov_names = 1;
end

if verbose
    disp('Filtering ...')
end

tic

%% Initialization stuff

if do_cov_names
cov_names = cov_struct.Cov_names;
end

num_filters = length(filt_struct);
cov_num_in = cell(num_filters); filt_func = cell(num_filters); filt_shift = cell(num_filters);
new_cov_names = cell(1,1);

for filter_num = 1:num_filters
    cov_num_in{filter_num} = filt_struct(filter_num).cov_num_in;
    filt_func{filter_num} = filt_struct(filter_num).filt_func;
    filt_shift{filter_num} = filt_struct(filter_num).filt_shift;
    filt_names{filter_num} = filt_struct(filter_num).filt_name;
end


% Keep track of original indices
orig_size = size(old_X,2);
new_idx = 1:orig_size;
new_cov_names = cell(1,orig_size);

%% Filter

% For each filter, do the deed
for filt_num = 1:num_filters
    % For each original (unfiltered) covariate to filter
    for cov_num = 1:length(cov_num_in{filt_num})
        % Find the right covariate to filter
        read_idx = find(new_idx == cov_num_in{filt_num}(cov_num));
        if filt_shift{filt_num} < 0
            % Tack on zeros at the end
            to_filter = [old_X(:,read_idx); zeros(abs(filt_shift{filt_num}),1)];
%                 to_filter = old_X(:,read_idx);
            % Filter
            filt_temp = filt_func{filt_num};
            temp = filter(filt_temp,1,full(to_filter));
%             temp = filtfilt(filt_temp,1,full(to_filter));
            % Shift back the desired amount
            temp = temp((abs(filt_shift{filt_num})+1):end);
            % Insert into matrix
            old_X = [old_X(:,1:read_idx), temp, old_X(:,(read_idx+1):end)];
        else
            to_filter = old_X(:,read_idx);
            filt_temp = [zeros(abs(filt_shift{filt_num}),1); filt_func{filt_num}];
            % Filter
            temp = filter(filt_temp,1,full(to_filter));
%             temp = filtfilt(filt_temp,1,full(to_filter));
            old_X = [old_X(:,1:read_idx), temp, old_X(:,(read_idx+1):end)];
        end
        % Increment variable to account for filtered cov just inserted
        new_idx = [new_idx(1:read_idx), inf, new_idx((read_idx+1):end)];
        if do_cov_names
        temp_new_name = [cov_names{cov_num_in{filt_num}(cov_num)} '_' filt_names{filt_num}];
        new_cov_names = horzcat(new_cov_names(1:read_idx), {temp_new_name}, new_cov_names((read_idx+1):end));
%         new_cov_names{read_idx + 1} = [cov_names{cov_num_in{filt_num}(cov_num)} '_' filt_names{filt_num}];
        end
    end
end

% Account for this stuff (new indices) in the cov_struct
% This is kind of stupid and complicated but it works
if do_cov_names
filt_cov_struct = cov_struct;
else
filt_cov_struct = [];
end

for orig_idx = 1:orig_size
    idx1 = find(new_idx==orig_idx);
    
    if orig_idx == orig_size; 
        idx2 = length(new_idx)+1; 
    else
        idx2 = find(new_idx==(orig_idx+1));
    end
    
    num_inf = idx2 - idx1 - 1;
    filt_idx = ((1:num_inf) + idx1 - new_idx(idx1))*1i;
    if do_cov_names
    struct_path = struct_find(cov_struct.Covariates, orig_idx);
    filt_cov_struct.Covariates = struct_set_value(filt_cov_struct.Covariates,struct_path{1},filt_idx*-1i);
    end
end

    
good_cols = isinf(new_idx);
new_X = sparse(old_X(:,good_cols));
if do_cov_names
new_cov_names = new_cov_names(good_cols);
filt_cov_struct.Cov_names = new_cov_names;
end

if verbose
disp(['Filtered in ' num2str(toc) ' seconds'])
end

end
