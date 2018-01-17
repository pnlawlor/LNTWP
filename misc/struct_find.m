function [ locations ] = struct_find( struct, item, current_path )
% Recursive search of a struct for whatever might be at the end of that
% struct. Or the rainbow. You never know.

if nargin < 3
    current_path = [];
end

num_results = 1;
locations = {};
f = fields(struct);
current_path = [current_path, 1];

% For all fields in a struct
for j = 1:length(f);

    current_path(end) = j;

    % Check to see if that field is also a struct
    if isstruct(struct.(f{j}))
        % If so, keep drilling
        temp = struct_find(struct.(f{j}), item, current_path);
        if ~isempty(temp)
            locations(num_results:(num_results+length(temp)-1)) = temp; 
            num_results = num_results + 1;
        end
    % Do this later
    elseif iscell(struct.(f{j}))
        continue
    else
        % If not, check to see if it's the thing you're looking for
        if struct.(f{j}) == item % do this recursively?
            locations{num_results} = current_path;
            num_results = num_results + 1;
        else
            continue 
        end
    end


end

