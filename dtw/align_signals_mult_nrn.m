function [ acc_cost_matrix, cost_matrix, path, path_matrix ] = align_signals_mult_nrn( FR, spikes, transition_priors, weights )

%% Initialization

if nargin < 3
    transition_priors = [1, 1; 0, 1];
end

num_nrn = size(FR,2);
FR_cell = FR;
spikes_cell = spikes;
FR = cell2mat(FR);
spikes = cell2mat(spikes);


%% Fill out cost matrix

% This is the point at which the cost matrices are combined across neurons.
% This simply adds their cost matrices and thus does not assume a weighted
% linear combination which could be the most general case

FR = flipud(FR);

cost_matrix = score_function_vec(FR,spikes);
cost_matrix = sum(cost_matrix,3);

%% Fill out cost matrix

% acc_cost_matrix = dynamic_function( cost_matrix, transition_priors);

weight_type = 0; % Type of weighting across neurons. Default is 0.

if weight_type == 0  

    [acc_cost_matrix, pointer_matrix] = dynamic_function2(cost_matrix, transition_priors);
    
elseif weight_type == 1 % Weight by avg FR
    
    alpha = 1; % 1 = Vanilla weighting. 0 = no weighting. 

    weights_rep = repmat(weights,[size(FR,1) 1 size(spikes,1)]);
    weights_rep = permute(weights_rep,[1 3 2]);
    weights_rep = log(weights_rep + eps);
    
    mm = cost_matrix + alpha*weights_rep;
    
    acc_cost_matrix = dynamic_function2(sum(mm,3), transition_priors);
    
    cost_matrix(:,:,end+1) = sum(mm,3);
    
elseif weight_type == 2 % Weight by avg FR and trial FR
    
    alpha = 1; % 1 = Vanilla weighting. 0 = no weighting. 
    beta = .1; % 1 = entirely avg FR. 0 = entirely within-trial FR
    
    % Across trial FR weights
    weights_rep = repmat(weights,[size(FR,1) 1 size(spikes,1)]);
    weights_rep = permute(weights_rep,[1 3 2]); % Reshaping for vectorization
    weights_rep = log(weights_rep);
    
    % Within-trial FR weights
    weights2 = mean(spikes,1);
    weights2_rep = repmat(weights2,[size(FR,1) 1 size(spikes,1)]);
    weights2_rep = permute(weights2_rep,[1 3 2]); % Reshaping for vectorization
    weights2_rep = log(weights2_rep + eps);
    
    mm = cost_matrix + alpha*( ...
        beta*weights_rep + (1-beta)*weights2_rep);
    
    acc_cost_matrix = dynamic_function2(sum(mm,3), transition_priors);
    
    cost_matrix(:,:,end+1) = sum(mm,3);
end 

%% Get path

[path, path_matrix] = get_path(pointer_matrix, acc_cost_matrix);

%% Collapse cost matrix

cost_matrix = sum(cost_matrix,3);

%% Score function vectorized - Poisson

    function [ score ] = score_function_vec (FR, spikes)
        
        nn = size(FR,2);
        
        FR_rep = repmat(FR,[1 1 size(spikes,1)]);
        FR_rep = permute(FR_rep,[1 3 2]);
        
        spikes_rep = repmat(spikes,[1 1 size(FR,1)]);
        spikes_rep = permute(spikes_rep,[3 1 2]);
        
        score = spikes_rep.*log(FR_rep) - FR_rep - log(factorial(spikes_rep));
        
    end

%% Dynamic function

    function [acc_cost_matrix, pointer_matrix] = dynamic_function2(cost_matrix, transition_priors)
        
        padding_needed = max(size(transition_priors));
        
        acc_cost_matrix = -inf(size(cost_matrix,1) + padding_needed,size(cost_matrix,2) + padding_needed); % make it bigger in each dimension for padding
        pointer_matrix = cell(size(cost_matrix,1) + padding_needed,size(cost_matrix,2) + padding_needed);
        
        num_rows = size(cost_matrix,1);
        num_cols = size(cost_matrix,2);
        
        % Pad cost matrix
        cost_matrix_pad = acc_cost_matrix;
        cost_matrix_pad(1:num_rows,(padding_needed+1):end) = cost_matrix;
        
        acc_cost_matrix((num_rows+1),padding_needed) = 0; % place to start that's not -inf
         
        for col_num = 1:num_cols
            pad_col_num = col_num + padding_needed;
            
            for row_num = num_rows:-1:1
                pad_row_num = row_num;
            
                % Reset off-corner entry; was only to initiate corner
                if ~(col_num==1 && row_num==num_rows)
                    acc_cost_matrix((num_rows+1),padding_needed) = -inf;
                end
                
                % Function that calculates accumulated cost
                [acc_cost_matrix(pad_row_num,pad_col_num), pointer_matrix{pad_row_num,pad_col_num}] = calc_acc_entry2( ...
                            acc_cost_matrix, cost_matrix_pad, transition_priors, [pad_row_num, pad_col_num]);
                
            end
        end
    
        % Remove padding
        acc_cost_matrix(:,1:padding_needed) = []; 
        acc_cost_matrix((num_rows+1):(num_rows+padding_needed),:) = []; 
        pointer_matrix(:,1:padding_needed) = [];
        pointer_matrix((num_rows+1):(num_rows+padding_needed),:) =  []; 
        
        % Remove padding from individual entries in pointer matrix
        pointer_matrix = cellfun(@(x)x-[0,padding_needed],pointer_matrix,'UniformOutput',false);
        
    end

%% Calculate entry

    function [ entry, pointer ] = calc_acc_entry2(acc_cost_matrix, cost_matrix, transition_priors, coord)
        
        % Find connected nodes
        [connected_nodes, transition_scores] = which_are_connected(coord, transition_priors, size(cost_matrix));
        num_connected_nodes = size(connected_nodes,1);
        temp_scores = nan(num_connected_nodes,1); % initialize
        
        % Calculate accumulated scores for each connected node
        for node = 1:num_connected_nodes
            conn_node_coord = connected_nodes(node,:);
            transition_score = transition_scores(node);
            
            temp_scores(node) = acc_cost_matrix(conn_node_coord(1),conn_node_coord(2)) ...
                + transition_score ... 
                + cost_matrix(coord(1),coord(2)); % total score = accumulated cost of connected node + transition + current node llhd
        end
        
        entry = max(temp_scores);
        [~,idx_max] = max(temp_scores);
        pointer = connected_nodes(idx_max,:); % save pointer to best preceding node
        
    end

%% Function for: are two nodes connected?

    function [ conn_node_coords, transition_scores ] = which_are_connected( node_coord, transition_priors, matrix_size )
        
        % Go through all possible transitions
        conn_node_coords = cell(size(transition_priors));
        for i = size(transition_priors,1):-1:1
            for j = 1:size(transition_priors,2)
                new_coord(1) = node_coord(1) + (size(transition_priors,1)-i);
                new_coord(2) = node_coord(2) - (j-1);
                
                conn_node_coords{i,j} = new_coord;
            end
        end
        
        conn_node_coords = cell2mat(reshape(conn_node_coords,[],1));
        transition_scores = reshape(transition_priors,[],1);
        
        % Get rid of self transition
        conn_node_coords(size(transition_priors,1),:) = [];
        transition_scores(size(transition_priors,1)) = [];
        
    end
   
%% Function to get path
    function [path, path_matrix] = get_path(pointer_matrix, acc_cost_matrix)
        
        path_matrix = zeros(size(acc_cost_matrix));
        
        num_rows = size(acc_cost_matrix,1);
        num_cols = size(acc_cost_matrix,2);
        
        % initialize coordinates
        curr_coord = [1, num_cols];
        end_coord = [num_rows, 1];
        counter = 2; % first is already written
        
        path = [1, num_cols];
        path_matrix(num_rows,1) = 1;
        path_matrix(1,num_cols) = 1;
        
        while ~all(curr_coord == end_coord)
            next_step = pointer_matrix{curr_coord(1), curr_coord(2)};
            path(counter,1) = next_step(1);
            path(counter,2) = next_step(2);
            path_matrix(curr_coord(1), curr_coord(2)) = 1;
            curr_coord = next_step;
            counter = counter + 1;
        end
        
    end        

end

