function [ acc_cost_matrix, cost_matrix, path, path_matrix ] = align_signals( FR, spikes, transition_priors )

%% Initialization

if nargin < 3
    transition_priors = [1, 1; 0, 1];
end

% transition_priors = transition_priors./sum(sum(transition_priors));

cost_matrix = nan(length(FR),length(spikes));
% global_path = nan();

%% Fill out cost matrix

FR = flipud(FR);

for iFR = 1:length(FR)
    for iSp = 1:length(spikes)
        cost_matrix(iFR,iSp) = score_function(FR(iFR),spikes(iSp));
    end
end

%% Fill out accumulated cost matrix

acc_cost_matrix = dynamic_function( cost_matrix, transition_priors);

%% Get path

[path, path_matrix] = get_path(acc_cost_matrix, transition_priors);

%% Score function - Poisson
    % This might only be for the likelihood
    
    function [ score ] = score_function( FR , spikes )
        
        score = sum(            ...
            spikes.*log(FR) - FR - log(factorial(spikes))       ... % From the Poisson model of spiking $\frac{\lambda^n}{n!}exp(-\lambda)$
            );
        
    end

%% Score function - Linear

    function [ score ] = score_function2( FR , spikes )
        
        score = sum(            ...
            -20*(FR - spikes).^2        ... % From the Poisson model of spiking $\frac{\lambda^n}{n!}exp(-\lambda)$
            );
        
    end

%% Dynamic programming function

    function [ acc_cost_matrix ] = dynamic_function( cost_matrix, transition_priors)
        % Current strategy is to fill out the "outer-most" layers first
        
        acc_cost_matrix = -inf(size(cost_matrix));
        
        num_layers = max(size(cost_matrix));
        
        for layer_num = 1:num_layers
            
            i_start = max(size(cost_matrix,2) - (layer_num - 1),1); % To make sure we stop at edge of matrix
            j_start = min(layer_num,size(cost_matrix,1));
            j_end = size(cost_matrix,1); 
            
            % For "spikes" layer
            if layer_num < size(cost_matrix,2)
                for i = i_start:-1:1
                    acc_cost_matrix(j_start,i) = calc_acc_entry(acc_cost_matrix, cost_matrix, transition_priors, [j_start i]);
                end
            end
            
            % For "FR" layer
            if layer_num < size(cost_matrix,1)
                for j = j_start:j_end
                    acc_cost_matrix(j,i_start) = calc_acc_entry(acc_cost_matrix, cost_matrix, transition_priors, [j i_start]);
                end
            end
            
        end
        
    end

%% Function to calculate entry in accumulated cost matrix

    function [ entry ] = calc_acc_entry(acc_cost_matrix, cost_matrix, transition_priors, coord)
        if (coord(1)==1) && (coord(2)==size(cost_matrix,2))
            % Assign top right corner of accumulated cost matrix
            entry = cost_matrix(coord(1),coord(2));
        else
            % Find connected nodes
            [connected_nodes, transition_scores] = which_are_connected(coord, transition_priors, size(cost_matrix));
            num_connected_nodes = size(connected_nodes,1);
            temp_scores = nan(num_connected_nodes,1);
            
            % Calculate accumulated scores for each connected node
            for node = 1:num_connected_nodes
                conn_node_coord = connected_nodes(node,:);
                transition_score = transition_scores(node);
                temp_scores(node) = acc_cost_matrix(conn_node_coord(1),conn_node_coord(2)) + log(transition_score) + cost_matrix(coord(1),coord(2)); % total score = accumulated cost of connected node + transition + current node llhd
            end
            
            % Compare best new scores for this node with anything existing
            if max(temp_scores) > acc_cost_matrix(coord(1),coord(2))
                entry = max(temp_scores);
            else
                entry = acc_cost_matrix(coord(1),coord(2));
            end
        end
    end

%% Function for: are two nodes connected?
% Build in ranking of nodes too 

    function [ conn_node_coords, transition_scores ] = which_are_connected( node_coord, transition_priors, matrix_size )
        % Initialize
        num_possible_trans = numel(transition_priors)-1;
        conn_node_coords = nan(num_possible_trans,2);
        transition_scores = nan(num_possible_trans,1);
        
        % Go through all possible transitions
        conn_node_coords = cell(size(transition_priors));
        for i = size(transition_priors,1):-1:1
            for j = 1:size(transition_priors,2)
                new_coord(1) = node_coord(1) - (size(transition_priors,1)-i);
                new_coord(2) = node_coord(2) + (j-1);
                
                conn_node_coords{i,j} = new_coord;
            end
        end
        
        conn_node_coords = cell2mat(reshape(conn_node_coords,[],1));
        transition_scores = reshape(transition_priors,[],1);
        
        % Get rid of self transition
        conn_node_coords(size(transition_priors,1),:) = [];
        transition_scores(size(transition_priors,1)) = [];
        
        % Clean up others
        good_entries = find((conn_node_coords(:,1) <= matrix_size(1)).*(conn_node_coords(:,1) > 0).*(conn_node_coords(:,2) <= matrix_size(2)).*(conn_node_coords(:,2) > 0));
        conn_node_coords = conn_node_coords(good_entries,:);
        transition_scores = transition_scores(good_entries);
        transition_scores(transition_scores==0) = transition_scores(transition_scores==0) - Inf;
    end

%% Function to get path

    function [ path, path_matrix ] = get_path( acc_matrix, transition_priors )
        % Initialize
        longest_possible_path = size(acc_matrix,1) + size(acc_matrix,2);
        path_matrix = zeros(size(acc_matrix));
        path = nan(longest_possible_path,2);
        allowed_trans = double(transition_priors > 0);
        allowed_trans(allowed_trans==0)= -Inf;
        
        % Pad matrix to account for edge cases
        pad_size = size(transition_priors);
        acc_matrix = padarray(acc_matrix,pad_size,-Inf);
        path(1,1) = size(acc_matrix,1) - pad_size(1);
        path(1,2) = 1 + pad_size(2);
        
        start_coord = [size(acc_matrix,1) 1] + [-pad_size(1) pad_size(2)]; % Start location 
        end_coord = [1 size(acc_matrix,2)] + [pad_size(1) -pad_size(2)]; % End location
        curr_coord = start_coord;
        counter = 2;
        
        while ~all(curr_coord == end_coord)
            % Get scores of allowed transitions
            allowed_cols = curr_coord(2) - 1 + [1:(size(transition_priors,2))]; allowed_cols = sort(allowed_cols,'ascend');
            allowed_rows = curr_coord(1) + 1 - [1:(size(transition_priors,1))]; allowed_rows = sort(allowed_rows,'ascend');
            scores = acc_matrix(allowed_rows,allowed_cols).*allowed_trans;
            scores(size(transition_priors,1),1) = -Inf;
            [~,temp] = max(scores(:));
            [best_transition(1), best_transition(2)] = ind2sub(size(scores),temp);
%             best_transition(2) = find(max(scores,[],1),1);
%             best_transition(1) = find(max(scores,[],2),1);
            
            % Change current location to best allowed transition
            curr_coord(2) = curr_coord(2) + best_transition(2) - 1;
            curr_coord(1) = curr_coord(1) - (size(transition_priors,1) - best_transition(1));
%             disp(curr_coord)
            
            if curr_coord(2) > size(acc_matrix,2) || curr_coord(1) > size(acc_matrix,1)
                disp('Theres an issue')
            end
            
            % Save in path
            path(counter,1) = curr_coord(1);
            path(counter,2) = curr_coord(2);
            
            % Update counter
            counter = counter + 1;
            
        end
        
        % Subtract pad coord
        path = path - repmat(pad_size,[size(path,1) 1]);
        
        % Fix any NaNs in path
        idx_last_value = find(~isnan(path(:,1)),1,'last');
        idx_nan = isnan(path(:,1));
        path(idx_nan,:) = repmat(path(idx_last_value,:),[sum(idx_nan) 1]);
        
        % Make path matrix
        for i = 1:size(path,1)
            path_matrix(path(i,1),path(i,2)) = 1;
        end
        
    end

end

