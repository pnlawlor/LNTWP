function [ warp_matrix ] = generate_warp( dimensions, transition_prob, run_len)

    %% Initialize
    warp_matrix = zeros(dimensions); % Initialize warp matrix
    transition_prob = transition_prob / sum(sum(transition_prob)); % Normalize probabilities to sum to 1
    transition_prob = flipud(transition_prob);
    prob_reshape = reshape(transition_prob,[],1);
    prob_reshape = prob_reshape(prob_reshape>0);
    
    % Note allowed steps
    [row_allowed, col_allowed] = find(transition_prob>0); % Gives subscripts of allowed probabilities
    row_allowed = -(row_allowed - 1); % Subtract 1 to get "difference". So that e.g., "go over 1" corresponds to second column.
    col_allowed = col_allowed - 1; % " ", but negative bc we need to 
    
    %% Define nether region
    nether_region = ones(dimensions);
    step_angles = atan2(-row_allowed,col_allowed);
    [~,steepest_step_idx] = max(step_angles);
    steepest_step = [row_allowed(steepest_step_idx), col_allowed(steepest_step_idx)];
    
    current_loc = [dimensions(1), 1];
    
    for col = 1:dimensions(2)
       
        % Set everything above it to 0
        nether_region(1:(current_loc(1)-1),col) = zeros((current_loc(1)-1),1);
        
        % Take steeptest step, update location
        current_loc = current_loc + steepest_step;
        
    end
    
    nether_region_temp = flipud(nether_region);
%     nether_region = tril(nether_region_temp) + tril(nether_region_temp)';
    nether_region = tril(nether_region_temp) + rot90(tril(nether_region_temp),2);
    nether_region = flipud(logical(nether_region));
    
    %% Simulate 
    
    % Start at lower left corner
    current_loc = [dimensions(1),1];
    warp_matrix(current_loc(1),current_loc(2)) = 1;
    step_num = 1;
    current_run = 0;
    
    % While current location is not at top right corner
    end_location = [1, dimensions(2)];
    
    while ~all(current_loc == end_location)
        
        % Check which steps are allowed given current position (NB edges)
        next_loc = repmat(current_loc,length(row_allowed),1) + [row_allowed,col_allowed]; % Add theoretical transitions to current location to see theoretical endpoints
        rows_good = (next_loc(:,1)>0).*(next_loc(:,1)<=dimensions(1)); % Allowed rows
        cols_good = (next_loc(:,2)>0).*(next_loc(:,2)<=dimensions(2)); % Allowed columns
        trans_good = logical(rows_good.*cols_good); % Both rows and columns that are allowed
        next_loc = next_loc(trans_good,:); % Limit to allowable transitions
        
        prob_reshape2 = prob_reshape(trans_good);

        nether_region_linear = nether_region(sub2ind(size(nether_region),next_loc(:,1),next_loc(:,2)));

        trans_good2 = logical(nether_region_linear); % Both rows and columns that are allowed % FIX THIS
        trans_good_idx = find(trans_good2);
        
        % Determine if run is allowed
        if (run_len > 0) && (step_num ~= 1)
            % Determine location of would-be run
            run_loc(1) = current_loc(1) + previous_step(1);
            run_loc(2) = current_loc(2) + previous_step(2);
            
            % If the step is within the matrix
            if (run_loc(1)<=0) || (run_loc(1)>dimensions(1)) || (run_loc(2)<=0) || (run_loc(2)>dimensions(2))
                run_allowed = 0;
            else
                % Check to see if the step is in the nether region
                run_allowed = logical(nether_region(run_loc(1),run_loc(2)));
            end
        else
            run_allowed = 0;
        end
        
        % Draw run if desired
        if (run_len > 0) && (current_run < run_len) && run_allowed
            
            % Use previous step
            row_step = previous_step(1);
            col_step = previous_step(2);
            
            current_run = current_run + 1; % Update current run length
            
        % Otherwise, draw random step given current position
        else
            trans_prob_good_reduced = prob_reshape2(trans_good2);
            trans_prob_good_reduced = trans_prob_good_reduced / sum(trans_prob_good_reduced); % Normalize
            draw_from_this = cumsum(trans_prob_good_reduced); % Generate "CDF" over allowed transitions for sampling
            draw_from_this = [0; draw_from_this]; % Tack on a 0 so that #intervals = #transitions
            draw = rand; % Random number between 0 and 1
            draw_result = find(histcounts(draw,draw_from_this)); % Outputs: of the allowable ones, which was chosen
            draw_idx = trans_good_idx(draw_result); % Finds the "draw_result"th allowable transition

            row_step = row_allowed(draw_idx);
            col_step = col_allowed(draw_idx);
            
            current_run = 0; % Update current run to 0 since this step was random.
        end
        
        % Update current location
        current_loc = current_loc + [row_step, col_step]; % Add transition to current location
        warp_matrix(current_loc(1),current_loc(2)) = 1;

        % Keep track of previous step
        previous_step = [row_step, col_step];
        step_num = step_num + 1;
        
    end

end

