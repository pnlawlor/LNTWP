% Simulate single neuron data for time warping project

%% User inputs

num_trials = 20;
bin_size = 10; % is ms
dt = bin_size/1000;
dir_pref = pi/4; % in radians
dir_pref_width = pi/8; % In radians, the allowable range of movement angles around the preferred direction
bins_per_trial = 50;
width_range = 1:3;
baseline_term = 3; % in Hz
tuning_strength = 4; 

measured_onset = 20; % When the movement started (bin#)
onset_range = 10:30; % When the neural response to movement started (bin#)

visualize_data = 1;

%% Initialize

firing_rates = repmat({zeros(bins_per_trial,1)},num_trials,1);
spikes = firing_rates;
movement_data = repmat({zeros(bins_per_trial,2)},num_trials,1);
movement_data_for_sim = repmat({zeros(bins_per_trial,2)},num_trials,1);
X_cell = cell(num_trials,1);
movement_angle = cell(num_trials,1);

% Narrower basis set
% basis_temp = getBasis('rcos',10,20)';
% basis_temp = basis_temp(:,(end-2):end);

% Wider basis set
basis_temp = getBasis('rcos',10,40)';
basis_temp = basis_temp(:,(end-4):2:end);

basis = cell(3,1);
for i = 1:size(basis_temp,2)
    basis{i} = basis_temp(basis_temp(:,i)>0,i);
end

% Determine weights for preferred direction
w_x = tuning_strength*cos(dir_pref);
w_y = tuning_strength*sin(dir_pref);

% Form weight vector
beta = [baseline_term, w_x, w_y]';

%% Simulate data

% for each trial
for trial_num = 1:num_trials
    % Choose direction of reach
    movement_angle{trial_num} = dir_pref_width*2*(rand-.5) + dir_pref;
    cov_movement = [cos(movement_angle{trial_num}) sin(movement_angle{trial_num})];
    
    % Choose onset time of neural activity
    onset = onset_range(randi(length(onset_range))); % Choose random onset of neural activity
%     onset = 10; % Fixed onset
    
    % Choose width of temporal basis function
    width = width_range(randi(length(width_range)));
%     width = 3; % Fixed width

    % Filter movement data with temporal basis function
    movement_data{trial_num}(measured_onset,:) = cov_movement; % What is actually measured in the experiment
    movement_data_for_sim{trial_num}(onset,:) = cov_movement; % So that simulated neural data is wonky
    x_temp = filter(basis{width},1,movement_data_for_sim{trial_num}(:,1));
    y_temp = filter(basis{width},1,movement_data_for_sim{trial_num}(:,2));
    
    % Simulate data for this trial
    X_cell{trial_num} = [(1+x_temp*0),x_temp,y_temp];
    firing_rates{trial_num} = exp(X_cell{trial_num}*beta);
    spikes{trial_num} = poissrnd(firing_rates{trial_num}*dt);    
    
end

%% Visualize data

if visualize_data
    
    for tr = 1:num_trials
        figure
        plot(firing_rates{tr}*dt), hold on
        plot(spikes{tr}), hold off
        axis([0 length(firing_rates{tr}) 0 10])
        title(['Trial: ' num2str(tr) ', angle: ' num2str(movement_angle{tr}*180/pi)]);
    end
    
end