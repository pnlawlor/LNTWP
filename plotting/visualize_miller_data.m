%% Load scripts

addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Scripts\Pat'))
addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Scripts_TW'))

%% Choose data
Brian_data = 0;

%% Load data

if Brian_data
% fpath = 'C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Source_data\Miller_lab\Processed\Brian_uncertainty\';
% load([fpath 'Brian_uncertainty_correct_only_v3.mat'])
% load([fpath 'PD_fit_v2.mat'])
else
% New path    
fpath = 'C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Source_data\Miller_lab\Processed\Matt\UnWindowed\bin_sizes\'; % Monkey 1

% Old path
% fpath = 'C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Source_data\Miller_lab\Processed\Matt\'; % Monkey 1

% By trial
% load([fpath 'MRT_CO_v3_by_trial.mat']) % Individual trials (<= 4 reaches)

% Individual reaches
% New files
%  Monkey 1
% load([fpath 'M1_RT_dt5ms.mat']) % Monkey 1; dt 5ms
% load([fpath 'M1_RT_dt10ms.mat']) % Monkey 1; dt 10ms
% load([fpath 'M1_RT_dt20ms.mat']) % Monkey 1; dt 20ms
% load([fpath 'M1_RT_dt40ms.mat']) % Monkey 1; dt 40ms

% Monkey 2
load([fpath 'M2_RT_dt10ms.mat']) % Monkey 1; dt 10ms

% Old files
% load([fpath 'MRT_CO_v2.mat']) % Monkey 1
% load([fpath 'RT_MrT_f1_v1.mat']); Data = Data2; % Monkey 2 file 1
% load([fpath 'RT_MrT_f2_v1.mat']); Data = Data2; % Monkey 2 file 2
% load([fpath 'RT_MrT_f3_v1.mat']); Data = Data2; % Monkey 2 file 3

end

% num_units_M1 = size(Data(1).neural_data_M1{1},1);
num_units_PMd = size(Data(1).neural_data_PMd{1},1);
num_trials_total = length(Data(1).neural_data_PMd);

%% User inputs (Brian) - visualize by trial
if Brian_data
% block_num = 2;
% dt = .01;
% 
% target_loc = [1:8]*2*pi/8;
% 
% trials_desired = 1:num_trials_total;
% trials_bad = [];
% % trials_bad = [1 19 28 32 37 43 45 59 75 125 145 160 202 206 213 217 233 243];
% trials = setdiff(trials_desired,trials_bad);
% 
% neurons_M1 = 1:num_units_M1;    
% neurons_PMd = 1:100;
% sort_nrn_by_PD = 0;
% 
% thresh_speed = 8;
% min_speed = 12;
% 
% time_window = [-10:20]; % With respect to speed threshold crossing
% % time_window = nan;
% % time_window = [0:150];
% 
% from_cue_on = 0;
% from_move_start = 1;
% 
% save_data = 0;
% plot_figs = 1;
end

%% User inputs (Matt) - visualize by trial
if ~Brian_data

block_num = 1;
dt = .01;

target_loc = [-7:2:7]*pi/8; % -pi:pi for Matt's data
% target_loc = [-4]*pi/8; % -pi:pi for Matt's data

trials_desired = 1:num_trials_total;
% trials_bad = [48]; % Bad trials: monkey 1
trials_bad = [];

trials = setdiff(trials_desired,trials_bad);

% neurons_M1 = 1:num_units_M1;    
neurons_PMd = 1:num_units_PMd;
sort_nrn_by_PD = 0;

thresh_speed = 12;
min_speed = 12; % old value was 12

time_window_sec = [-0.6,0.4]; % With respect to speed threshold crossing
time_window_bin = round(time_window_sec/dt); % Convert to bins
time_window = [time_window_bin(1):time_window_bin(2)]; % Make array

enforce_window_len = 0; % whether window must be same size as time_window

from_cue_on = 0;
from_move_start = 1;

save_data = 1;
plot_figs = 0;

min_rxn_time = 0.1; % old value was .1
wind_acc = 0.1; % amount of time in sec after threshold crossing that speed must not decrease (filters bad reaches)
wind_acc_tol = 0.5; % percentage of bind within wind_acc that must be accelerating (filters bad reaches)

end

%% Visualize by trial

% Choose trials with a given target location
flag_speed = 0;
trials_temp = zeros(num_trials_total,1);
trials_temp(trials) = 1;

if Brian_data
% %     target_locs_rep = repmat(Data(block_num).block_info(:,2),[1 length(target_loc)]);
%     target_locs_rep = repmat(Data(block_num).block_info(trials_desired,2),[1 length(target_loc)]);
% %     cond = abs(target_locs_rep - repmat(target_loc,[num_trials_total 1]));
%     cond = abs(target_locs_rep - repmat(target_loc,[length(trials_desired) 1]));
%     trials_temp_target_loc = any(cond<.05,2);
% %     trials = find(trials_temp.*trials_temp_target_loc);
%     trials = find(trials_temp_target_loc);
else
    temp = cell2mat(Data.reach_dir);
    target_locs_rep = repmat(temp,[1 length(target_loc)]);
    cond = abs(target_locs_rep - repmat(target_loc,[num_trials_total 1]));
    trials_temp_target_loc = any((cond<pi/8),2); % .39 is 2pi/8/2. I.e., half of a 2pi/8 bin in either direction
    trials = find(trials_temp.*trials_temp_target_loc);    
end

% Initialize storage variables
% spikes_M1 = cell(length(trials),length(neurons_M1));
spikes_PMd = cell(length(trials),length(neurons_PMd));
movement_cov = cell(length(trials),1);
kinematics = cell(length(trials),1);
go_cue = cell(length(trials),1);
outer_target = cell(length(trials),1);
target_dir = cell(length(trials),1);
target_dir_num = cell(length(trials),1);

if sort_nrn_by_PD
%     [~,idx_sort_M1] = sort(PD_M1,'ascend');
    [~,idx_sort_PMd] = sort(PD_PMd,'ascend');
else
%     idx_sort_M1 = 1:num_units_M1;
    idx_sort_PMd = 1:num_units_PMd;
end

for trial_idx = 1:length(trials)
    trial_num = trials(trial_idx);
    num_bins = size(Data(block_num).neural_data_PMd{trial_num}(neurons_PMd,:),2);
    
    x_pos = Data(block_num).kinematics{trial_num}(:,1);
    y_pos = Data(block_num).kinematics{trial_num}(:,2);
    x_vel = Data(block_num).kinematics{trial_num}(:,3);
    y_vel = Data(block_num).kinematics{trial_num}(:,4);
    x_acc = Data(block_num).kinematics{trial_num}(:,5);
    y_acc = Data(block_num).kinematics{trial_num}(:,6);
    speed = sqrt(x_vel.^2 + y_vel.^2);
    
    % Get movement onset
    if from_move_start
        if any(speed > thresh_speed)
            if ~Brian_data
                idx_target_on = find(Data.target_on{trial_num}); % When the target came on (in bins)
                bins_after_tgt = ceil(min_rxn_time/dt); % Min rxn time (in bins)
                idx_possible = idx_target_on + bins_after_tgt; % define earliest bin meeting rxn time criteria
                bins_possible = speed*0; % initialize variable
                bins_possible(idx_possible:end) = 1; % set which pins are candidates for inclusion
                bins_thresh_cross = (speed > thresh_speed).*(bins_possible); % bins that both exceed spd threshold and are after target onset + rxn time
                
                if sum(bins_thresh_cross) < 1 % if nothing meets critera, skip trial
                    continue
                end               
                
                bin_align = find(bins_thresh_cross,1); % Set bin to align data to
%                 disp(['Trial: ' num2str(trial_num) ', bin align: ' num2str(bin_align)])
            else
%                 bins_thresh_cross = (speed > thresh_speed);
%                 bin_align = find(speed > thresh_speed,1); % Set bin to align data to
            end
        else
            continue
            flag_speed = 1;
        end
    end
    
    % If you want from cue on
    if from_cue_on
        bin_align = 31;
    end
        
    
    % Change plotting time window
    if ~any(isnan(time_window))
        time_window_trial = time_window + bin_align; % For one reach, define data window relative to alignment bin
%         time_window_trial = 1:num_bins; % For whole trial
        
        % Do some error checking
        time_window_trial(time_window_trial<1) = []; % set any bins earlier than time window start to empty
        time_window_trial(time_window_trial>num_bins) = []; % set any bins later than time window end to empty
    else
        time_window_trial = 1:num_bins;
    end
    
    % Exclusion condition: speed must (mostly) monotonically increase for a fixed
    % amount of time (wind_acc)
    speed_wind_acc = speed(bin_align:(bin_align + round(wind_acc/dt)));
    perc_dec = sum((diff(speed_wind_acc)<0))/length(speed_wind_acc);
%     disp(['Trial: ' num2str(trial_num) ' . Percentage deceleration: ' num2str(perc_dec)])
    if  perc_dec > wind_acc_tol
        continue
    end
    
    % Exclusion condition: reaches that don't meet minimum speed
%     if max(speed(time_window_trial)) < min_speed
    if max(speed_wind_acc) < min_speed
        continue
    end
    
    % Open plotting window
    if plot_figs
    figure('units','normalized','outerposition',[0 0 .4 1])
%     figure
    end
    
    % Spatial kinematics
    if plot_figs
        subplot(6,4,[1 5 9 13 17 21])
        plot(x_pos(time_window_trial),y_pos(time_window_trial))
        axis([-12 12 -12 12])
        axis square
    end
    
    % X and Y kinematic
    % 1 = xpos, 3 = xvel, 5 = xacc
    if plot_figs
        subplot(6,4,[2 3 4])
        plot(x_vel(time_window_trial),'b'); hold on
        plot(y_vel(time_window_trial),'r'); hold off
        ylim([-30 30])
        xlim([0 length(time_window_trial)])
        title(['Trial: ' num2str(trial_num)])
        if flag_speed
            ylabel(['Speed threshold not reached!'])
        end
    end
    
    % Speed
%     if Brian_data
%     time_temp = Data(block_num).block_info(trial_num,6) - Data(block_num).block_info(trial_num,5);
%     bin_go_cue = round(time_temp/dt);
%     end
    
    if plot_figs
        subplot(6,4,[6 7 8])
        plot(speed(time_window_trial),'LineWidth',2)

        % Threshold line
        line([0 length(time_window_trial)],[thresh_speed thresh_speed],'Color','g')
        % Threshold crossing line
        thresh_cross_to_plot = find(bins_thresh_cross(time_window_trial),1);
        line([thresh_cross_to_plot thresh_cross_to_plot],[-30 30],'Color','g')


        % Not sure what ths is for - when aligning to cue?
        if isnan(time_window)
%             line([bin_align bin_align],[-30 30],'Color','g')
        else
            % Where to put the crossing line depends on the alignmnet
            if time_window(1) < -100
%                 line([bin_align bin_align],[-30 30],'Color','g') % For when time_window(1) is big and negative
            else
    %             line([(-time_window(1)+1) (-time_window(1)+1)],[-30 30],'Color','g') % For when time_window(1) is within the trial
            end
        end

        % Go cue line
        if Brian_data
    %     line([(bin_go_cue - time_window_trial(1) - 1) (bin_go_cue - time_window_trial(1) - 1)],[-30 30],'Color','m')
        end    

        ylim([-5 30])
        xlim([0 length(time_window_trial)])
        if flag_speed
            ylabel(['Speed threshold not reached!'])
        end
    end
    
    % M1 neural data
    if plot_figs
    subplot(6,4,[10 11 12 14 15 16])
%     imagesc(Data(block_num).neural_data_M1{trial_num}(idx_sort_M1(neurons_M1),time_window_trial),[0 3])
    xlabel('Time (10ms bins)'), ylabel('M1 N#')
    
    % PMd neural data
    subplot(6,4,[18 19 20 22 23 24])
    imagesc(Data(block_num).neural_data_PMd{trial_num}(idx_sort_PMd(neurons_PMd),time_window_trial),[0 3])
    xlabel('Time (10ms bins)'), ylabel('PMd N#')
    end
    
    if save_data
        % Save everything
%         spikes_M1(trial_idx,:) = mat2cell(Data(block_num).neural_data_M1{trial_num}(neurons_M1,time_window_trial)',[length(time_window_trial)], [ones(1,length(neurons_M1))]);
%         spikes_M1(trial_idx,:) = mat2cell(Data(block_num).neural_data_M1{trial_num}(idx_sort_M1(neurons_M1),time_window_trial)',[length(time_window_trial)], [ones(1,length(neurons_M1))]);
        spikes_PMd(trial_idx,:) = mat2cell(Data(block_num).neural_data_PMd{trial_num}(neurons_PMd,time_window_trial)',[length(time_window_trial)], [ones(1,length(neurons_PMd))]);
%         spikes_PMd(trial_idx,:) = mat2cell(Data(block_num).neural_data_PMd{trial_num}(idx_sort_PMd(neurons_PMd),time_window_trial)',[length(time_window_trial)], [ones(1,length(neurons_PMd))]);
        
        if Brian_data
%         temp = zeros(size(Data(block_num).neural_data_M1{trial_num},2),2);
%         temp(bin_align,:) = [cos(Data(block_num).block_info(trial_num,2)) sin(Data(block_num).block_info(trial_num,2))];
        else
            temp = zeros(size(Data(block_num).neural_data_PMd{trial_num},2),3); % initialize variable
            temp(bin_align,:) = [1 cos(Data.reach_dir{trial_num}) sin(Data.reach_dir{trial_num})]; % set value of 'bin_align' to indicate when movement started and its direction
        end
        
        movement_cov{trial_idx} = temp(time_window_trial,:); % include movement data from time_window
        kinematics{trial_idx} = [x_pos(time_window_trial) y_pos(time_window_trial), ... % include kinematics from time_window
                                    x_vel(time_window_trial) y_vel(time_window_trial), ...
                                    x_acc(time_window_trial) y_acc(time_window_trial)];
        if Brian_data
%         temp2 = temp(:,1)*0;
%         temp2(bin_go_cue) = 1;
%         go_cue{trial_idx} = temp2(time_window_trial);
        end

        % Write reach direction to file
        temp3 = repmat(temp(:,1)*0,[1 3]);
        temp3(1,1) = 1;
        if Brian_data
%         temp3(1,2:3) = [cos(Data(block_num).block_info(trial_num,2)) sin(Data(block_num).block_info(trial_num,2))];
%         target_dir{trial_idx} = Data(block_num).block_info(trial_num,2);
%         target_dir_num{trial_idx} = round(target_dir{trial_idx}/(pi/4));
        else
            temp3(1,2:3) = [cos(Data.reach_dir{trial_num}) sin(Data.reach_dir{trial_num})];
            target_dir{trial_idx} = Data.reach_dir{trial_num};
            target_dir_num{trial_idx} = round(target_dir{trial_idx}/(pi/8));
        end
        
        outer_target{trial_idx} = temp3(time_window_trial,:);
        
    end
end

if save_data
    
    if enforce_window_len
        full_trial_len = abs(time_window(1)) + abs(abs(time_window(end))) + 1;
        good_trials = cellfun( ...
            @(x)(size(x,1)==full_trial_len), ...
            spikes_PMd(:,1), ...
            'UniformOutput',false);
        good_trials = cell2mat(good_trials);
    else    
        good_trials = cellfun( ...
            @isempty, ...
            spikes_PMd(:,1), ...
            'UniformOutput',false);
        good_trials = ~cell2mat(good_trials);
    end
    
    % Write info to file
%     Data_win(block_num).spikes_M1 = spikes_M1(good_trials,:);
    Data_win(block_num).spikes_PMd = spikes_PMd(good_trials,:);
    Data_win(block_num).movement_cov = movement_cov(good_trials);
    Data_win(block_num).kinematics = kinematics(good_trials);
    Data_win(block_num).go_cue = go_cue(good_trials);
    Data_win(block_num).outer_target = outer_target(good_trials);
    Data_win(block_num).target_dir_num = target_dir_num(good_trials);
    Data_win(block_num).target_dir = target_dir(good_trials);
end

%% User inputs - visualize by neuron

block_num = 1;

target_loc = [-7:2:7]*pi/8; % -pi:pi for Matt's data
% target_loc = [-7:2:-3]*pi/8;
% target_loc = [1:8]*2*pi/8; % 0:2pi for Brian's data

min_reach_len = 0;

trials_desired = 1:num_trials_total;
% trials_bad = [1 19 28 32 37 43 45 59 75 125 145 160 202 206 213 217 233 243]; % Brian
% trials_bad = [23 32 34 69 139 122 110 195 196 194 186 180 164 213 202 284 306 382 368 437 423 466 48 36 12 5 150 146 132 126 107 182 179 232 224 203 291 350 331 400 355 353 477 469]; % Matt
trials_bad = [];
trials_to_use = setdiff(trials_desired,trials_bad);

neurons_M1 = [];    
neurons_PMd = [41:60];
thresh_speed = 8;

from_cue_on = 0;
from_move_start = 1;

% time_window = [1:30]; % With respect to speed threshold crossing
time_window = [-40:20];
% time_window = nan;

box_reach_st = [-100 100 -100 100]; % To define where reaches must start (e.g., center). [x1 x2 y1 y2], x2>x1,y2>y1; For Matt, x [-11,7], y [-7, 11]
% box_reach_st = [-4.5 .5 -.5 4.5]; % This is the "center box"
% box_reach_st = [-6 2 -2 6]; % This is the "center box"

min_rxn_time = .1;
dt = .01;

plot_figs = 1;
sort_trials = 1;

%% Visualize spatial tuning curves

% Choose trials with a given target location
trials_temp = zeros(num_trials_total,1);
trials_temp(trials_to_use) = 1;

% Choose reaches into the right directions
if Brian_data
    target_locs_rep = repmat(Data(block_num).block_info(:,2),[1 length(target_loc)]);
    cond = abs(target_locs_rep - repmat(target_loc,[num_trials_total 1]));
    trials_temp_target_loc = any(cond<.05,2);
else
    target_locs_temp = cell2mat(Data.reach_dir);
    target_locs_rep = repmat(target_locs_temp,[1 length(target_loc)]);    
    cond = abs(target_locs_rep - repmat(target_loc,[num_trials_total 1]));
    trials_temp_target_loc = any(cond<=pi/8,2);
end

% Choose reaches starting within a certain box
if Brian_data
    trials = find(trials_temp.*trials_temp_target_loc);
else
    % If you want first reach only
%     reach_nums = cell2mat(Data.reach_num);
%     trials_reach1 = reach_nums == 1;
    
    % If you want reaches of a certiain length
    reach_lens = cell2mat(Data.reach_len);
    trials_reach_len = reach_lens > min_reach_len;
    
    trials_box = reach_by_st_loc(Data(block_num),box_reach_st);
%     trials = find(trials_temp.*trials_temp_target_loc.*trials_box.*trials_reach1);
    trials = find(trials_temp.*trials_temp_target_loc.*trials_box.*trials_reach_len);
%     trials = find(trials_temp.*trials_temp_target_loc);
    
end

neurons_all = [neurons_M1 neurons_PMd];

% Initialize
neural_data_cell = cell(length(neurons_all),1);
PD = nan(length(neurons_all),1);

for nrn_idx = 1:length(neurons_all)
    nrn_num = neurons_all(nrn_idx);
    
    if nrn_idx > length(neurons_M1)
        region = 'PMd';
    else
        region = 'M1';
    end
    
    switch region
        case 'PMd'
            for trial_idx = 1:length(trials)
                trial_num = trials(trial_idx);
                
                Data_struct(trial_idx).data = full(Data(block_num).neural_data_PMd{trial_num}(nrn_num,:));
                Data_struct(trial_idx).align_bin = zeros(1,length(Data_struct(trial_idx).data));
                
                if Brian_data
                    Data_struct(trial_idx).target_dir = Data(block_num).block_info(trial_num,2);
                else
                    Data_struct(trial_idx).target_dir = Data(block_num).reach_dir{trial_num}; 
                end
                
                vel_x = Data(block_num).kinematics{trial_num}(:,3);
                vel_y = Data(block_num).kinematics{trial_num}(:,4);
                speed = smooth(sqrt(vel_x.^2 + vel_y.^2));
                
                if from_move_start
                    if any(speed > thresh_speed)
                        if ~Brian_data
                            idx_target_on = find(Data.target_on{trial_num}); % When the target came on (in bins)
                            bins_after_tgt = ceil(min_rxn_time/dt); % Min rxn time (in bins)
                            idx_possible = idx_target_on + bins_after_tgt;
                            bins_possible = speed*0;
                            bins_possible(idx_possible:end) = 1;
                            bins_thresh_cross = (speed > thresh_speed).*(bins_possible);
                            if sum(bins_thresh_cross) < 1
                                continue
                            end
                            bin_align = find(bins_thresh_cross,1); % For when speed threshold is reached; Need to push forward to around rxn time so as not to get end of current movement
                        else
                            bin_align = find(speed > thresh_speed,1); % For when speed threshold is reached; Need to push forward to around rxn time so as not to get end of current movement
                        end
                    else
                        continue
                    end
                end
                
                if from_cue_on
%                     bin_align = 1;
                    if ~Brian_data
                        bin_align = find(Data.target_on{trial_num});
                    else
                        bin_align = 1;
                    end
                end
                
                Data_struct(trial_idx).align_bin = bin_align;
            end
        case 'M1'
            for trial_idx = 1:length(trials)
                trial_num = trials(trial_idx);
                
                Data_struct(trial_idx).data = full(Data(block_num).neural_data_M1{trial_num}(nrn_num,:));
                Data_struct(trial_idx).align_bin = zeros(1,length(Data_struct(trial_idx).data));
                
                if Brian_data
                    Data_struct(trial_idx).target_dir = Data(block_num).block_info(trial_num,2);
                else
                    Data_struct(trial_idx).target_dir = Data(block_num).reach_dir{trial_num};
                end
                
                vel_x = Data(block_num).kinematics{trial_num}(:,3);
                vel_y = Data(block_num).kinematics{trial_num}(:,4);
                speed = smooth(sqrt(vel_x.^2 + vel_y.^2));
                
                if from_move_start
                    if any(speed > thresh_speed)
                        if ~Brian_data
                            idx_target_on = find(Data.target_on{trial_num}); % When the target came on (in bins)
                            bins_after_tgt = ceil(min_rxn_time/dt); % Min rxn time (in bins)
                            idx_possible = idx_target_on + bins_after_tgt;
                            bins_possible = speed*0;
                            bins_possible(idx_possible:end) = 1;
                            bins_thresh_cross = (speed > thresh_speed).*(bins_possible);
                            if sum(bins_thresh_cross) < 1
                                continue
                            end
                            bin_align = find(bins_thresh_cross,1); % For when speed threshold is reached; Need to push forward to around rxn time so as not to get end of current movement
                        else
                            bin_align = find(speed > thresh_speed,1); % For when speed threshold is reached; Need to push forward to around rxn time so as not to get end of current movement
                        end
                    else
                        continue
                    end
                end
                
                if from_cue_on
                    if ~Brian_data
                        bin_align = find(Data.target_on{trial_num});
                    else
                        bin_align = 1;
                    end
                end
                
                Data_struct(trial_idx).align_bin = bin_align;
            end
    end
    
    % Align data on whichever signal (target on, movement start, etc)
    [neural_data_aligned, aligned_bin] = align_data(Data_struct);
    
    % Group neural data by target direction
    neural_data_grouped = cell(length(target_loc),1);
    neural_data_mean = neural_data_grouped;
    x_fit = cell(length(target_loc),1);
    
    if plot_figs
    figure
%     figure('units','normalized','outerposition',[0 0 .95 1])
    subplot(1,2,2)
    end
    
    for tgt_idx = 1:length(target_loc)
        if Brian_data
            trials_correct_target = find(abs(cond(trials,tgt_idx)) < .05); % Note that trials indexes are a little weird here. This is a subset of the already selected trials.
        else
            trials_correct_target = find(abs(cond(trials,tgt_idx)) <= pi/8); % Note that trials indexes are a little weird here. This is a subset of the already selected trials.
        end
        
        if isnan(time_window)
            neural_data_grouped{tgt_idx} = neural_data_aligned(trials_correct_target,:);
        else
            time_window_trial = time_window + aligned_bin;
            time_window_trial(time_window_trial<1) = [];
            time_window_trial(time_window_trial>size(neural_data_aligned,2)) = [];
%             neural_data_grouped{tgt_idx} = neural_data_aligned(trials_correct_target,time_window + aligned_bin);
            neural_data_grouped{tgt_idx} = neural_data_aligned(trials_correct_target,time_window_trial);
        end
        
        neural_data_mean{tgt_idx} = 100*nanmean(neural_data_grouped{tgt_idx},2);
        
        % Add some jitter
        x_plot = repmat(target_loc(tgt_idx),length(trials_correct_target),1) + .1*(rand(length(trials_correct_target),1)-.5);
        x_fit{tgt_idx} = repmat(target_loc(tgt_idx),length(trials_correct_target),1);
        
        if plot_figs
        scatter(x_plot,neural_data_mean{tgt_idx},'r')
        hold on
        end
    end
    
    FR_mean = nanmean(nanmean(neural_data_aligned))*100;
    
    % Arrange data for fit
    x_fit_all = cell2mat(x_fit);
    neural_data_all = cell2mat(neural_data_mean);
    
    % Fit
    X = [cos(x_fit_all) sin(x_fit_all)];
    [fit_param, dev, stats] = glmfit(X,neural_data_all,'poisson');
    PD(nrn_idx) = atan2(fit_param(3),fit_param(2));
    
    % Prepare plot for fit function
    if Brian_data
        x_plot_temp = [0:.01:2*pi]';
    else
        x_plot_temp = [-pi:.01:pi]';
    end
    x_plot_cos = cos(x_plot_temp);
    x_plot_sin = sin(x_plot_temp);
    y = glmval(fit_param, ...
        [x_plot_cos, x_plot_sin], ...
        'log');
    
    % Plot tuning data
    if plot_figs
    plot(x_plot_temp,y,'b')
    
    if Brian_data
        axis([0 7 0 60])
    else
        axis([-3.5 3.5 0 60])
    end
    
    xlabel('Target direction (rad)')
    ylabel('Mean FR (Hz)')
    title([region ' neuron: ' num2str(nrn_num) ', center: ' num2str(PD(nrn_idx)*180/pi) ])
    end
    
    % Plot rasters
    % Sort neural data by target direction
    if Brian_data
        [dir_reach,idx_sorted_trials] = sort(Data(block_num).block_info(trials,2));
    else
        [dir_reach,idx_sorted_trials] = sort(cell2mat(Data.reach_dir(trials)));
        [len_reach,idx_sorted_trials_len] = sort(cell2mat(Data.reach_len(trials)));
        dir_reach_unsorted = cell2mat(Data.reach_dir(trials));
        len_reach_unsorted = cell2mat(Data.reach_len(trials));
    end
    
    % Apply sort to trials if necessary
    if isnan(time_window)
        if sort_trials
            neural_data_aligned = neural_data_aligned(idx_sorted_trials,:);
        end
        neural_data_cell{nrn_idx} = neural_data_aligned;
    else
        if sort_trials
            neural_data_aligned = neural_data_aligned(idx_sorted_trials,time_window_trial); % Sort trials by reach direction
%           neural_data_aligned = neural_data_aligned(idx_sorted_trials_len,time_window_trial); % Sort trials by reach length
        end
        neural_data_cell{nrn_idx} = neural_data_aligned;
    end
    
    if plot_figs
    subplot(1,2,1)
    imagesc(neural_data_aligned,[0 3])
%     imagesc(circshift(neural_data_aligned,20),[0 3])
    
    if isnan(time_window)
        title([region ' neuron: ' num2str(nrn_num) '. Mean FR: ' num2str(FR_mean) ', aligned to bin: ' num2str(aligned_bin)])
    else
        title([region ' neuron: ' num2str(nrn_num) '. Mean FR: ' num2str(FR_mean) ', aligned to bin: ' num2str(-time_window(1))])
    end
    
    hold off
    end
    
    disp(['Completed: ' region ' neuron ' num2str(nrn_num)])

end
