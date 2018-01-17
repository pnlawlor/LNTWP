% Takes data files from Miller lab, and outputs cleaned files:
% 1) binned and aligned appropriately
% 2) desired neurons

%% Add paths

addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Scripts\Pat'))
addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Scripts_TW'))

%% Choose data
brian_data = 0;

%% Load data

if brian_data
%     fpath = 'C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Source_data\Miller_lab\Raw\Brian\';
%     
%     % Load neural data file
%     load([fpath 'alldays0806.mat'])
%     % load([fpath 'raw_spike_matrices.mat']) % To avoid doing the getSpkMat every time, which takes forever
%     % Load kinematics file
%     load([fpath '0806_kin.mat'])
else
%     fpath = 'C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Source_data\Miller_lab\Raw\Matt\data_for_pat_raw\'; % Matt RT; monkey 1
%     load([fpath 'RT_VR_BL_2014-01-16.mat']) % Matt RT; monkey 1
%     load([fpath 'CO_VR_BL_2014-03-03.mat']) % Matt CO
    
    fpath = 'C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Source_data\Miller_lab\Raw\Matt\MrT\'; % Matt RT; monkey 2
    load([fpath 'MrT_RT_VR_BL_2013-09-04.mat']) % Matt RT; monkey 2 file 1
%     load([fpath 'MrT_RT_VR_BL_2013-09-06.mat']) % Matt RT; monkey 2 file 2
%     load([fpath 'MrT_RT_VR_BL_2013-09-10.mat']) % Matt RT; monkey 2 file 3

    alldays.tt = trial_table;
end

%% User inputs - Brian data

if brian_data
    
% num_trials = 100;
% correct_trials_only = 1;
% dt = .01; % in seconds. .01 = 10ms
% FR_min = 5; % in Hz
% blocks = 1;%:3;
% 
% id_success = 32;
% id_failure = 34;
% 
% col_start_time = 5;
% col_end_time = 7;
% col_reach_dir = 2;
% 
% time_end_extend = .4;
% time_before_extend = -.3; % needs to be negative to go backwards

end

%% User inputs - Matt data

if ~brian_data

correct_trials_only = 0; % All included trials are correct
dt = .01; % in seconds. .01 = 10ms, .02 = 20ms, etc
blocks = 1;

col_start_time = 1; % RT task
col_end_time = 22; % RT task

% col_start_time = 8; % CO task
% col_end_time = 13; % CO task
% col_reach_dir = 3; % CO task

time_end_extend = .4;
time_before_extend = -.3; % needs to be negative to go backwards

end

%% Initialize

if brian_data

    % Get block start and end times
%     bl_st(1) = alldays(1).tt(1,col_start_time);
%     bl_end(1) = alldays(1).tt(end,col_end_time);
%     
%     bl_st(2) = alldays(2).tt(1,col_start_time);
%     bl_end(2) = alldays(2).tt(end,col_end_time);
%     
%     bl_st(3) = alldays(3).tt(1,col_start_time);
%     bl_end(3) = alldays(3).tt(end,col_end_time);
%     
%     % Get number of trials
%     if correct_trials_only
%         bl_trials{1} = (alldays(1).tt(:,8)==id_success);
%         bl_trials{2} = (alldays(2).tt(:,8)==id_success);
%         bl_trials{3} = (alldays(3).tt(:,8)==id_success);
%     else
%         bl_trials{1} = ones(size(alldays(1).tt,1),1);
%         bl_trials{2} = ones(size(alldays(2).tt,1),1);
%         bl_trials{3} = ones(size(alldays(3).tt,1),1);
%     end

else

    bl_trials{1} = ones(size(alldays.tt,1),1);
    
end

%% Extract data

if ~brian_data

    % For M1 units
%     for i = 1:length(M1.units)
%         M1_units{i} = M1.units(i).ts;
%     end
    
    % For PMd units
    for i = 1:length(PMd.units)
        PMd_units{i} = PMd.units(i).ts;
    end
    
    kin.pos(:,1) = cont.t;
    kin.pos(:,2:3) = cont.pos;
    kin.vel(:,1) = cont.t;
    kin.vel(:,2:3) = cont.vel;
    kin.acc(:,1) = cont.t;
    kin.acc(:,2:3) = cont.acc;
end


% Extract neural data
tic

% num_units_M1 = length(M1_units);
num_units_PMd = length(PMd_units);

% max_spk_time = max(cellfun(@max,[M1_units; PMd_units]));
% max_spk_time = max(kin.pos(:,1)) - 1;
max_spk_time = max(kin.pos(:,1));
    
edges = 0:dt:(max_spk_time + dt); % Defines edges for histogram

% New way - histcounts()
% neural_data_temp_M1 = nan(num_units_M1,length(edges)-1);
neural_data_temp_PMd = nan(num_units_PMd,length(edges)-1);

% for nrn_num = 1:num_units_M1
%     neural_data_temp_M1(nrn_num,:) = histcounts(M1_units{nrn_num},edges);
%     disp(['M1 unit: ' num2str(nrn_num) '/' num2str(num_units_M1)])
% end

for nrn_num = 1:num_units_PMd
    neural_data_temp_PMd(nrn_num,:) = histcounts(PMd_units{nrn_num},edges);
    disp(['PMd unit: ' num2str(nrn_num) '/' num2str(num_units_PMd)])
end

% IMPORTANT: 
% Spikes collected from 0 ms onward; Kinematics from 1000ms onward
% This gets rid of the spikes before kinematics  are collected
% I.e., with 10ms bins, 1000ms = 100 bins * 10ms/bin

num_bins_discard = round(1/dt); % number of bins to discard before kinematics are colletcted
bin_start = num_bins_discard + 1; % bin to start keeping track of data

% neural_data_temp_M1 = neural_data_temp_M1(:,bin_start:end);
neural_data_temp_PMd = neural_data_temp_PMd(:,bin_start:end);

disp(['Neural data extracted: ' num2str(toc) ' seconds'])

% Extract 
tic
kin_ts = downsample(kin.pos(:,1),round(1000*dt));  % This doesn't need to be filtered (it's just timestamps)
x_pos = decimate(kin.pos(:,2),round(1000*dt)); % Sampling rate of kinematics is 1 kHz. thus, decimate so bin size is same as dt
y_pos = decimate(kin.pos(:,3),round(1000*dt));
x_vel = decimate(kin.vel(:,2),round(1000*dt));
y_vel = decimate(kin.vel(:,3),round(1000*dt));
x_acc = decimate(kin.acc(:,2),round(1000*dt));
y_acc = decimate(kin.acc(:,3),round(1000*dt));
disp(['Kinematic data extracted: ' num2str(toc) ' seconds'])

% Make timestamps
ts = edges(1:(end-1));
ts = ts(bin_start:end);

%% Make cells

for block_idx = 1:length(blocks)
    block_num = blocks(block_idx);
    trials = find(bl_trials{block_num}==1);
    
    % Initialize
    num_trials = sum(bl_trials{block_num});
    Data(block_num).kinematics = cell(num_trials,1);
%     Data(block_num).neural_data_M1 = cell(num_trials,1);
    Data(block_num).neural_data_PMd = cell(num_trials,1);
    Data(block_num).block_info = alldays(block_num).tt(trials,:);
    Data(block_num).trials = bl_trials{block_num};
    
%     num_units_M1 = length(M1_units);
    num_units_PMd = length(PMd_units);
    
    % Arrange
    for trial_idx = 1:length(trials)
        trial_num = trials(trial_idx);
        
        disp(['Writing data for trial: ' num2str(trial_num)])
        
        tr_start = alldays(block_num).tt(trial_num,col_start_time) + time_before_extend;
        tr_end = alldays(block_num).tt(trial_num,col_end_time) + time_end_extend;
        tr_bins = logical((ts>=tr_start).*(ts<=tr_end));
        
        % Take neural data from between tr_start and tr_end
%         Data(block_num).neural_data_M1{trial_idx} = neural_data_temp_M1(:,tr_bins);
        Data(block_num).neural_data_PMd{trial_idx} = neural_data_temp_PMd(:,tr_bins);
        
        % Take kinematics from between tr_start and tr_end
        Data(block_num).kinematics{trial_idx} = [x_pos(tr_bins), y_pos(tr_bins), ...
                                    x_vel(tr_bins), y_vel(tr_bins), ...
                                    x_acc(tr_bins), y_acc(tr_bins), ...
                                    kin_ts(tr_bins)];
                                
        % For CO task
%         Data(block_num).reach_dir{trial_idx,1} = alldays(block_num).tt(trial_num,col_reach_dir);
%         Data(block_num).reach_len{trial_idx,1} = 10;
%         temp = zeros(sum(tr_bins),1); temp(1 + ceil(abs(time_before_extend)/dt)) = 1;
%         Data(block_num).target_on{trial_idx,1} = temp;
                                
        % Get timestamps (imposed by me) to make sure they match those of
        % the kinematics
        Data(block_num).timestamps{trial_idx,1} = ts(tr_bins)';
        
    end
end

%% Further arrange data - For data with multiple reaches per trial; break them up
% Note that the data extracted here will be further refined/windowed in
% visualize_miller_data.m, and during the analysis script. 
% I.e., this isn't the final data. This data window is larger/more
% permissive.

% More user inputs
min_reach_len = 2; % in cm
% hold_time = .4;
time_before_cue = -.3; % Amount of time before target comes on to keep track of data (in sec)
max_reach_time = round(1.4/dt) + ceil(abs(time_before_cue)/dt); % max time for reach
spd_thresh = 8; % in cm/sec
buff = .3; % Velocity often non-zero when cue comes on. Look forward at least this much to find end of reach (in sec)
end_buff = .3; % Allows reach end to be a little later than the official end of trial. Just to be a little more permissive (larger data window)
pd_lag = .096; % Photodetector wasn't used, so "cue on" is the command signal, not the detection signal. *Average* lag in Miller lab is 96 ms

% Initialization
Data2 = Data;
idx = 1;

for tr = 1:num_trials
        reaches = find(~isnan(Data.block_info(tr,[3 8 13 18])));
        num_reach = length(reaches);
%         num_bins = size(Data.neural_data_M1{tr},2);
%         tr_end = Data.block_info(tr,22);
        tr_end = Data.block_info(tr,col_end_time);

        ts = Data.kinematics{tr}(:,end); % Based upon kinematic time stamps

        x_vel = Data.kinematics{tr}(:,3);
        y_vel = Data.kinematics{tr}(:,4);
        x_pos = Data.kinematics{tr}(:,1);
        y_pos = Data.kinematics{tr}(:,2);
        
        spd = sqrt(x_vel.^2 + y_vel.^2);

        for reach_idx = 1:num_reach
            reach = reaches(reach_idx);

            idx_cue_on = 2 + 5*(reach - 1); % in trial table, when target comes on
            
            % Correct for command signal lag
            cue_on = Data.block_info(tr,idx_cue_on);
            cue_on = cue_on + pd_lag;
            
            % Correction: for reaches 2-4, target is displayed 100ms before
            % time in trial table
            if reach > 1
                cue_on = cue_on - .1; % get the actual time from trial table; subtract 100 ms for correction
            end

            
            wind_st = cue_on + time_before_cue; % when the data window starts

            [~,idx_cue_on2] = min(abs(ts-cue_on)); % find time bin when cue comes on

            % If reach wasn't completed or was invalid for some reason
            % (Columns referenced here will be NaN if reach is invalid)
            if isnan(Data.block_info(tr,idx_cue_on+1)) || isnan(Data.block_info(tr,idx_cue_on+2))
                continue
            end

            % Determine end time for each reach
            if reach < 4 % For reaches 1-3 on a given trial
                % Find time of: min velocity before next go cue
                idx_cue_on_next = 2 + 5*(reach);
                cue_on_next = Data.block_info(tr,idx_cue_on_next);
                
                % For reach to end: (this is permissive) 
                % 1) slow (falls below speed threshold), and
                % 2) a certain minimum time after cue onset must have elapsed (buff) 
                % 3) end has to be before the next reach starts plus buff
                cond_reach_end = logical((spd < spd_thresh).*(ts > (cue_on + buff)).*(ts < (cue_on_next + buff)));
                
                % If conditions not met, skip
                if sum(cond_reach_end) == 0
                    continue
                end
                
                % Within times meeting conditions, find one w min speed
                % Old way
%                 [~,idx_spd_min] = min(spd(cond_reach_end)); % Within condition-meeting window, idx of slowest movement
%                 idx_reach_end = find(cond_reach_end==1,1) + idx_spd_min - 1; % reach end bin = first bin meeting cond + bin (within cond bins) w slowest speed
%                 reach_end = ts(idx_reach_end);
                
                % New way
                cond_reach_end_nan = 1.*cond_reach_end; % initialize and convert to scalar array
                cond_reach_end_nan(cond_reach_end_nan < 1) = nan;
                [~,idx_reach_end] = min(spd.*cond_reach_end_nan);
                reach_end = ts(idx_reach_end);
                
                
            else
                % Find time of min velocity after last reach
                cond_reach_end = logical((spd < spd_thresh).*(ts > (cue_on + buff)).*(ts < (tr_end + end_buff))); % Buffer is there bc velocity is likely 0 when go cue on. Make sure it's after movement started.
                
                % If conditions not met, skip
                if sum(cond_reach_end) == 0
                    continue
                end
                
                % Old way
%                 [~,idx_spd_min] = min(spd(cond_reach_end));
%                 idx_reach_end = find(cond_reach_end==1,1) + idx_spd_min - 1; 
%                 reach_end = ts(idx_reach_end);
                
                % New way
                cond_reach_end_nan = 1.*cond_reach_end; % initialize and convert to scalar array
                cond_reach_end_nan(cond_reach_end_nan < 1) = nan;
                [~,idx_reach_end] = min(spd.*cond_reach_end_nan);
                reach_end = ts(idx_reach_end);

            end

            
            % Define window of time to save
%             wind_trial = logical((ts>=reach_st).*(ts<=reach_end));
            wind_reach = logical((ts>=wind_st).*(ts<=reach_end));
            
            
            % Add meta-data
            Data2.trial_num{idx,1} = tr;
            Data2.reach_num{idx,1} = reach;
            Data2.reach_st{idx,1} = cue_on;
            Data2.cue_on{idx,1} = cue_on;
            Data2.reach_end{idx,1} = reach_end;
            Data2.reach_pos_st{idx,1} = Data.kinematics{tr}(idx_cue_on2,1:2);
            Data2.reach_pos_end{idx,1} = Data.kinematics{tr}(idx_reach_end,1:2);
            delta_pos = Data2.reach_pos_end{idx,1} - Data2.reach_pos_st{idx,1};
            [Data2.reach_dir{idx,1}, Data2.reach_len{idx,1}] = cart2pol(delta_pos(1),delta_pos(2));
            
            idx_target_on = 1 + ceil(abs(time_before_cue)/dt); % Target is on in bin 1 unless extra time before is added
            temp = zeros(sum(wind_reach),1); temp(idx_target_on) = 1;
            Data2.target_on{idx,1} = temp; 

            % Copy stuff
            Data2.kinematics{idx} = Data.kinematics{tr}(wind_reach,:);
%             Data2.neural_data_M1{idx} = Data.neural_data_M1{tr}(:,wind_reach);
            Data2.neural_data_PMd{idx} = Data.neural_data_PMd{tr}(:,wind_reach);
            Data2.block_info = Data.block_info;
            Data2.time_window{idx,1} = wind_reach;
            Data2.timestamps{idx,1} = Data.timestamps{tr}(wind_reach);

            if (Data2.reach_len{idx} < min_reach_len) || (sum(wind_reach) > max_reach_time)
                continue
            end

            idx = idx + 1;

        end
    
    
end

%% Plot some stuff

make_plots = 0;

if make_plots
    block_num = 1;
    num_trials = length(Data(block_num).neural_data_M1);
    
    for tr_num = 201:300
        figure
        plot(Data(1).kinematics{tr_num}(:,1),Data(1).kinematics{tr_num}(:,2)); %hold on
        axis([-12 12 -12 12])
        axis square
        %     figure
        %     plot(Data(1).kinematics{tr_num}(:,5),Data(1).kinematics{tr_num}(:,6)); %hold on
        %     axis([-20 20 -20 20])
        %     axis square
        
        %     figure
        %     subplot(6,1,1)
        %     plot(Data(1).kinematics{tr_num}(:,1))
        %     subplot(6,1,2)
        %     plot(Data(1).kinematics{tr_num}(:,2))
        %     subplot(6,1,3)
        %     plot(Data(1).kinematics{tr_num}(:,3))
        %     subplot(6,1,4)
        %     plot(Data(1).kinematics{tr_num}(:,4))
        %     subplot(6,1,5)
        %     plot(Data(1).kinematics{tr_num}(:,5))
        %     subplot(6,1,6)
        %     plot(Data(1).kinematics{tr_num}(:,6))
    end
    %hold off
end

