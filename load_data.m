if model_num == 8 || model_num == 10 || model_num == 11
    if real_data
        if ~quest
            data_path = './source_data/processed/';
        end

            % Data files - basic
            data_fname = 'thresholds/MM_S1_8cms.mat'; % Monkey M, session 1. 10 ms time bins. 8 cm/s velocity threshold for defining reach start.
%             data_fname = 'thresholds/MT_S1_8cms.mat'; % Monkey T, session 1. 10 ms time bins. 8 cm/s velocity threshold for defining reach start. 
%             data_fname = 'thresholds/MT_S2_8cms.mat'; % Monkey T, session 2. 10 ms time bins. 8 cm/s velocity threshold for defining reach start.
%             data_fname = 'thresholds/MT_S3_8cms.mat'; % Monkey T, session 3. 10 ms time bins. 8 cm/s velocity threshold for defining reach start.

            % Data files - different time bin sizes
%             data_fname = 'bin_sizes\MM_S1_5ms.mat'; % Monkey M, session 1. 5 ms time bins. 8 cm/s velocity threshold for defining reach start.
%             data_fname = 'bin_sizes\MM_S1_10ms.mat'; % Monkey M, session 1. 10 ms time bins. 8 cm/s velocity threshold for defining reach start.
%             data_fname = 'bin_sizes\MM_S1_20ms.mat'; % Monkey M, session 2. 20 ms time bins. 8 cm/s velocity threshold for defining reach start.
%             data_fname = 'bin_sizes\MM_S1_40ms.mat'; % Monkey M, session 3. 40 ms time bins. 8 cm/s velocity threshold for defining reach start.

            % Data files - different velocity thresholds
%             data_fname = 'thresholds\MM_S1_6cms'; % Monkey M, session 1. 10 ms time bins. 6 cm/s velocity threshold for defining reach start.
%             data_fname = 'thresholds\MM_S1_8cms'; % Monkey M, session 1. 10 ms time bins. 8 cm/s velocity threshold for defining reach start.
%             data_fname = 'thresholds\MM_S1_10cms'; % Monkey M, session 1. 10 ms time bins. 10 cm/s velocity threshold for defining reach start.
%             data_fname = 'thresholds\MM_S1_12cms'; % Monkey M, session 1. 10 ms time bins. 12 cm/s velocity threshold for defining reach start.
        
        load([data_path data_fname])
        disp(['Loaded: ' data_fname])
        
        spikes_all_preselect = Data_win.spikes_PMd; % PMd
        spikes_all_preselect = remove_duplicate_neurons(spikes_all_preselect,0.4); % Removes highly correlated units. Default correlation threshold is 0.4. 
        
        num_nrn = size(spikes_all_preselect,2);
        
        outer_target_on = Data_win.outer_target;
        kinematics = Data_win.kinematics;
    else
        
        if ~quest
            data_path = '.\simulation_data\';
        else
%             data_path = '/projects/b1024/pat/TW/Source_data/'; % Quest
        end
        

        if strcmp(warp_type,'none')
            data_fname = 'SimData_250nrn_400tr-11-17-16--08-29.mat'; % SimData; 250nrns, 400tr, 010 warp, run0
        elseif strcmp(warp_type,'weak')
            data_fname = 'SimData_250nrn_400tr-11-17-16--08-32.mat'; % SimData; 250nrns, 400tr, 111 warp, run0
        elseif strcmp(warp_type,'strong')
            data_fname = 'SimData_250nrn_400tr-11-17-16--08-34.mat'; % SimData; 250nrns, 400tr, 111 warp, run3
        end
        
        load([data_path data_fname])
        disp(['Loaded: ' data_fname])
        
        bins_per_trial = size(Data_win.spikes_PMd{1,1},1);
        num_nrn = size(Data_win.spikes_PMd,2);
        num_trials = size(Data_win.spikes_PMd,1);
        spikes_all_preselect = Data_win.spikes_PMd; % PMd (simulated)
        
    end

end

max_num_trials = size(spikes_all_preselect,1);
movement_data_all = Data_win.movement_cov;