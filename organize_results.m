%% Simulations: Locate simulation results files

sim_results_fpath = '.\simulation_data\results\';

results_folders = {'prior111\','prior121\'};
num_folders = length(results_folders);
files_to_use = 1:15;

%% Simulations: Collect summary statistics from simulation results files for Warp Condition x Num Nrn data


% Create another struct with all of the results
% Struct(i,j,k) gives the ith prior, jth warp condition (no, weak, strong) and kth
% number of neurons (10,20, 50, 100, 200)


for fold_num = 1:num_folders
    
    sim_results_fnames = dir([sim_results_fpath results_folders{fold_num}]);
    sim_results_fnames = sim_results_fnames(3:end); % First 2 are '.' and '..'
    
    for fnum_idx = 1:length(files_to_use)
        fnum = files_to_use(fnum_idx);
        
        load([sim_results_fpath results_folders{fold_num} sim_results_fnames(fnum).name])
        
        % Determine which warp condition this is
        weak = strfind(sim_results_fnames(fnum).name,'weak');
        strong = strfind(sim_results_fnames(fnum).name,'strong');
        
        if weak
            row = 2;
        elseif strong
            row = 3;
        else
            row = 1;
        end
        
        % Determine how many  neurons there are
        num_nrn = Results(1).description.num_nrn;
        
        if num_nrn == 10
            col = 1;
        elseif num_nrn == 20
            col = 2;
        elseif num_nrn == 50
            col = 3;
        elseif num_nrn == 100
            col = 4;
        elseif num_nrn == 200
            col = 5;
        end
        
        SummaryStats(fold_num,row,col).Results = Results;
        SummaryStats(fold_num,row,col).File_name = sim_results_fnames(fnum).name;
        
    end
end

%% Real data: Collect effect sizes for different conditions

% Recent
file_loc = '.\results\conditions\M1\';
% file_loc = '.\results\conditions\M2\';

files = dir(file_loc);
files = files(3:end); % First 2 are . and ..

num_files = length(files);

for fnum = 1:num_files
    
    % Load file
    fname = files(fnum).name;
    load([file_loc fname])
    
    % Stick it in struct
    EffectSizes(fnum).Results = Results;
    EffectSizes(fnum).file_name = fname;
    
end

%% Real data: Collect effect sizes for different velocity thresholds

file_loc = '.\results\thresholds\';

files = dir(file_loc);
files = files(3:end); % First 2 are . and ..

num_files = length(files);

for fnum = 1:num_files
    
    % Load file
    fname = files(fnum).name;
    load([file_loc fname])
    
    % Stick it in struct
    EffectSizes(fnum).Results = Results;
    EffectSizes(fnum).file_name = fname;
    
end

%% Real data: Collect effect size for different time bin sizes

file_loc = '.\results\bin_sizes\';

files = dir(file_loc);
files = files(3:end); % First 2 are . and ..

num_files = length(files);

for fnum = 1:num_files
    
    % Load file
    fname = files(fnum).name;
    load([file_loc fname])
    
    % Stick it in struct
    EffectSizes(fnum).Results = Results;
    EffectSizes(fnum).file_name = fname;
    
end
