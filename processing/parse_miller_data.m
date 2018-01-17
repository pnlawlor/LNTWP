function [ neural_data, kine_data ] = parse_miller_data( data, kinematics )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% tt has trial info. 6th column is go cue, 7th is trial end. still need
% kine data from pavan. need movement start for each trial 

dt = .01;
num_bins_before = 50;
num_bins_after = 150;

num_trials = ; % something
num_units = length(data);

neural_data = cell(num_trials,num_units);
% neural_data_temp = cell(num_units,1);
move_start = nan(num_trials,1);
move_dir = nan(num_trials,1);
kine_data = cell(num_trials,1);

    
neural_data_temp = getSpkMat(data,dt);
    

for trial_num = 1:num_trials
    move_start = ;
    move_start_bin = round(move_start);
    move_dir = ;
    
    kine_data{trial_num} = kinematics((move_start_bin-num_bins_before):(move_start_bin+num_bins_after),:);
    
    for unit_num = 1:num_units
        neural_data{trial_num,unit_num} = ... 
            neural_data_temp{trial_num}((move_start_bin-num_bins_before):(move_start_bin+num_bins_after));
    end
end


end

