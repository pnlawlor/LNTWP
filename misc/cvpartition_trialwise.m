function [ cv_struct ] = cvpartition_trialwise( trial_numbers, desired_nCV )
% This mimics the functionality of Matlab's cvpartition() which divides
% data into (stratified) CV folds, but ""_trialwise does not split trials into different
% CV folds. This is important for fitting time series because data points
% within a trial are "less indpendent" than data points across trials.
% Thus, one should not train on some data points from a given trial and
% test on the remaining data points from that same trial. 
% This function attempts to create CV folds of the same size without
% splitting trials. This may mean that the CV folds are not exactly the
% same size.
%
% Properties needed:
%   N:              Number of observations
%   NumTestSets:    Number of test sets
%   TestSize:       Size of each test set
%   TrainSize:      Size of each training set
%   Type:           Type of partition
%
% Methods needed: (need to do this?)
%   disp:           Diisplay cvpartition obect
%   display:        Display cvpartition object
%   repartition:    Repartition data for cross-validation
%   test:           Test indices for cross-validation
%   training:       Training indices for cross-validation

%% Initialization

% See if there are same num bins/trial
trial_numbers = reshape(trial_numbers,[],1);
trial_numbers(length(trial_numbers)+1) = max(trial_numbers) + 1; % add one at the end. necessary to check if the last trial is the same length as the others.
[unique_trials, idx]  = unique(sort(trial_numbers));
unique_trials = unique_trials(1:(end-1)); % get rid of the artificially added one
trial_numbers = trial_numbers(1:(end-1)); % same thing
bin_diff = diff(idx);
uniform_trial_size = (numel(unique(bin_diff))==1);

memory_intensive = 0;

%% Divide up data into CV folds

% In the unusual case you don't want to cross-validate - make them the same
if desired_nCV == 1
    cv_struct(1).test = logical(0*trial_numbers + 1);
    cv_struct(1).training = cv_struct(1).test;
    return
end

if uniform_trial_size
    
    trials_per_CV = length(unique_trials)/desired_nCV;
    shuff_trial_order = reshape(unique_trials(randperm(length(unique_trials))),[],1);
    
    tr_idx = 0;
    for CV_num = 1:desired_nCV
        if CV_num == 1
            % Shuffle trials into CV folds
%             CV_tr_nums = shuff_trial_order(1:ceil(trials_per_CV));
            CV_tr_nums = unique_trials(1:desired_nCV:end);
            % Choose bins with the right trial numbers
            cv_struct(CV_num).test = any(repmat(trial_numbers,[1 length(CV_tr_nums)]) == repmat(CV_tr_nums',[length(trial_numbers) 1]),2);
            cv_struct(CV_num).training = ~cv_struct(CV_num).test;
%             tr_idx = tr_idx + length(CV_tr_nums);
        else
            % Shuffle trials into CV folds
%             CV_tr_nums = shuff_trial_order((tr_idx+1):(tr_idx+floor(trials_per_CV)));
            CV_tr_nums = unique_trials(CV_num:desired_nCV:end);
            % Choose bins with the right trial numbers
            cv_struct(CV_num).test = any(repmat(trial_numbers,[1 length(CV_tr_nums)]) == repmat(CV_tr_nums',[length(trial_numbers) 1]),2);
            cv_struct(CV_num).training = ~cv_struct(CV_num).test;
%             tr_idx = tr_idx + length(CV_tr_nums);
        end
    end
    
else
%     error('This has not been coded yet')
    num_trials = length(unique(trial_numbers));
    
    for CV_num = 1:desired_nCV
        
        trials_this_CV = CV_num:desired_nCV:num_trials;
        
        if memory_intensive
            trial_numbers_rep = repmat(trial_numbers,1,length(trials_this_CV));
            trials_this_CV_rep = repmat(trials_this_CV,length(trial_numbers),1);
            cond = any(trial_numbers_rep==trials_this_CV_rep,2);
        else
            cond = zeros(size(trial_numbers,1),1);
            for tr = 1:length(trials_this_CV)
                cond(trial_numbers==trials_this_CV(tr)) = 1;
            end
        end
        
        cv_struct(CV_num).test = logical(cond);
        cv_struct(CV_num).training = logical(~cond);
        
%         cv_struct(CV_num).test = any(repmat(trial_numbers,[1 length(CV_tr_nums)]) == repmat(CV_tr_nums',[length(trial_numbers) 1]),2);
%         cv_struct(CV_num).training = ~cv_struct(CV_num).test;
        
    end
    
end


end

