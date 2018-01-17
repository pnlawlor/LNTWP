function [ predictions, fit_parameters, fit_info, pseudo_R2] = fit_poiss_GLM( X_data, y_data, num_CV, dt, reg_strength, alpha, type, bins_per_trial, num_BS )
% Fit Poisson GLM

%% Initialize

if nargin < 7
    bins_per_trial = 700;
    alpha = 1;
    disp('Using L1 regularization')
    num_BS = 1000;
elseif nargin < 6
    % Set regularization strength
    alpha = 1;
    bins_per_trial = 700;
    disp('Using L1 regularization')
    num_BS = 1000;
elseif nargin < 5
    reg_strength = false;
    bins_per_trial = 700;
    disp('Fitting regularization parameter')
    num_BS = 1000;
elseif nargin < 4
    dt = .01;
    bins_per_trial = 700;
    disp('dt not supplied, using dt = .01 (10 ms)')
    num_BS = 1000;
elseif nargin < 3
    num_CV = 10;
    bins_per_trial = 700;
    disp('Using 10 cross validation folds')
    num_BS = 1000;
elseif nargin < 2
    bins_per_trial = 700;
    error('You must provide X and y data')
end

fit_info = cell(1,num_CV);
fit_parameters = cell(1,num_CV);
predictions = repmat({nan(1,length(y_data))},num_CV,1); % One cell for each CV, each cell initialized with NaN
X_data = full(X_data);

% Whether to calculate pseudoR2 for each fold or over all test data
old_pr2 = 0;

if old_pr2
    pseudo_R2 = nan(num_CV,3);
else
    pseudo_R2 = nan(1,3);
end

%% User inputs

old_format = (length(bins_per_trial)==1);

if old_format
    num_trials = length(y_data)/bins_per_trial;
    trial_numbers = reshape(repmat((1:num_trials),[bins_per_trial 1]),[],1);
else
    num_trials = length(bins_per_trial);
    trial_numbers = [];
    for tr = 1:num_trials
        trial_numbers = [trial_numbers; repmat(tr,bins_per_trial(tr),1)];
    end
end

use_par_bootci = 0;
num_lambda = 200;
num_internal_CV = 10;

% For Poisson
link_f = 'log';
prob_family = 'poisson';

% For Bernoulli
% link_f = 'logit';
% prob_family = 'binomial';

if strcmp(prob_family,'binomial')
    y_data = logical(y_data>0);
end

if use_par_bootci
    bootci_options = statset('UseParallel',true);
else
    bootci_options = statset('UseParallel',false);
end

%% Initialize cross-validation

% cv_partition = cvpartition(y_data,'k',num_CV);      % Data structure that keeps track of training/testing indices
cv_partition = cvpartition_trialwise(trial_numbers,num_CV); % My version of cvpartition which allocates partitions by trials

% If glmnet, initialize some stuff
if strmatch(type,'glmnet')
    opt.alpha = alpha;
    opt.intr = true;
    
    % Auto determination of lambda or not
    if any(reg_strength == false)
        opt.nlambda = num_lambda;
        opt.lambda = logspace(-2,0,num_lambda);
    else
        opt.lambda = reg_strength;
    end
    options = glmnetSet(opt);
end

%% Fit

for CV = 1:num_CV
    % Get indices
%     idx_train = cv_partition.training(CV);
%     idx_test = cv_partition.test(CV);
    idx_train = cv_partition(CV).training;
    idx_test = cv_partition(CV).test;
    
    % Fit
    switch type
        case 'lassoglm'
            % Fit using LassoGLM; cross validate to find optimal penalty param
            
            if (CV == 1) && (any(reg_strength == false) || length(reg_strength)>1)
                [fit_temp, fit_info_temp] = lassoglm(...
                    X_data(:,:),y_data(:),...
                    prob_family, ...           % Link function
                    'offset', repmat(log(dt),size(y_data(:))), ...
                    'CV', 10, ...             % Number of CVs to determine regularization parameter
                    'Alpha', alpha, ...          % Type of regularization (L2 = small; L1 = 1)
                    'Lambda', reg_strength ...
                    );

                % Append intercept
                idx_min_dev = fit_info_temp.IndexMinDeviance;
                reg_strength = fit_info_temp.Lambda(idx_min_dev);
            end
            
            [fit_temp, fit_info{CV}] = lassoglm(...
                    X_data(idx_train,:),y_data(idx_train),...
                    prob_family, ...           % Link function
                    'offset', repmat(log(dt),size(y_data(idx_train))), ...
                    'Alpha', alpha, ...          % Type of regularization (L2 = small; L1 = 1)
                    'Lambda', reg_strength ...
                    );
                
            fit_parameters{CV} = [fit_info{CV}.Intercept; fit_temp];
        
        case 'glmfit'
            % Fit using GLMFit, no regularization
            [fit_parameters{CV}, ~, fit_info{CV}] = glmfit(...
                                                         X_data(idx_train,:), y_data(idx_train), ...
                                                         prob_family, ...
                                                         'offset', repmat(log(dt),size(y_data(idx_train))) ...
                                                         );
                                                     
        case 'glmnet'
            % Set size of offset correctly
            opt.offset = repmat(log(dt),size(y_data(idx_train)));
            
            try
                % If first CV fold, find best lambda if not supplied
                if (CV == 1) && (any(reg_strength == false) || length(reg_strength)>1)
                    % Cross validate (automatically) to find best lambda                

                    if any(reg_strength) == false
                        opt.nlambda = num_lambda;
    %                   opt.lambda = logspace(-3,0,num_lambda);
    %                   opt.lambda = [];
                        options = glmnetSet(opt);
                    else
                        opt.lambda = reg_strength;
                        options = glmnetSet(opt);
                    end

    %                 tic
                    % Use only training data to find best lambda
    %                 cv_fit = cvglmnet(...
    %                                     X_data(idx_train), y_data(idx_train), ...
    %                                     prob_family, options,[],num_internal_CV);
    %                                 
                    % Use all data to find best lambda
                    cv_fit = cvglmnet(...
                                        X_data(idx_train,:), y_data(idx_train), ...
                                        prob_family, options,[],num_internal_CV);

    %                 toc
                    best_lambda = cv_fit.lambda_min; % lambda of min deviance
    %                 best_lambda = cv_fit.lambda_1se; % lambda within 1se of min deviance lambda. to avoid overfitting, more conservative
                    opt.lambda = best_lambda;
                    opt.nlambda = [];
                end

                % Set options (single lambda)
                options = glmnetSet(opt);

                % Then do regular train/test fitting
                fit_info{CV} = glmnet(X_data(idx_train,:),y_data(idx_train),prob_family,options);
                fit_parameters{CV} = glmnetCoef(fit_info{CV});
            
            % If glmnet doesn't work for some reason
            catch ee
                
                disp('glmnet broke. Using lassoglm')
                
                if (CV == 1) && (any(reg_strength == false) || length(reg_strength)>1)
                [fit_temp, fit_info_temp] = lassoglm(...
                    X_data(:,:),y_data(:),...
                    prob_family, ...           % Link function
                    'offset', repmat(log(dt),size(y_data(:))), ...
                    'CV', 10, ...             % Number of CVs to determine regularization parameter
                    'Alpha', alpha, ...          % Type of regularization (L2 = small; L1 = 1)
                    'Lambda', reg_strength ...
                    );

                % Append intercept
                idx_min_dev = fit_info_temp.IndexMinDeviance;
                reg_strength = fit_info_temp.Lambda(idx_min_dev);
            end
            
            [fit_temp, fit_info{CV}] = lassoglm(...
                    X_data(idx_train,:),y_data(idx_train),...
                    prob_family, ...           % Link function
                    'offset', repmat(log(dt),size(y_data(idx_train))), ...
                    'Alpha', alpha, ...          % Type of regularization (L2 = small; L1 = 1)
                    'Lambda', reg_strength ...
                    );
                
            fit_parameters{CV} = [fit_info{CV}.Intercept; fit_temp];
                
            end
    end
    
    
    % Predict on test data

    predictions{CV}(idx_test) = glmval(...
        fit_parameters{CV},...
        X_data(idx_test,:),...
        link_f,...
        'constant','on',...
        'offset', repmat(log(dt),size(y_data(idx_test))) ...
        );
    
    if old_pr2
    % Calculate pseudo R2
    pseudo_R2(CV,1) = compute_pseudo_R2(...
                                    y_data(idx_test),...        % actual y values
                                    predictions{CV}(idx_test)',...         % y predictions
                                    mean(y_data(idx_train))...  % mean of training set y values
                                    );
    
    pseudo_R2(CV,2:3) = pseudo_R2(CV,1)*[1 1];
                                
    disp(['CV ' num2str(CV) ' completed. Pseudo R2 value: ' num2str(pseudo_R2(CV,1)) ' [' num2str(pseudo_R2(CV,2:3)) ']'])
    end
end

if ~old_pr2
    pred_combined = nansum(cell2mat(predictions),1);
    
    pseudo_R2(1,1) = compute_pseudo_R2(...
                                y_data,...        % actual y values
                                pred_combined',...  % y predictions
                                mean(y_data)...  % mean of y values
                                );
    if num_BS == 1
        pseudo_R2(1,2:3) = pseudo_R2(1,1)*[1 1];
    elseif num_BS > 1
        pseudo_R2(1,2:3) = bootci(num_BS, ...
                                {@compute_pseudo_R2, ...
                                y_data, ...
                                pred_combined', ...
                                mean(y_data)}, ...
                                'Options',bootci_options);
    end

    disp(['Pseudo R2 value: ' num2str(pseudo_R2(1,1)) ' [' num2str(pseudo_R2(1,2:3)) ']'])

end

end

