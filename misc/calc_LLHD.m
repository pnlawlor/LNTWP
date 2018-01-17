function [ LLHD_total, LLHD_mean ] = calc_LLHD( data, preds, opt )

% parse args
if iscell(preds)
    format = 'pat_format'; % cell array format
else
    format = 'array_format';
end

if ~isfield(opt,'type')
    opt.type = 'gaussian';
    disp(['No likelihood function provided. Using Gaussian.'])
end

% initialize
n_vars = size(data,2);
n_dp = size(data,1);
preds_mat = nan(n_dp, n_vars);

% form inputs
if strcmp(format,'pat_format')
    for var = 1:n_vars
        temp = cell2mat(preds{var});
        temp = nansum(temp,1)';
        preds_mat(:,var) = temp;
    end
    
    preds_old = preds;
    preds = preds_mat;
    
end

% check if data and preds are the same size
if size(data) ~= size(preds)
    error('Data and prediction matrices do not have the same size.')
end

% calculate LLHD matrix
switch opt.type
    case 'gaussian'
        error('Not coded yet')
    case 'poisson'
        LLHD_mat = data.*log(preds + eps) - preds - log(factorial(data));
end

% calculate total and mean LLHD
LLHD_total = sum(sum(LLHD_mat));
LLHD_mean = LLHD_total/(n_vars*n_dp);
    
end

