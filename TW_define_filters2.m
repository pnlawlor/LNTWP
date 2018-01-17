function [ filt_struct ] = TW_define_filters2( num_nrn, num_bf, plot_filt, dt )
% Defines how to filter the covariates for the Associative Learning
% analysis projects
% Order of original X: neural data, CS on, CS off, US on, US off, blink
% amplitude
% NOTE: All original covariates will be removed. In order to keep an
% original covariate as-is, explicitly filter it with a 1

if nargin < 3
    plot_filt = false;
    dt = .01;
    disp(['Using 10 ms time bins by default'])
end

%% Initialization

idx = 1;

%% Assign movement filters

time_window = 0.2; % 200 ms; temporal extent of spike history terms

% Get basis
rcos_basis = getBasis('rcos',num_bf,round(time_window/dt))'; % Ian's function
% [~,rcos_basis] = rcosbasis(1:80,[20 20 20 20],[3 5 7.5 10]); % 1 = 200ms total width, 2 = 300, 3 = 400ms;

% Loop through for each temporal basis function
for i = 1:num_bf
    filt_struct(idx).filt_name = ['Spike history ' num2str(i)];
    filt_struct(idx).cov_num_in = [1:num_nrn];
    filt_struct(idx).filt_func = rcos_basis(:,i);
    filt_struct(idx).filt_shift = 0;
    idx = idx + 1;
end


%% Plot filters

if plot_filt
    for filt_num = 1:length(filt_struct)
        figure()
        plot(filt_struct(filt_num).filt_func)
        title(filt_struct(filt_num).filt_name)
    end
end

end

