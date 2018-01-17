function [ filt_struct ] = TW_define_filters( plot_filt, model_num, dt )
% Defines how to filter the covariates for the Associative Learning
% analysis projects
% Order of original X: neural data, CS on, CS off, US on, US off, blink
% amplitude
% NOTE: All original covariates will be removed. In order to keep an
% original covariate as-is, explicitly filter it with a 1

if nargin < 3
    plot_filt = false;
    dt = .01;
    disp(['Using bin size of 10 ms'])
end

%% Initialization

idx = 1;

%% Assign movement filters

% Get basis
% rcos_basis = getBasis('rcos',5,20)';

% rcos_center = 20;
rcos_center = round(0.2/dt); % 200 ms
basis_length = round(0.8/dt); % 800 ms
width1 = round(.03/dt); % 30 ms
width2 = round(.05/dt); % 50 ms
width3 = round(.075/dt); % 75 ms
width4 = round(.1/dt); % 100 ms

if model_num == 3 || model_num == 9
%     [~,rcos_basis] = rcosbasis(1:80,[20 20 20 20],[3 5 7.5 10]); % 1 = 200ms total width, 2 = 300, 3 = 400ms;
    [~,rcos_basis] = rcosbasis(1:basis_length,[rcos_center rcos_center rcos_center rcos_center],[width1 width2 width3 width4]); 
    rcos_basis = round(rcos_basis); % To make into steps
else
    [~,rcos_basis] = rcosbasis(1:basis_length,[rcos_center rcos_center rcos_center rcos_center],[width1 width2 width3 width4]); 
%     [~,rcos_basis] = rcosbasis(1:80,[20 20 20 20],[3 5 7.5 10]); % 1 = 200ms total width, 2 = 300, 3 = 400ms;
end

if model_num == 0 % M1 movement tuning; not done yet
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Early wide';
    filt_struct(idx).cov_num_in = [1 2];
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = -40;
    filt_struct(idx).filt_shift = -round(0.4/dt); % 400 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Mid wide';
    filt_struct(idx).cov_num_in = [1 2];
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = -20;
    filt_struct(idx).filt_shift = -round(0.2/dt); % 200 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Late wide';
    filt_struct(idx).cov_num_in = [1 2];
    filt_struct(idx).filt_func = rcos_basis(:,3);
    filt_struct(idx).filt_shift = 0;
    idx = idx + 1;
    
elseif model_num == 1 || model_num == 7   % PMd visual response with tuning
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 13';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -20;
    filt_struct(idx).filt_shift = -round(0.2/dt); % -200 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 12';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -15;
    filt_struct(idx).filt_shift = -round(0.15/dt); % -150 ms
    idx = idx + 1; 
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 1';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -10;
    filt_struct(idx).filt_shift = -round(0.1/dt); % -100 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 2';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -5;
    filt_struct(idx).filt_shift = -round(0.05/dt); % -50 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 3';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
    filt_struct(idx).filt_shift = 0;
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 4';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = 5;
    filt_struct(idx).filt_shift = round(0.05/dt); % 50 ms
    idx = idx + 1; 
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 5';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = 10;
    filt_struct(idx).filt_shift = round(0.1/dt); % 100 ms
    idx = idx + 1;    
    
        % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 6';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = 15;
    filt_struct(idx).filt_shift = round(0.15/dt); % 150 ms
    idx = idx + 1;    
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 7';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = 20;
    filt_struct(idx).filt_shift = round(0.2/dt); % 200 ms
    idx = idx + 1;    
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 8';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
    filt_struct(idx).filt_shift = 25;
    filt_struct(idx).filt_shift = round(0.25/dt); % 250 ms
    idx = idx + 1;    
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 9';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = 30;
    filt_struct(idx).filt_shift = round(0.3/dt); % 300 ms
    idx = idx + 1;   
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Early narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -15;
    filt_struct(idx).filt_shift = -round(0.15/dt); % -150 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Mid narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -5;
    filt_struct(idx).filt_shift = -round(0.05/dt); % -50 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Late narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = 5;
    filt_struct(idx).filt_shift = round(0.05/dt); % 50 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very late narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = 15;
    filt_struct(idx).filt_shift = round(0.15/dt); % 150 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very late narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = 25;
    filt_struct(idx).filt_shift = round(0.25/dt); % 250 ms	
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Early medium';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = -10;
    filt_struct(idx).filt_shift = -round(0.1/dt); % -100 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Mid medium';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
    filt_struct(idx).filt_shift = 0;
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Late medium';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = 10;
    filt_struct(idx).filt_shift = round(0.1/dt); % 100 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very late medium';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = 20;
    filt_struct(idx).filt_shift = round(0.2/dt); % 200 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very late medium';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = 300;
    filt_struct(idx).filt_shift = round(0.3/dt); % 300 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Early wide';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
%     filt_struct(idx).filt_shift = -10;
    filt_struct(idx).filt_shift = -round(0.1/dt); % -100 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Mid wide';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
    filt_struct(idx).filt_shift = 0;
    idx = idx + 1; 
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Late wide';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
%     filt_struct(idx).filt_shift = 10;
    filt_struct(idx).filt_shift = round(0.1/dt); % 100 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very late wide';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
%     filt_struct(idx).filt_shift = 20;
    filt_struct(idx).filt_shift = round(0.2/dt); % 200 ms
    idx = idx + 1;
    
elseif model_num == 2 % PMd visual response without tuning
    
        % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 13';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -20;
    filt_struct(idx).filt_shift = -round(0.2/dt); % -200 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 12';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -15;
    filt_struct(idx).filt_shift = -round(0.15/dt); % -150 ms
    idx = idx + 1; 
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 1';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -10;
    filt_struct(idx).filt_shift = -round(0.1/dt); % -100 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 2';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -5;
    filt_struct(idx).filt_shift = -round(0.05/dt); % -50 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 3';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
    filt_struct(idx).filt_shift = 0;
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 4';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = 5;
    filt_struct(idx).filt_shift = round(0.05/dt); % 50 ms
    idx = idx + 1; 
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 5';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = 10;
    filt_struct(idx).filt_shift = round(0.1/dt); % 100 ms
    idx = idx + 1;    
    
        % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 6';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = 15;
    filt_struct(idx).filt_shift = round(0.15/dt); % 150 ms
    idx = idx + 1;    
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 7';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = 20;
    filt_struct(idx).filt_shift = round(0.2/dt); % 200 ms
    idx = idx + 1;    
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 8';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = 25;
    filt_struct(idx).filt_shift = round(0.25/dt); % 250 ms
    idx = idx + 1;    
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 9';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = 30;
    filt_struct(idx).filt_shift = round(0.3/dt); % 300 ms
    idx = idx + 1;   
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Early narrow';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -15;
    filt_struct(idx).filt_shift = -round(0.15/dt); % -150 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Mid narrow';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -5;
    filt_struct(idx).filt_shift = -round(0.05/dt); % -50 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Late narrow';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = 5;
    filt_struct(idx).filt_shift = round(0.05/dt); % 50 ms	
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very late narrow';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = 15;
    filt_struct(idx).filt_shift = round(0.15/dt); % 150 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very late narrow';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = 25;
    filt_struct(idx).filt_shift = round(0.25/dt); % 250 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Early medium';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = -10;
    filt_struct(idx).filt_shift = -round(0.1/dt); % -100 ms	
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Mid medium';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
    filt_struct(idx).filt_shift = 0;
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Late medium';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = 10;
    filt_struct(idx).filt_shift = round(0.1/dt); % 100 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very late medium';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = 20;
    filt_struct(idx).filt_shift = round(0.2/dt); % 200 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very late medium';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = 30;
    filt_struct(idx).filt_shift = round(0.3/dt); % 300 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Early wide';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
%     filt_struct(idx).filt_shift = -10;
    filt_struct(idx).filt_shift = -round(0.1/dt); % -100 ms	
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Mid wide';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
    filt_struct(idx).filt_shift = 0;
    idx = idx + 1; 
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Late wide';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
%     filt_struct(idx).filt_shift = 10;
    filt_struct(idx).filt_shift = round(0.1/dt); % 100 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very late wide';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
%     filt_struct(idx).filt_shift = 20;
    filt_struct(idx).filt_shift = round(0.2/dt); % 200 ms	
    idx = idx + 1;
    
elseif model_num == 3

    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 13';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -20;
    filt_struct(idx).filt_shift = -round(0.2/dt); % -200 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 12';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -15;
    filt_struct(idx).filt_shift = -round(0.15/dt); % -150 ms
    idx = idx + 1; 
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 1';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -10;
    filt_struct(idx).filt_shift = -round(0.1/dt); % -100 ms	
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 2';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -5;
    filt_struct(idx).filt_shift = -round(0.05/dt); % -50 ms	
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 3';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
    filt_struct(idx).filt_shift = 0;
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 4';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = 5;
    filt_struct(idx).filt_shift = round(0.05/dt); % 50 ms
    idx = idx + 1; 
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 5';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = 10;
    filt_struct(idx).filt_shift = round(0.1/dt); % 100 ms
    idx = idx + 1;    
    
        % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 6';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = 15;
    filt_struct(idx).filt_shift = round(0.15/dt); % 150 ms
    idx = idx + 1;    
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 7';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = 20;
    filt_struct(idx).filt_shift = round(0.2/dt); % 200 ms	
    idx = idx + 1;    
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 8';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = 25;
    filt_struct(idx).filt_shift = round(0.25/dt); % 250 ms
    idx = idx + 1;    
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow 9';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = 30;
    filt_struct(idx).filt_shift = round(0.3/dt); % 300 ms
    idx = idx + 1;   
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Early narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -15;
    filt_struct(idx).filt_shift = -round(0.15/dt); % -150 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Mid narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -5;
    filt_struct(idx).filt_shift = -round(0.05/dt); % -50 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Late narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = 5;
    filt_struct(idx).filt_shift = round(0.05/dt); % 50 ms	
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very late narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = 15;
    filt_struct(idx).filt_shift = round(0.15/dt); % 150 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very late narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = 25;
    filt_struct(idx).filt_shift = round(0.25/dt); % 250 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Early medium';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = -10;
    filt_struct(idx).filt_shift = -round(0.1/dt); % -100 ms	
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Mid medium';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
    filt_struct(idx).filt_shift = 0;
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Late medium';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = 10;
    filt_struct(idx).filt_shift = round(0.1/dt); % 100 ms	
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very late medium';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = 20;
    filt_struct(idx).filt_shift = round(0.2/dt); % 200 ms	
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very late medium';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = 30;
    filt_struct(idx).filt_shift = round(0.3/dt); % 300 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Early wide';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
%     filt_struct(idx).filt_shift = -10;
    filt_struct(idx).filt_shift = -round(0.1/dt); % -100 ms	
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Mid wide';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
    filt_struct(idx).filt_shift = 0;
    idx = idx + 1; 
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Late wide';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
%     filt_struct(idx).filt_shift = 10;
    filt_struct(idx).filt_shift = round(0.1/dt); % 100 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very late wide';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
%     filt_struct(idx).filt_shift = 20;
    filt_struct(idx).filt_shift = round(0.2/dt); % 200 ms
    idx = idx + 1;    

elseif model_num == 4 || model_num == 5 || model_num == 6 
    % Movement filter 1
    filt_struct(idx).filt_name = 'Early';
    filt_struct(idx).cov_num_in = [3 4 5 6]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -50; % -50 shifts it back by ~30
    filt_struct(idx).filt_shift = -round(0.5/dt); % -50 shifts it back by ~30
    idx = idx + 1;  
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Early';
    filt_struct(idx).cov_num_in = [3 4 5 6]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -45; % -50 shifts it back by ~30
    filt_struct(idx).filt_shift = -round(0.45/dt); % -50 shifts it back by ~30
    idx = idx + 1;  
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Early';
    filt_struct(idx).cov_num_in = [3 4 5 6]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -40; % -50 shifts it back by ~30
    filt_struct(idx).filt_shift = -round(0.4/dt); % -50 shifts it back by ~30
    idx = idx + 1;  
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Early';
    filt_struct(idx).cov_num_in = [3 4 5 6]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
    filt_struct(idx).filt_shift = -35; % -50 shifts it back by ~30
    filt_struct(idx).filt_shift = -round(0.35/dt); % -50 shifts it back by ~30
    idx = idx + 1;      
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Early';
    filt_struct(idx).cov_num_in = [3 4 5 6]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -30; % -50 shifts it back by ~30
    filt_struct(idx).filt_shift = -round(0.3/dt); % -300 ms
    idx = idx + 1;  
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Early';
    filt_struct(idx).cov_num_in = [3 4 5 6]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -25; % -50 shifts it back by ~30
    filt_struct(idx).filt_shift = -round(0.25/dt); % -250 ms
    idx = idx + 1;  

% ===== Model 8 =====
elseif model_num == 8
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -50;
    filt_struct(idx).filt_shift = -round(0.5/dt); % -500 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -45;
    filt_struct(idx).filt_shift = -round(0.45/dt); % -450 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -40;
    filt_struct(idx).filt_shift = -round(0.4/dt); % -400 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -35;
    filt_struct(idx).filt_shift = -round(0.35/dt); % -350 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -30;
    filt_struct(idx).filt_shift = -round(0.3/dt); % -300 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -25;
    filt_struct(idx).filt_shift = -round(0.25/dt); % -250 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -20;
    filt_struct(idx).filt_shift = -round(0.2/dt); % -200 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -15;
    filt_struct(idx).filt_shift = -round(0.15/dt); % -150 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -10;
    filt_struct(idx).filt_shift = -round(0.1/dt); % -100 ms
    idx = idx + 1;

    % Movement filter 1
    filt_struct(idx).filt_name = 'Narrow';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -50;
    filt_struct(idx).filt_shift = -round(0.5/dt); % -500 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Narrow';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -40;
    filt_struct(idx).filt_shift = -round(0.4/dt); % -400 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Narrow';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -30;
    filt_struct(idx).filt_shift = -round(0.3/dt); % -300 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Narrow';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -20;
    filt_struct(idx).filt_shift = -round(0.2/dt); % -200 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Narrow';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -10;
    filt_struct(idx).filt_shift = -round(0.1/dt); % -100 ms	
    idx = idx + 1;

    % Movement filter 1
    filt_struct(idx).filt_name = 'Medium';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = -50;
    filt_struct(idx).filt_shift = -round(0.5/dt); % -500 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Medium';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = -35;
    filt_struct(idx).filt_shift = -round(0.35/dt); % -350 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Medium';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = -20;
    filt_struct(idx).filt_shift = -round(0.2/dt); % -200 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Wide';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
%     filt_struct(idx).filt_shift = -50;
    filt_struct(idx).filt_shift = -round(0.5/dt); % -500 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Wide';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
%     filt_struct(idx).filt_shift = -30;
    filt_struct(idx).filt_shift = -round(0.3/dt); % -300 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Wide';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
%     filt_struct(idx).filt_shift = -10;
    filt_struct(idx).filt_shift = -round(0.1/dt); % -100 ms	
    idx = idx + 1;

% ===== Model 9 =====    
elseif model_num == 9 
%     % Movement filter 1
%     filt_struct(idx).filt_name = 'Very narrow';
%     filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
% %     filt_struct(idx).cov_num_in = [1]; % sin and cos
%     filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -45;
%     idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -46;
    filt_struct(idx).filt_shift = -round(0.46/dt); % 460 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -39;
    filt_struct(idx).filt_shift = -round(0.39/dt); % 390 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -32;
    filt_struct(idx).filt_shift = -round(0.32/dt); % 320 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -25;
    filt_struct(idx).filt_shift = -round(0.25/dt); % -250 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -18;
    filt_struct(idx).filt_shift = -round(0.18/dt);
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -11;
    filt_struct(idx).filt_shift = -round(0.11/dt);
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Narrow';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -41;
    filt_struct(idx).filt_shift = -round(0.41/dt);
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Narrow';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -30;
    filt_struct(idx).filt_shift = -round(0.3/dt); % -300 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Narrow';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -19;
    filt_struct(idx).filt_shift = -round(0.19/dt);
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Medium';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = -40;
    filt_struct(idx).filt_shift = -round(0.4/dt); % -400 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Medium';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = -25;
    filt_struct(idx).filt_shift = -round(0.25/dt); % -250 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Wide';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
%     filt_struct(idx).filt_shift = -40;
    filt_struct(idx).filt_shift = -round(0.4/dt); % -400 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Wide';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
%     filt_struct(idx).filt_shift = -19;
    filt_struct(idx).filt_shift = -round(0.19/dt);
    idx = idx + 1;
    
   
elseif model_num == 10 % ===== Model 10: Untuned =====
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1]; % 1 = untuned, 2 & 3 = sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -50;
    filt_struct(idx).filt_shift = -round(0.5/dt); % -500 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -45;
    filt_struct(idx).filt_shift = -round(0.45/dt); % -450 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -40;
    filt_struct(idx).filt_shift = -round(0.4/dt); % -400 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -35;
    filt_struct(idx).filt_shift = -round(0.35/dt); % -350 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -30;
    filt_struct(idx).filt_shift = -round(0.3/dt); % -300 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -25;
    filt_struct(idx).filt_shift = -round(0.25/dt); % -250 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -20;
    filt_struct(idx).filt_shift = -round(0.2/dt); % -200 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -15;
    filt_struct(idx).filt_shift = -round(0.15/dt); % -150 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -10;
    filt_struct(idx).filt_shift = -round(0.1/dt); % -100 ms	
    idx = idx + 1;

    % Movement filter 1
    filt_struct(idx).filt_name = 'Narrow';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -50;
    filt_struct(idx).filt_shift = -round(0.5/dt); % -500 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Narrow';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -40;
    filt_struct(idx).filt_shift = -round(0.4/dt); % -400 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Narrow';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -30;
    filt_struct(idx).filt_shift = -round(0.3/dt); % -300 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Narrow';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -20;
    filt_struct(idx).filt_shift = -round(0.2/dt); % -200 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Narrow';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -10;
    filt_struct(idx).filt_shift = -round(0.1/dt); % -100 ms		
    idx = idx + 1;

    % Movement filter 1
    filt_struct(idx).filt_name = 'Medium';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = -50;
    filt_struct(idx).filt_shift = -round(0.5/dt); % -500 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Medium';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = -35;
    filt_struct(idx).filt_shift = -round(0.35/dt); % -350 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Medium';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = -20;
    filt_struct(idx).filt_shift = -round(0.2/dt); % -200 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Wide';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
%     filt_struct(idx).filt_shift = -50;
    filt_struct(idx).filt_shift = -round(0.5/dt); % -500 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Wide';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
%     filt_struct(idx).filt_shift = -30;
    filt_struct(idx).filt_shift = -round(0.3/dt); % -300 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Wide';
    filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
%     filt_struct(idx).filt_shift = -10;
    filt_struct(idx).filt_shift = -round(0.1/dt); % -100 ms	
    idx = idx + 1;
    
    % Kinematics filter 1
    filt_struct(idx).filt_name = 'Kinematics: filtered with narrow rcos';
    filt_struct(idx).cov_num_in = [8];
    filt_struct(idx).filt_func = rcos_basis(:,1);
    filt_struct(idx).filt_shift = 0 - rcos_center;
    idx = idx + 1;
    
    % Kinematics filter 2
    filt_struct(idx).filt_name = 'Kinematics: filtered with narrow rcos';
    filt_struct(idx).cov_num_in = [8];
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -15 - rcos_center;
    filt_struct(idx).filt_shift = -round(0.15/dt) - rcos_center;
    idx = idx + 1;
    
% ===== Model 11 =====
elseif model_num == 11
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -50;
    filt_struct(idx).filt_shift = -round(0.5/dt); % -500 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -45;
    filt_struct(idx).filt_shift = -round(0.45/dt); % -450 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -40;
    filt_struct(idx).filt_shift = -round(0.4/dt); % -400 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -35;
    filt_struct(idx).filt_shift = -round(0.35/dt); % -350 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -30;
    filt_struct(idx).filt_shift = -round(0.3/dt); % -300 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -25;
    filt_struct(idx).filt_shift = -round(0.25/dt); % -250 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -20;
    filt_struct(idx).filt_shift = -round(0.2/dt); % -200 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -15;
    filt_struct(idx).filt_shift = -round(0.15/dt); % -150 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Very narrow';
    filt_struct(idx).cov_num_in = [1 2 3]; % sin and cos
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -10;
    filt_struct(idx).filt_shift = -round(0.1/dt); % -100 ms	
    idx = idx + 1;

    % Movement filter 1
    filt_struct(idx).filt_name = 'Narrow';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -50;
    filt_struct(idx).filt_shift = -round(0.5/dt); % -500 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Narrow';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -40;
    filt_struct(idx).filt_shift = -round(0.4/dt); % -400 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Narrow';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -30;
    filt_struct(idx).filt_shift = -round(0.3/dt); % -300 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Narrow';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -20;
    filt_struct(idx).filt_shift = -round(0.2/dt); % -200 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Narrow';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,2);
%     filt_struct(idx).filt_shift = -10;
    filt_struct(idx).filt_shift = -round(0.1/dt); % -100 ms	
    idx = idx + 1;

    % Movement filter 1
    filt_struct(idx).filt_name = 'Medium';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = -50;
    filt_struct(idx).filt_shift = -round(0.5/dt); % -500 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Medium';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = -35;
    filt_struct(idx).filt_shift = -round(0.35/dt); % -350 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Medium';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,3);
%     filt_struct(idx).filt_shift = -20;
    filt_struct(idx).filt_shift = -round(0.2/dt); % -200 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Wide';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
%     filt_struct(idx).filt_shift = -50;
    filt_struct(idx).filt_shift = -round(0.5/dt); % -500 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Wide';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
%     filt_struct(idx).filt_shift = -30;
    filt_struct(idx).filt_shift = -round(0.3/dt); % -300 ms
    idx = idx + 1;
    
    % Movement filter 1
    filt_struct(idx).filt_name = 'Wide';
    filt_struct(idx).cov_num_in = [1 2 3];
%     filt_struct(idx).cov_num_in = [1]; % sin and cos
    filt_struct(idx).filt_func = rcos_basis(:,4);
%     filt_struct(idx).filt_shift = -10;
    filt_struct(idx).filt_shift = -round(0.1/dt); % -100 ms	
    idx = idx + 1;
    
    % Kinematics filter 1
    filt_struct(idx).filt_name = 'Kinematics: filtered with narrow rcos';
    filt_struct(idx).cov_num_in = [4 5 6 7 8];
    filt_struct(idx).filt_func = rcos_basis(:,1);
    filt_struct(idx).filt_shift = 0 - rcos_center;
    idx = idx + 1;
    
    % Kinematics filter 2
    filt_struct(idx).filt_name = 'Kinematics: filtered with narrow rcos';
    filt_struct(idx).cov_num_in = [4 5 6 7 8];
    filt_struct(idx).filt_func = rcos_basis(:,1);
%     filt_struct(idx).filt_shift = -15 - rcos_center;
    filt_struct(idx).filt_shift = -round(0.15/dt) - rcos_center;
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

