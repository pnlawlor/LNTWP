%% User inputs 1

% Detect if running on Quest supercomputer
questdir = '/home/pnl429/';
if isdir(questdir)
    quest = 1;
else
    quest = 0;
end

% Choose model for GLM
model_num = 8; % 8 = full, 10 = untuned + kinematics, 11 = full + kinematics

% Choose properties of data to use
FR_thresh_Hz = 2; 
dt = .01;
FR_thresh = FR_thresh_Hz*dt;

real_data = 0; % 1 = real data from monkey; 0 = simulated data
num_nrn_desired = 10;

%% Load scripts

rng('default')

if ~quest
    addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Scripts\Pat'))
    addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Projects\FEF\Scripts\Pavan_GLM_code_v2\FEF\utilities'))
    addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Scripts_TW'))
else
    addpath(genpath('/projects/b1024/pat/ALFC/Scripts/General/general_scripts/'))
    addpath(genpath('/projects/b1024/pat/TW/Scripts_TW/'))
end
    
%% Choose brain region/model
% 1 = PMd visual with tuning
% 2 = PMd visual without tuning
% 3 = PMd visual with tuning, block basis functions
% 4 = PMd movement with tuning
% 5 = PMd movement with tuning, not tied to kinematics
% 6 = PMd movement with tuning, tied to kinetmatics; Matt data
% 7 = PMd "visual" with tuning, tied to cue onset; Matt data
% 8 = PMd movement with tuning, aligned but not identical to kinematics
% 9 = PMd movement with tuning, aligned but not identical to kinematics, boxcars
% 10 = PMd movement withOUT tuning (i.e., untuned model), aligned but not identical to kinematics
% 11 = PMd movement with rcos bumps AND kinematics. For the haters.

%% Load data

load_data

%% User inputs 2

rng('default')

transition_priors = [0 1; 0 2; 0 1];
transition_prior_scale = 1;

transition_priors = transition_priors/sum(sum(transition_priors));
transition_priors = log(transition_priors);
transition_priors = transition_prior_scale*transition_priors;

fit_method = 'lassoglm'; % Options: glmfit, lassoglm, glmnet
num_CV = 2;
num_WV = 10; 
num_BS_GLM = 1;

tol_type = 'relTol';
conv_thresh = 1e-3;

max_iter = 15; % Max number of EM iterations

alpha = .01; % L1 vs L2; alpha small = L2; alpha 1 = L1
lambda = .05; % Regularization strength

trials = 1:max_num_trials;
neurons = 1:num_nrn;
% neurons = 1:num_nrn_desired; % For simulations

do_shuffle = 0; % Shuffle control: apply learned warps to wrong trials

include_spk_history = 0; % Include spike history terms. This adds an extra iteration. 
num_spk_history_bf = 5; % Number of temporal basis functions to use for spike history

num_rep = 10; % Only necessary for simulated data
Results_fpath = 'C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Data\Temp\'; % home computer
Results_fpath = ''; % quest

% Some initialization
if real_data
    Results_fname_datatype = 'RealData';
else
    Results_fname_datatype = 'SimData';
end

warp_type = 'strongwarp';
Results_fname_base = [Results_fname_datatype '_' num2str(num_nrn_desired) 'nrn_' num2str(length(trials)) 'tr' '_' warp_type];
datetimestr = datestr(now,'mm-dd-yy--HH-MM');

%% Start repetitions (for simulations)

for rep = 1:num_rep

%% Initialize

initialize 

%% Process

process

%% Fit model

fit_model

%% ==== Visualization and Analysis Tools ====

%% Visualize raw data

trs = []; % trs = 1:20;
nrns = []; %idx_sort(1);
iter = 1;

for nrn = nrns
    for tr = trs
        figure
        subplot(2,1,1)
        plot(predictions{tr,nrn,1})
        xlim([1 length(predictions{tr,nrn,iter})])
        
        subplot(2,1,2)
        imagesc(spikes_warped{tr,nrn,iter}')

        title(['Neuron: ' num2str(nrn) ', trial: ' num2str(tr)])
    end
end

%% Visualize covariates and GLM fit

trs = []; %81:90;
nrns = []; %53;
iter = 1;

y_all = cell2mat(spikes);

% Visualize GLM fits
for nrn = nrns
    for tr = trs
        for i = 1:length(iter)
        figure
        subplot(1,3,1)
        imagesc(spikes{tr,nrn},[0,2])
        subplot(1,3,2)
        FR_mean = mean(y_all(:,nrn));
%         imagesc(predictions{tr,nrn,1},[(FR_mean - .03) (FR_mean + .03)])
        imagesc(predictions{tr,nrn,i},[0,1])
        subplot(1,3,3)
        imagesc(X_cell{tr,1})
        end
    end
end

%% Visualize PSTHs and predictions - NOT DONE YET

trs = [];
nrns = [];
iter = 2;

for nrn_idx = 1:length(nrns)
    nrn = nrns(nrn_idx);
    temp_data = nan(length());
    for tr = trs
        
        
                            
    end
    
    figure    
    subplot(2,1,1)
        
    subplot(2,1,2)
end

%% Visualize warps

trs = [];
nrns = [];
iter = [];
cost_mat = 0;
num_wv = 1;

if cost_mat
    for nrn_idx = 1:length(nrns)
        nrn = nrns(nrn_idx);
        for tr = trs
            for wv = 1:num_wv
            
            min_val = prctile(reshape(zscore(cost_matrix{tr,wv,iter}(:,:,end)),[],1),1,1);
            max_val = prctile(reshape(zscore(cost_matrix{tr,wv,iter}(:,:,end)),[],1),99,1);
            
%             min_val = -1000;
%             max_val = -950;
            
            plot_signals_and_matrix(predictions{tr,nrn,iter}, ...
                spikes_warped{tr,nrn,1}, ...
                ... zscore(cost_matrix{tr,wv,iter}(:,:,1)), ...
                cost_matrix{tr,wv,iter}(:,:,1), ...
                warp_path_matrix{tr,wv,iter})%, ...
%                 [min_val, max_val]);
            
            end
        end
    end
else
    for nrn_idx = 1:length(nrns)
        nrn = nrns(nrn_idx);
        for tr = trs
            for wv = 1:num_wv
            
            min_val = prctile(acc_matrix{tr,wv,iter}(:,end),1,1);
            max_val = prctile(acc_matrix{tr,wv,iter}(:,1),99,1);
            
            plot_signals_and_matrix(predictions{tr,nrn,1}, ...
                spikes_warped{tr,nrn,1}, ...
                acc_matrix{tr,wv,iter}(:,:), ...
                warp_path_matrix{tr,wv,iter});%, ...
%                 [min_val, max_val]);
            end
        end
    end
end

%% Visualize warped covariates

trs = [];
iter1 = 1;
iter2 = 2;
crange = [0 2];

for tr_idx = 1:length(trs);
    tr = trs(tr_idx);

    figure
    subplot(1,2,1)
    imagesc(X_cell{tr,1},crange)
    subplot(1,2,2)
    imagesc(X_cell{tr,2},crange)    
    
end

%% Visualize summed cost functions and spike rasters

% trs = [39 92 142]; %100-120 are good
trs = [];
iter = 4;
wv = 1;
reach_dir = [-2*pi/8];

for tr_idx = 1:length(trs)
    tr = trs(tr_idx);
    if abs(Data_win.target_dir{tr}-reach_dir) > pi/8
%         continue
    end
    figure('units','normalized','outerposition',[.4 0 .4 1])
    subplot(2,1,1)
%     temp = sum(cost_matrix{tr,iter}(:,:,1:195),3);
    temp = acc_matrix{tr,wv,iter};
    temp(logical(warp_path_matrix{tr,wv,iter})) = prctile(median(temp),10);
    imagesc(temp)
    
%     temp2 = median(median(temp))*warp_path_matrix{tr,wv,iter};
%     imagesc(zscore(temp))
%     imagesc(temp + temp2)
    axis square
    title(['Summed cost matrix, trial: ' num2str(tr)])
    
    subplot(2,1,2)
    temp = cell2mat(spikes_warped(tr,:,1))';
    imagesc(temp,[0 2])
    axis square
    title('Spikes rasters (all neurons)')
end

%% Visualize summed cost functions, spike rasters, and unwarped predictions

% trs = [39 92 142]; %100-120 are good
trs = [];
iter = 10;
wv = 1;
reach_dir = [-2*pi/8];

for tr_idx = 1:length(trs)
    tr = trs(tr_idx);
    if abs(Data_win.target_dir{tr}-reach_dir) > pi/8
%         continue
    end
    figure('units','normalized','outerposition',[.1 0 .6 1])
    subplot(2,2,1)
%     temp = sum(cost_matrix{tr,iter}(:,:,1:195),3);
    temp = acc_matrix{tr,wv,iter};
    temp(logical(warp_path_matrix{tr,wv,iter})) = prctile(median(temp),10);
    imagesc(temp)
    
%     temp2 = median(median(temp))*warp_path_matrix{tr,wv,iter};
%     imagesc(zscore(temp))
%     imagesc(temp + temp2)
    axis square
    title(['Summed cost matrix, trial: ' num2str(tr)])
    
    subplot(2,2,3)
    temp = cell2mat(spikes_warped(tr,:,1))';
    imagesc(temp,[0 2])
    axis square
    title('Spikes rasters (all neurons)')
    
    subplot(2,2,2)
    temp2 = cell2mat(predictions(tr,:,1))';
    temp4 = nan(2*size(temp2,1),size(temp2,2));
    temp4(1:2:end,:) = temp2;
    temp4(2:2:end,:) = temp;
    imagesc(temp2,[0 .5])
    axis square
    title('Spike predictions before un-warp')
    
    subplot(2,2,4)
    temp3 = cell2mat(predictions(tr,:,iter))';
    temp5 = nan(2*size(temp3,1),size(temp3,2));
    temp5(1:2:end,:) = temp3;
    temp5(2:2:end,:) = temp;
    imagesc(temp3,[0 .5])
    axis square
    title('Spike predictions after un-warp')
    
end

%% Visualize how predicted raster has changed

trs = []; 
iter = max_iter;
wv = 1;
reach_dir = [-4:4]*pi/8; % Good ones: -1pi/8, 

for tr_idx = 1:length(trs)
    tr = trs(tr_idx);
    
    if abs(Data_win.target_dir{tr}-reach_dir) > pi/8
        continue
    end
    
    figure('units','normalized','outerposition',[.1 0 .3 1],'Name',['Trial: ' num2str(tr)])
%     subplot(2,2,1)
% %     temp = sum(cost_matrix{tr,iter}(:,:,1:195),3);
%     temp = acc_matrix{tr,wv,iter};
%     temp(logical(warp_path_matrix{tr,wv,iter})) = prctile(median(temp),10);
%     imagesc(temp)
%     
% %     temp2 = median(median(temp))*warp_path_matrix{tr,wv,iter};
% %     imagesc(zscore(temp))
% %     imagesc(temp + temp2)
%     axis square
%     title(['Summed cost matrix, trial: ' num2str(tr)])
    
    colormap(bone)

    subplot(4,1,2)
    temp = cell2mat(spikes_warped(tr,:,1))';
    imagesc(temp,[0 2])
    axis square
%     title('Spikes rasters (all neurons)')
    
    subplot(4,1,1)
    
    temp2 = cell2mat(predictions(tr,:,1))';
%     temp2 = cell2mat(predictions_SH(tr,:,1))'; % SH without time warping
    
    temp4 = nan(2*size(temp2,1),size(temp2,2));
    temp4(1:2:end,:) = temp2;
    temp4(2:2:end,:) = temp;
    imagesc(temp2,[0 .3])
    axis square
%     title('Spike predictions before un-warp')
    
    subplot(4,1,3)
    
    temp3 = cell2mat(predictions(tr,:,iter))'; % Normal
%     temp3 = cell2mat(predictions_SH(tr,:,1))'; % SH without time warping
    
    temp5 = nan(2*size(temp3,1),size(temp3,2));
    temp5(1:2:end,:) = temp3;
    temp5(2:2:end,:) = temp;
    imagesc(temp3,[0 .3])
    axis square
%     title('Spike predictions after un-warp')

    colormap(parula)
    subplot(4,1,4)
    imagesc(temp3-temp2,[-.2 .2])
    axis square
    
end

%% Visualize predicted rasters

nrns = [];

if real_data
    [~,idx_sort] = sort(cell2mat(Data_win.target_dir));
else
    [~,idx_sort] = sort(Data_win.target_dir);
end

iter1 = 1;
iter2 = 10;
iter3 = 11;

for nrn_idx = 1:length(nrns)
    nrn = nrns(nrn_idx);
    
    figure
    subplot(2,2,1)
    imagesc(reshape(cell2mat(spikes_all(idx_sort,nrn)),bins_per_trial(1),[])',[0 2])
    subplot(2,2,2)
    imagesc(reshape(cell2mat(predictions(idx_sort,nrn,iter1)),bins_per_trial(1),[])',[0 2])
    subplot(2,2,3)
    imagesc(reshape(cell2mat(predictions(idx_sort,nrn,iter2)),bins_per_trial(1),[])',[0 2])
    subplot(2,2,4)
    imagesc(reshape(cell2mat(predictions(idx_sort,nrn,iter3)),bins_per_trial(1),[])',[0 2])
    
end

%% See if time warping helped Pseudo R2

do_this = 0;
num_bs = 1000;

if do_this

    new_plot = 1;
    nrns = 1:num_nrn;
    
    color = [1 0 0];
    ax = [-.01 .25 -.01 .25];
    
    iter1 = 1;
    iter2 = 3;

    % Plot option 1
    if ~new_plot
        
        temp = cell2mat(pseudo_R2(:,1));
        
        data1 = cell2mat(pseudo_R2(:,iter1));
        data2 = cell2mat(pseudo_R2(:,iter2));
        
        % data1 = cell2mat(pseudo_R2_SH(:,1)); % Spike history without time warping
        
        [fp,~,fi] = glmfit(data1,data2);
        bs_fit = bootci(num_bs,{@glmfit,data1,data2});
        slope_ci = bs_fit(:,2)';
    
        figure
        
        scatter(data1(:,1),data2(:,1))
        hold on
        dd = 0:.01:max(data2(:,1));
        plot(dd,dd)
        plot(dd,fp(2)*dd + fp(1))
        hold off
        
        xlabel('PR2 (before warp)')
        ylabel('PR2 (after warp)')
        title(['Slope: ' num2str(fp(2)) '  [' num2str(slope_ci) ']'])
        
        axis([0 1.2*max(dd) 0 1.2*max(dd)])
        
        axis square
    else
    % Plot option 2
    options = statset('UseParallel',false);
    PR2_x = nan(length(nrns),3);
    PR2_y = nan(length(nrns),3);
    
        for nrn_idx = 1:length(nrns)
            disp(['Analyzing neuron: ' num2str(nrn_idx) '/' num2str(length(nrns))])

            y = y_all(:,nrn_idx);

            % To compare warped and not warped
            % Choose predictions for comparison (the "better" one)
            pred_2 = nansum(cell2mat(predictions_combined{nrn_idx,iter2}))';

            % Choose predictions for comparison (the "worse" one)
            pred_1 = nansum(cell2mat(predictions_combined{nrn_idx,iter1}))'; % Normal time warping (last iteration)
%                 pred_1 = nansum(cell2mat(predictions_combined_SH{nrn_idx,1}))'; % For comparison with SH without time warping (for effect size of time warping)
            %     pred_1 = nansum(cell2mat(predictions_unt_combined{nrn_idx,iter1}))'; % For comparison with untuned model (for effect size of tuning)


            PR2_x(nrn_idx,1) = compute_pseudo_R2(y,pred_1,mean(y));
            PR2_y(nrn_idx,1) = compute_pseudo_R2(y,pred_2,mean(y));
            
%             PR2_x(nrn_idx,2:3) = bootci(1,{@compute_pseudo_R2,y,pred_1,mean(y)},'Options',options);
%             PR2_y(nrn_idx,2:3) = bootci(1,{@compute_pseudo_R2,y,pred_2,mean(y)},'Options',options);
            
            % fast version - no bootstrapping or invoking bootci
            PR2_x(nrn_idx,2:3) = PR2_x(nrn_idx,1)*[1,1];
            PR2_y(nrn_idx,2:3) = PR2_y(nrn_idx,1)*[1,1];
            
        end
        
        scatter_better_2D(PR2_x(:,1),PR2_x(:,2:3),PR2_y(:,1),PR2_y(:,2:3),ax,color);
        
    end
end

%% Plot test LLHD over iterations

do_this = 0;
iters = 1:iter;

if do_this
    figure
    plot(LLHD_mean(iters))
    title('Test LLHD')
    xlabel('Iteration #')
    ylabel('Mean LLHD')
end

%% Relative Pseudo R2 between conditions

plot_this = 0;

nrns = 1:num_nrn;
do_bs = 0;
num_bs = 1000;
do_untuned_model = 0;

comp1 = '1'; % options: '1', 'SH'
comp2 = 'conv'; % options: conv, conv+SH

% -------- to do: automatic iteration choosing

last_iter = find(~isnan(LLHD_mean),1,'last');

if strcmp(comp2,'conv+SH') && include_spk_history % if spk history was used (last iter)
    iter2 = last_iter;
elseif strcmp(comp2,'conv') && include_spk_history % if spk history was used and you don't want that one
    iter2 = last_iter - 1;
elseif strcmp(comp2,'conv') && ~include_spk_history
    iter2 = last_iter;
end

if strcmp(comp1,'SH')
    iter1 = 'SH';
elseif strcmp(comp1,'1')
    iter1 = 1;
end

disp(['Comparing: ' comp2 ' and ' comp1 '.'])

RPR2 = nan(length(nrns),3); 

% Set parallel options
options = statset('UseParallel',true);

for nrn_idx = 1:length(nrns)
    disp(['Analyzing neuron: ' num2str(nrn_idx) '/' num2str(length(nrns))])
    
    y = y_all(:,nrn_idx);
    
    % Choose predictions for comparison (the "better" one)
    pred_2 = nansum(cell2mat(predictions_combined{nrn_idx,iter2}))';

    % Choose predictions for comparison (the "worse" one)
    if strcmp(comp1,'1')
        pred_1 = nansum(cell2mat(predictions_combined{nrn_idx,iter1}))'; % Normal time warping (last iteration; for effect size of SH) 
    elseif strcmp(comp1,'SH')
        pred_1 = nansum(cell2mat(predictions_combined_SH{nrn_idx,1}))'; % For comparison with SH without time warping (for effect size of time warping)
    end
    
    if do_untuned_model
%         fname = 'Predictions_PMd_fit42_untuned.mat'; % Monkey 1
%         fpath = 'C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Data\2-11-16\'; % Monkey 1
        
        fname = 'Predictions_PMd_fit44_untuned.mat'; % Monkey 2
        fpath = 'C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Data\2-24-16\'; % Monkey 2
        
        vars_to_load = 'predictions_combined';
        
        temp = load([fpath fname],vars_to_load);
        predictions_unt_combined = temp.predictions_combined;
        pred_1 = nansum(cell2mat(predictions_unt_combined{nrn_idx,iter2}))'; % For comparison with untuned model (for effect size of tuning)
    end
    
    % Calculate RPR2
    RPR2(nrn_idx,1) = compute_rel_pseudo_R2(y,pred_1,pred_2);
    
    % Calculate error bars if desired
    if do_bs
        RPR2(nrn_idx,2:3) = bootci(num_bs,{@compute_rel_pseudo_R2,y,pred_1,pred_2},'Options',options);
        U(nrn_idx) = RPR2(nrn_idx,3) - RPR2(nrn_idx,1);
        L(nrn_idx) = RPR2(nrn_idx,1) - RPR2(nrn_idx,2);
    else
        RPR2(nrn_idx,2:3) = RPR2(nrn_idx,1)*[1 1];
        U(nrn_idx) = 0;
        L(nrn_idx) = 0;
    end
    
    % To compare tuned and untuned (for calibration to PR2)
%     pred_2 = nansum(cell2mat(predictions_combined{nrn_idx,1}))';
%     pred_1 = nansum(cell2mat(predictions_combined_untuned{nrn_idx,1}))';
%     RPR2_untuned(nrn_idx,1) = compute_rel_pseudo_R2(y,pred_1,pred_2);
    
%     scatter(nrn_idx,RPR2_untuned(nrn_idx,1)); hold on
    
end

% Average of warp/un-warp
sem = std(RPR2(:,1)) / sqrt(num_nrn);
mean_effect = mean(RPR2(:,1));
median_effect = median(RPR2(:,1));

% Significance testing of effect


% Average of tuned/untuned
% sem = std(RPR2_untuned(:,1)) / sqrt(num_nrn);
% scatter(num_nrn+7,mean(RPR2_untuned(:,1)))
% errorbar(num_nrn+7,mean(RPR2_untuned(:,1)),2*sem)

if plot_this
    figure
%     scatter(nrn_idx,RPR2(nrn_idx,1)); hold on
    scatter(nrns,RPR2(nrns,1)); hold on
    errorbar(nrns,RPR2(nrns,1),L,U,'o'); hold on
    scatter(num_nrn+5,mean(RPR2(:,1))); hold on
    errorbar(num_nrn+5,mean(RPR2(:,1)),2*sem,'o'); hold on
    
    line([0 length(nrns)+10],[0 0]); hold on
    
    % line([0 length(nrns)],[mean_effect mean_effect])
    axis([0 length(nrns)+10, -.01 .05])
    title(['Median RPR2: ' num2str(median_effect)])
    hold off
end

% Save results
if quest
    Results(rep).RPR2 = RPR2;
    Results(rep).mean_effect = mean_effect;
    Results(rep).mean_effect_sem = sem;
    Results(rep).median_effect = median_effect;
else
    Results.RPR2 = RPR2;
    Results.mean_effect = mean_effect;
    Results.mean_effect_sem = sem;
    Results.median_effect = median_effect;
end

%% Visually compare simulated warps with inferred warps

trials_to_use = [];
iter = 5;
rep = 1;

plot_this = 0;

for tr_idx = 1:length(trials_to_use)
    tr = trials_to_use(tr_idx);
    
    if plot_this
        figure
        subplot(4,5,[1:10])
    %     imagesc(Data_win.warp_matrices{tr})
        imagesc(Results(1).warp_matrices_sim{tr})
        axis square
    end
    
    for wv = 1:num_WV

%         mat_diff = compare_warp_matrices(Data_win.warp_matrices{tr},warp_path_matrix{tr,wv,iter});
        mat_diff = compare_warp_matrices(Results(1).warp_matrices_sim{tr},Results(rep).warp_matrices_inf{tr,wv,iter});
        
%         mat_diff_no_warp = compare_warp_matrices(flipud(eye(size(Data_win.warp_matrices{tr}))),warp_path_matrix{tr,wv,iter});
        mat_diff_no_warp = compare_warp_matrices(flipud(eye(size(Results(1).warp_matrices_sim{tr}))),Results(rep).warp_matrices_inf{tr,wv,iter});
        
        if plot_this
            subplot(4,5,10+wv)
            imagesc(Results(rep).warp_matrices_inf{tr,wv,iter})
            axis square
            
            title(['Distance: ' num2str(mat_diff) '; ' num2str(mat_diff_no_warp)])
        end
        
    end
end

%% Statistically compare simulated warps with inferred warps

last_iter = find(~isnan(LLHD_mean),1,'last');

iter = last_iter;
trials_to_use = 1:num_trials;

plot_this = 0;

% Initialize
warp_matrix_distances = nan(num_trials,num_WV);
warp_matrix_distances_shuffle = nan(num_trials,num_WV);
trials_rand = randperm(num_trials); % For shuffle trials


% for rep = 1:10
    
for tr_idx = 1:length(trials_to_use)
    tr = trials_to_use(tr_idx);
       
    for wv = 1:num_WV
%         warp_matrix_distances(tr,wv) = compare_warp_matrices(warp_path_matrix{tr,wv,iter},Data_win.warp_matrices{tr});
        warp_matrix_distances(tr,wv) = compare_warp_matrices(Results(rep).warp_matrices_inf{tr,wv,iter},Results(1).warp_matrices_sim{tr});
        
        tr_rand = trials_rand(tr);
%         warp_matrix_distances_shuffle(tr,wv) = compare_warp_matrices(warp_path_matrix{tr,wv,iter},Data_win.warp_matrices{tr_rand});
        warp_matrix_distances_shuffle(tr,wv) = compare_warp_matrices(Results(rep).warp_matrices_inf{tr,wv,iter},Results(1).warp_matrices_sim{tr_rand});
    end
end

% Using each WV as unique data points (not technically correct since each one is not indep)
% dist_reshape = reshape(warp_matrix_distances,[],1);
% dist_shuffle_reshape = reshape(warp_matrix_distances_shuffle,[],1);

% Averaging across WV instead of using each one individually
dist_reshape = mean(warp_matrix_distances,2);
dist_shuffle_reshape = mean(warp_matrix_distances_shuffle,2);

% Plot CIs
median_dist = median(dist_reshape);
median_dist_CI = bootci(1000,@median,dist_reshape);

median_dist_shuffle = median(dist_shuffle_reshape);
median_dist_shuffle_CI = bootci(1000,@median,dist_shuffle_reshape);

% Plot histograms
if plot_this
    figure
    h1 = histogram(dist_reshape); hold on
    h2 = histogram(dist_shuffle_reshape);

    figure 
    scatter(1,median_dist); hold on
    scatter([1 1],median_dist_CI); hold on
    scatter(2,median_dist_shuffle); hold on
    scatter([2 2],median_dist_shuffle_CI); hold off
    xlim([0 3])
    
%     figure
%     scatter(1,Results(rep).mean_dist); hold on
%     scatter([1 1],Results(rep).mean_dist_CI); hold on
%     scatter(2,Results(rep).mean_dist_shuffle); hold on
%     scatter([2 2],Results(rep).mean_dist_shuffle_CI); hold off
%     xlim([0 3])
end

% Save results
Results(rep).param_recovery.tau.median_dist_tau = median_dist;
Results(rep).param_recovery.tau.median_dist_tau_CI = median_dist_CI;
Results(rep).param_recovery.tau.median_dist_tau_shuffle = median_dist_shuffle;
Results(rep).param_recovery.tau.median_dist_tau_shuffle_CI = median_dist_shuffle_CI;

%% Statistically compare simulated beta and inferred beta

last_iter = find(~isnan(LLHD_mean),1,'last');

iter = last_iter;

plot_this = 0;

% Initialize
beta_distances = nan(num_nrn,num_CV);
beta_distances_noDTW = nan(num_nrn,num_CV);
beta_corr = nan(num_nrn,num_CV);
beta_corr_noDTW = nan(num_nrn,num_CV);
beta_corr_shuffle = nan(num_nrn,num_CV);


% Loop
for nrn_idx = 1:length(neurons)
    nrn_num = neurons(nrn_idx);

    for cv = 1:num_CV
        vec1 = Results(1).beta_sim(nrn_num,:);
        vec1 = vec1(2:end); % Disregard intercept; influenced by offset parameter, and isn't important part of RF anyway
        
        vec2 = Results(rep).beta_inf(nrn_num,:,cv);
        vec2 = vec2(2:end); % Disregard intercept
        
        vec3 = fit_parameters{nrn_num,1}{cv}';
        vec3 = vec3(2:end); % Disregard intercept
        
        idx_rand = randperm(length(vec1)); % Generate shuffle indexes
        vec4 = vec2(idx_rand); % Shuffle order of fit parameters
        
        beta_distances(nrn_num,cv) = calc_vector_distance(vec2,vec1); % angle between vectors, using definition of dot product; compare recovered to simulated
        beta_distances_noDTW(nrn_num,cv) = calc_vector_distance(vec3,vec1); % angle between vectors, using definition of dot product; compare GLM-only to simulated
        
        beta_corr(nrn_num,cv) = corr(vec1',vec2'); % linear correlation: simulated and recovered
        beta_corr_noDTW(nrn_num,cv) = corr(vec1',vec3'); % linear correlation: simulated and recovered without DTW
        beta_corr_shuffle(nrn_num,cv) = corr(vec1',vec4'); % linear correlation: simulated and shuffled recovered
        
    end
end

beta_distances_reshape = mean(beta_distances,2); % average over cv folds, so as not to treat them as indep data points
beta_distances_noDTW_reshape = mean(beta_distances_noDTW,2); % average over cv folds, so as not to treat them as indep data points
beta_corr_reshape = mean(beta_corr,2);
beta_corr_noDTW_reshape = mean(beta_corr_noDTW,2);
beta_corr_shuffle_reshape = mean(beta_corr_shuffle,2);

% Plot CIs
median_dist_beta = median(beta_distances_reshape);
median_dist_beta_CI = bootci(1000,@median,beta_distances_reshape);

median_dist_beta_noDTW = median(beta_distances_noDTW_reshape);
median_dist_beta_noDTW_CI = bootci(1000,@median,beta_distances_noDTW_reshape);

median_corr_beta = median(beta_corr_reshape);
median_corr_beta_CI = bootci(1000,@median,beta_corr_reshape);
median_corr_beta_noDTW = median(beta_corr_noDTW_reshape);
median_corr_beta_noDTW_CI = bootci(1000,@median,beta_corr_noDTW_reshape);
median_corr_beta_shuffle = median(beta_corr_shuffle_reshape);
median_corr_beta_shuffle_CI = bootci(1000,@median,beta_corr_shuffle_reshape);

% Plot histograms
if plot_this
    figure
    h1 = histogram(beta_distances_reshape)

    figure 
    scatter(1,median_dist_beta); hold on
    scatter([1 1],median_dist_beta_CI); hold on
    xlim([0 3])
    
%     figure
%     scatter(1,Results(rep).mean_dist); hold on
%     scatter([1 1],Results(rep).mean_dist_CI); hold on
%     scatter(2,Results(rep).mean_dist_shuffle); hold on
%     scatter([2 2],Results(rep).mean_dist_shuffle_CI); hold off
%     xlim([0 3])
end

% Save results
Results(rep).param_recovery.beta.median_dist_beta = median_dist_beta;
Results(rep).param_recovery.beta.median_dist_beta_CI = median_dist_beta_CI;
Results(rep).param_recovery.beta.median_dist_beta_noDTW = median_dist_beta_noDTW;
Results(rep).param_recovery.beta.median_dist_beta_noDTW_CI = median_dist_beta_noDTW_CI;
Results(rep).param_recovery.beta.median_corr_beta = median_corr_beta;
Results(rep).param_recovery.beta.median_corr_beta_CI = median_corr_beta_CI;
Results(rep).param_recovery.beta.median_corr_beta_noDTW = median_corr_beta_noDTW;
Results(rep).param_recovery.beta.median_corr_beta_noDTW_CI = median_corr_beta_noDTW_CI;
Results(rep).param_recovery.beta.median_corr_beta_shuffle = median_corr_beta_shuffle;
Results(rep).param_recovery.beta.median_corr_beta_shuffle_CI = median_corr_beta_shuffle_CI;

%% Look at predictions over iterations

trials_to_use = [];
nrn = [];
iterations = [];

for tr_idx = 1:length(trials_to_use)
    tr = trials_to_use(tr_idx);
    
    figure
    subplot(1,(length(iterations)+1),1)
    imagesc(spikes_warped{tr,nrn,1},[0 2]); hold on
    for iter_idx = 1:length(iterations);
        iter = iterations(iter_idx);
        subplot(1,(length(iterations)+1),iter+1)
        imagesc(predictions{tr,nrn,iter},[0 .7]); hold off
    end
    title(['Trial: ' num2str(tr)])
    
end

%% Save additional results

if ~quest
    save([Results_fpath Results_fname],'Results','-v7.3')
else
    save(Results_fname,'Results','-v7.3') % Save only Results

    % Break things up for Quest. Bastard.
    save([Results_fname '_LLHD'], ...
        'LLHD_test','LLHD_test_diag', ...
        '-v7.3') 
    
    save([Results_fname '_pred'], ...
        'predictions_combined','predictions', ...
        'predictions_combined_SH','predictions_SH', ...
        'predictions_combined2','predictions2', ...
        '-v7.3') 
    
    save([Results_fname '_pr2'], ...
        'pseudo_R2','pseudo_R2_SH', ...
        '-v7.3') 
    
    save([Results_fname '_warp'], ...
        'warp_path','warp_path_matrix', ...
        '-v7.3') 
    
end

disp(['Extended Results file saved: ' Results_fname])

%% End for loop for simulated data

end % For "repetition"