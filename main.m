%% User inputs 1

quest = 0;

% Choose model for GLM
model_num = 8; % 8 = basic, 10 = comprehensive without directional tuning, 11 = comprehensive

% Choose properties of data to use
FR_thresh_Hz = 2; 
dt = .01;
FR_thresh = FR_thresh_Hz*dt;

real_data = 1; % 1 = real data from monkey; 0 = simulated data
num_nrn_desired = 200;
warp_type = 'strong'; % If using simulated data, selects data file. Use 'strong', 'weak', or 'none'.

%% Load scripts

% Initialize seed
rng('default')

% Relative paths
addpath(genpath('./'))
    
%% Choose brain region/model
% 8 = Basic model: directional tuning, no reach velocity, no spike history
% 10 = Comprehensive model without tuning
% 11 = Comprehensive model: directional tuning, reach velocity and speed

%% Load data

load_data

%% User inputs 2

rng('default')

% User-chosen priors for transition probabilities in alignment matrix
% Bottom-left-most entry in the this matrix is the "starting location", 
% from which other transitions are defined. I.e., in the default transition 
% matrix [0,1; 0,1; 0,1]: entry in (3,1) (numbered from the top-left 
% corner) is the "starting location" and should always be set to 0. The 
% entry in (3,2) represents a horizontal step in the alignment matrix 
% (right 1, up 0). The entry in (2,2) represents a diagonal step 
% (right 1, up 1). The entry in (1,2) represents a steeper-than-diagonal 
% step (right 1, up 2). 

transition_priors = [0 1; 0 1; 0 1]; 
transition_prior_scale = 1; % Determines the relative weight of the prior relative to the neural data. Default 1. 
transition_priors = transition_priors/sum(sum(transition_priors)); % Normalizes the transition matrix so that they sum to 1. 
transition_priors = log(transition_priors); % Takes the log of the transition matrix; this is how it appears in the optimization equations.
transition_priors = transition_prior_scale*transition_priors; 

fit_method = 'lassoglm'; % GLM fitting method. Options: glmfit, lassoglm, glmnet. Default is lassoglm (regularized). Glmfit is unregularized. Glmnet is the fastest but requires downloading the package and can sometimes cause segfaults.
num_CV = 2; % Number of trial-wise cross-validation folds. Default is 2.
num_WV = 10; % Number of neuron-wise cross-validation folds. Default is 10.

tol_type = 'relTol'; % Rule that the code uses to terminate alternation for fitting LNP+DTW model. Default is relative tolerance, 'relTol'. Options: 'relTol','absTol','maxIter'
conv_thresh = 1e-4; % Tolerance for convergence detection when using 'relTol'. Equals dLLHD/LLHD. Default is 1e-4.

max_iter = 10; % Max number of altnerating iterations. Default is 10.

alpha = 0.01; % Regularization type (L1 vs L2). Alpha small = L2; alpha 1 = L1. Default is 0.01.
lambda = 0.05; % Regularization strength. Default is 0.05.

trials = 1:max_num_trials;
neurons = 1:num_nrn;

do_shuffle = 0; % Shuffle control: apply learned warps to wrong trials

include_spk_history = 0; % Include spike history terms. This adds an extra iteration. 
num_spk_history_bf = 5; % Number of temporal basis functions to use for spike history

% Simulation parameters
num_rep = 10; % Number of times to resample neurons and refit model. 
% neurons = 1:num_nrn_desired; % For simulations

% Some file save path initialization
if real_data
    Results_fname_datatype = 'RealData';
    Results_fname_base = [Results_fname_datatype '_' num2str(num_nrn_desired) 'nrn_' num2str(length(trials)) 'tr'];
else
    Results_fname_datatype = 'SimData';
    Results_fname_base = [Results_fname_datatype '_' num2str(num_nrn_desired) 'nrn_' num2str(length(trials)) 'tr' '_' warp_type];
end

Results_fpath = '.\results\';
datetimestr = datestr(now,'mm-dd-yy--HH-MM');

%% Start repetitions (for simulations only)

% for rep = 1:num_rep

%% Initialize

initialize 

%% Process

process

%% Fit model

fit_model

%% ==== Visualization and Analysis Tools ====

%% Visualize: spike raster, pre-warp predictions, warped predictions, difference between the two (Fig 7)

% This tools helps to visualize the impact of time warping on individual
% trials. In the manuscript, it generated Figure 7. For individual trials,
% it displays the LNP-only predictions (i.e., pre-warp; row 1), the spike 
% raster (row 2), and the LNP+DTW predictions (i.e., post-warp; row 3). It 
% also plots the difference in between the pre-warp and post-warp 
% predictions (row 4).

% Example trials
% Trials 28, 154, 276, and 315 are used for Fig 7 in the manuscript

trs = [28, 154, 315, 276];

reach_dir = [-4:4]*pi/8; % Use this to limit reaches to those towards a particular direction. Default is to use all reaches.

% Automatic iteration choosing
last_iter = find(~isnan(LLHD_mean),1,'last'); % Get last iteration number
iter = last_iter - do_shuffle - include_spk_history; % Don't include spike history or shuffle iterations

for tr_idx = 1:length(trs)
    tr = trs(tr_idx);
    
    if abs(Data_win.target_dir{tr}-reach_dir) > pi/8 % Discard reaches not in desired direction
        continue
    end
    
    figure('units','normalized','outerposition',[.1 0 .3 1],'Name',['Trial: ' num2str(tr)])

    % Spike raster
    ax2 = subplot(4,1,2);
    colormap(ax2,flipud(bone))
    temp = cell2mat(spikes_warped(tr,:,1))';
    imagesc(temp,[0 2])
    axis square
    title('Spikes')
    
    ax1 = subplot(4,1,1);
    colormap(ax1,flipud(bone))
    
    % Pre-warp predictions
    temp2 = cell2mat(predictions(tr,:,1))';
    temp4 = nan(2*size(temp2,1),size(temp2,2));
    temp4(1:2:end,:) = temp2;
    temp4(2:2:end,:) = temp;
    imagesc(temp2,[0 .3])
    axis square
    title('LNP predictions')

    % Post-warp predictions
    ax3 = subplot(4,1,3);
    colormap(ax3,flipud(bone))
    temp3 = cell2mat(predictions(tr,:,iter))';
    temp5 = nan(2*size(temp3,1),size(temp3,2));
    temp5(1:2:end,:) = temp3;
    temp5(2:2:end,:) = temp;
    imagesc(temp3,[0 .3])
    axis square
    title('LNP+DTW predictions')

    % Difference in predictions
    ax4 = subplot(4,1,4);
    colormap(ax4,parula)
    imagesc(temp3-temp2,[-.15 .15])
    axis square
    title('Difference in predictions')
    xlabel('Time (10 ms bins) relative to start of time window')
    ylabel('Neuron number')
    
end

%% Visualize: Preds vs Spikes with warp matrix (Fig 8)

trs = [28, 154, 315, 276]; % Four trials from manuscript (28,154,276,315)
wv = 1; % Warp validation to visualize. Default is 1.

% Automatic iteration choosing
last_iter = find(~isnan(LLHD_mean),1,'last'); % Get last iteration number
iter = last_iter - do_shuffle - include_spk_history; % Don't include spike history or shuffle iterations

for tr_idx = 1:length(trs)
    tr = trs(tr_idx);
    
    preds = cell2mat(predictions(tr,:,1))'; 
    spks = cell2mat(spikes_warped(tr,:,1))';
    warp_mat = warp_path_matrix{tr,wv,iter};
    cost_mat = sum(cost_matrix{tr,wv,2},3);
    cost_mat = zscore(cost_mat); % I z-score the alignment matrix so that it's easier to visualize. The fitting uses the non-z-scored version, however.
    
    figure('Position',[520 241 581 557],'Name',['Trial: ' num2str(tr)])
    
    % Predictions (left panel)
%     ax1 = subplot(4,4,[1 5 9]);
    ax1 = subplot(5,5,[1 2 6 7 11 12]);
    colormap(ax1,flipud(bone))
    imagesc(preds,[0 .3])
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    set(gca,'view',[-90 90])
    xlabel('Time (\tau)')
    ylabel('Nrn #')

    % Spikes (bottom panel)
%     ax2 = subplot(4,4,[14:16]);
    ax2 = subplot(5,5,[18:20, 23:25]);
    colormap(ax2,flipud(bone))
    imagesc(spks,[0 2])
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    xlabel('Time (t)')
    ylabel('Nrn #')
    
    % Warp matrix (middle panel)
%     ax3 = subplot(4,4,[2:4,6:8,10:12]);
    ax3 = subplot(5,5,[3:5,8:10,13:15]);
%     colormap([1 1 1; 255/255 105/255 180/255])    
%     imagesc(warp_mat) % Plot the alignment matrix using imagesc
    colormap(ax3,hot)
    imagesc(cost_mat) % Plot the cost matrix using imagesc. I
%     superimposed the two post-hoc; I couldn't figure out how to do
%     it programmatically. Just uncomment whichever one you want to see.
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    axis square
    
end

%% Visualize: Predicted rasters from multiple models

% This tool allows you to visualize model predictions for a single neuron
% across multiple types of models. I typically use this to compare the
% original spike raster with the LNP, LNP+DTW, and LNP+DTW+SH. 
% The trials here are sorted by reach direction, so that directional tuning
% is visually apparent.

nrns = [];

if real_data
    [~,idx_sort] = sort(cell2mat(Data_win.target_dir));
else
    [~,idx_sort] = sort(Data_win.target_dir);
end

iter1 = 1; % LNP only
iter2 = 4; % LNP+DTW
iter3 = 5; % LNP+DTW + SH

% *** add automatic iteration chooser

for nrn_idx = 1:length(nrns)
    nrn = nrns(nrn_idx);
    
    spks = reshape(cell2mat(spikes(idx_sort,nrn)),bins_per_trial(1),[])';
    pred1 = reshape(cell2mat(predictions(idx_sort,nrn,iter1)),bins_per_trial(1),[])';
    pred2 = reshape(cell2mat(predictions(idx_sort,nrn,iter2)),bins_per_trial(1),[])';
    pred3 = reshape(cell2mat(predictions(idx_sort,nrn,iter3)),bins_per_trial(1),[])';
    max1 = max(max(pred1));
    max2 = max(max(pred2));
    max3 = max(max(pred3));
    predmax = max([max1,max2,max3]);
    crange_preds = [0 predmax];
    
    figure
    subplot(2,2,1)
    imagesc(spks,[0 2])
    subplot(2,2,2)
    imagesc(pred1,crange_preds)
    subplot(2,2,3)
    imagesc(pred2,crange_preds)
    subplot(2,2,4)
    imagesc(pred3,crange_preds)
    
end

%% Plot test LLHD over iterations (Fig S6)

do_this = 0; % If =1, plots. If =0, skips.
iters = 1:max_iter;

if do_this
    figure
    plot(LLHD_mean(iters))
    title('Test LLHD')
    xlabel('Iteration #')
    ylabel('Mean LLHD')
%     ylim([-2e3 0])
%     xlim([0 15])
    
    clear reltol
    for it = 2:max(iters)
        LL0 = LLHD_mean(it-1);
        LL1 = LLHD_mean(it);
        reltol(it) = abs(LL1 - LL0) / min(abs(LL0), abs(LL1));
    end
    
    figure
    plot(reltol)
    title('Relative difference in LLHD')
    xlabel('Iteration #')
    ylabel('Relative difference in LLHD')
end

%% Pseudo R2 scatterplot (Fig 9)

% This tool plots a scatterplot of the PseudoR2 of one condition vs that of
% another. This is used in the manuscript to visualize that time warping
% improves the Pseudo R2 by comparing the PseudoR2 of the LNP+DTW model
% with that of the LNP model. You can use this to compare any two
% conditions, however. 

% Guide to alternation iteration numbers: 
% 1: LNP only, no DTW
% 2 through iter_conv: Successive alternation iterations until convergence detected
% iter_conv + 1: Converged LNP+DTW model with spike history included (if
% spike history option was selected; otherwise iter_conv is the last
% iteration).

do_this = 1; % Run this tool or not
do_bs = 0; % Bootstrap the PseudoR2 confidence interval or not
num_bs = 1000; % Number of bootstrap iterations if doing BS

if do_this

    new_plot = 1;
    nrns = 1:num_nrn;
    
    color = [1 0 0];
    ax = [-.01 .25 -.01 .25];
    
    comp1 = '1'; % Options: '1', 'SH'. Use '1' for comparing LNP+DTW with LNP
    comp2 = 'conv'; % Options: conv, conv+SH

    % Automatic iteration choosing
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


    % Plot option 1 (old version)
    if ~new_plot
        
        temp = cell2mat(pseudo_R2(:,1));
        
        data1 = cell2mat(pseudo_R2(:,iter1));
        data2 = cell2mat(pseudo_R2(:,iter2));
        
        % data1 = cell2mat(pseudo_R2_SH(:,1)); % Spike history without time warping
        
        [fp,~,fi] = glmfit(data1,data2); % Fit line to scatterplot
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
            pred_2 = nansum(cell2mat(predictions_combined{nrn_idx,iter2}),1)';

            % Choose predictions for comparison (the "worse" one)
            pred_1 = nansum(cell2mat(predictions_combined{nrn_idx,iter1}),1)'; % LNP
%             pred_1 = nansum(cell2mat(predictions_combined_SH{nrn_idx,1}))'; % LNP + spike history

            % Computer Pseudo-R2
            PR2_x(nrn_idx,1) = compute_pseudo_R2(y,pred_1,mean(y));
            PR2_y(nrn_idx,1) = compute_pseudo_R2(y,pred_2,mean(y));
            
            % If bootstrapping confidence intervals
            if do_bs
                PR2_x(nrn_idx,2:3) = bootci(num_bs,{@compute_pseudo_R2,y,pred_1,mean(y)},'Options',options);
                PR2_y(nrn_idx,2:3) = bootci(num_bs,{@compute_pseudo_R2,y,pred_2,mean(y)},'Options',options);
            else
                % Fast version - no bootstrapping or invoking bootci
                PR2_x(nrn_idx,2:3) = PR2_x(nrn_idx,1)*[1,1];
                PR2_y(nrn_idx,2:3) = PR2_y(nrn_idx,1)*[1,1];
            end
            
        end
        
        % Make scatterplot (custom script)
        scatter_better_2D(PR2_x(:,1),PR2_x(:,2:3),PR2_y(:,1),PR2_y(:,2:3),ax,color); hold on;
        xlabel('Pseudo-R2 (LNP)')
        ylabel('Pseudo-R2 (LNP+DTW)')
        
    end
end

%% Relative Pseudo R2 between conditions (Fig 10)

% This tool visualizes the marginal improvement in predictive power of one
% model over another (the Relative Pseudo R2, RPR2). A positive RPR2
% indicates that model 2 improves upon model 1. A negative RPR2 indicates
% that model 2 is worse than model 1. A value of zero indicates that model
% 2 and model 1 have the same predictive power. 
% More specifically, this tool plots the RPR2 for model 2 vs model 1 for
% each individual neuron as well as the population median. A number of
% different possible comparisons are available as detailed below. 

% Conditions used in comparisons:
% Time warping, basic model: (LNP_basic + DTW) - (LNP_basic)
% Shuffle: (LNP_basic + DTW_shuffled) - (LNP_basic)
% Time warping comprehensive model: (LNP_comp + DTW + SH) - (LNP_comp + SH)
% Spike history: (LNP_comp + DTW + SH) - (LNP_comp + DTW)
% Directional tuning: (LNP_comp + DTW + SH) - (LNP_comp_untuned + DTW + SH)

% Guide to alternation iteration numbers: 
% 1: LNP only, no DTW
% 2 through iter_conv: Successive alternation iterations until convergence detected
% iter_conv + 1: Converged LNP+DTW model with spike history included (if
% spike history option was selected; otherwise iter_conv is the last
% iteration).

plot_this = 1; % Plot this or not

nrns = 1:num_nrn; % Choose which neurons to analyze. Default is 1:num_nrn
do_bs = 0; % Bootstrap confidence intervals or not. 
num_bs = 1000; % Number of bootstrap iterations if doing BS.
do_untuned_model = 0; % If comparing comprehensive model to the untuned model. This requires fitting an untuned model (model_num = 10) and using the results file below.

comp1 = '1'; % Options: '1', 'SH', 'sec2last'
comp2 = 'conv+SH'; % Options: conv, conv+SH

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
elseif strcmp(comp1,'sec2last')
    iter1 = last_iter - 1;
end

disp(['Comparing: ' comp2 ' and ' comp1 '.'])

RPR2 = nan(length(nrns),3); 

% Set parallel options
options = statset('UseParallel',false);

for nrn_idx = 1:length(nrns)
    disp(['Analyzing neuron: ' num2str(nrn_idx) '/' num2str(length(nrns))])
    
    y = y_all(:,nrn_idx);
    
    % Choose predictions for comparison (the "better" one)
    pred_2 = nansum(cell2mat(predictions_combined{nrn_idx,iter2}),1)';

    % Choose predictions for comparison (the "worse" one)
    if strcmp(comp1,'1') || strcmp(comp1,'sec2last')
        pred_1 = nansum(cell2mat(predictions_combined{nrn_idx,iter1}),1)'; % LNP (in the case of '1') or LNP+DTW (in the case of 'sec2last'). But no spike history. 
    elseif strcmp(comp1,'SH')
        pred_1 = nansum(cell2mat(predictions_combined_SH{nrn_idx,1}),1)'; % LNP+SH (used to estimate effect size of time warping)
    end
    
    % If comparison condition 1 is the untuned model (used to estimate the
    % effect size of directional tuning)
    if do_untuned_model
        % Load data for either monkey
        fname = 'PMd_fit60.mat'; % Monkey 1
        fpath = 'C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Data\6-19-17\'; % Monkey 1
        
%         fname = 'PMd_fit64.mat'; % Monkey 2
%         fpath = 'C:\Users\pnlawlor\GoogleDrive\Research\Projects\Time_warping\Data\6-19-17\'; % Monkey 2
        
        vars_to_load = 'predictions_combined';
        
        temp = load([fpath fname],vars_to_load);
        predictions_unt_combined = temp.predictions_combined;
        iter_unt = find(~cellfun(@isempty,(predictions_unt_combined(nrn_idx,:))),1,'last');
        pred_1 = nansum(cell2mat(predictions_unt_combined{nrn_idx,iter_unt}),1)'; 
    end
    
    % Calculate RPR2
    RPR2(nrn_idx,1) = compute_rel_pseudo_R2(y,pred_1,pred_2);
    
    % Calculate confidence intervals if desired
    if do_bs
        RPR2(nrn_idx,2:3) = bootci(num_bs,{@compute_rel_pseudo_R2,y,pred_1,pred_2},'Options',options);
        U(nrn_idx) = RPR2(nrn_idx,3) - RPR2(nrn_idx,1);
        L(nrn_idx) = RPR2(nrn_idx,1) - RPR2(nrn_idx,2);
    else % Or if you don't want them. Much faster.
        RPR2(nrn_idx,2:3) = RPR2(nrn_idx,1)*[1 1];
        U(nrn_idx) = 0;
        L(nrn_idx) = 0;
    end  
end

% Average of warp/un-warp
sem = std(RPR2(:,1)) / sqrt(num_nrn);
mean_effect = mean(RPR2(:,1));
median_effect = median(RPR2(:,1));

if plot_this
    figure
%     scatter(nrn_idx,RPR2(nrn_idx,1)); hold on
    scatter(nrns,RPR2(nrns,1)); hold on
    errorbar(nrns,RPR2(nrns,1),L,U,'o'); hold on
    scatter(num_nrn+5,median_effect); hold on
    
    line([0 length(nrns)+10],[0 0]); hold on
    
    % line([0 length(nrns)],[mean_effect mean_effect])
    axis([0 length(nrns)+10, -.01 .05])
    title(['Median RPR2: ' num2str(median_effect)])
    xlabel('Neuron number')
    ylabel('RPR2')
    hold off
end

% Save results
if quest % If using simulated data (done on supercomputer cluster)
    Results(rep).RPR2 = RPR2;
    Results(rep).mean_effect = mean_effect;
    Results(rep).mean_effect_sem = sem;
    Results(rep).median_effect = median_effect;
else % If using real data
    Results.RPR2 = RPR2;
    Results.mean_effect = mean_effect;
    Results.mean_effect_sem = sem;
    Results.median_effect = median_effect;
end

%% Stretch vs shift (Fig S3)

% This tool helps to visualize the "gain-like" effect of time warping. For
% every neuron and reach, the firing rate gain due to time warping 
% (FR_LNP_DTW/FR_LNP) is plotted vs the magnitude of the time warp. 

do_this = 0;

nrns = 1:num_nrn;
trials_to_use = trials;

% Initialize
TW_strength = nan(length(trials_to_use),length(nrns));
FR_rel = nan(length(trials_to_use),length(nrns));
gain = nan(length(trials_to_use),length(nrns));

comp1 = '1'; % options: '1', 'SH'
comp2 = 'conv'; % options: conv, conv+SH

last_iter = find(~isnan(LLHD_mean),1,'last');

% Automatic iteration choosing
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


if do_this
    for nrn_idx = 1:length(nrns)
        nrn_num = nrns(nrn_idx);
        
        for tr_idx = 1:length(trials)
            tr_num = trials(tr_idx);
            
            % Get strength of TW
            mat1 = warp_path_matrix{tr_num,1,iter2};
            mat2 = flipud(eye(size(mat1,1)));
            TW_strength(tr_idx,nrn_idx) = compare_warp_matrices(mat1,mat2);
            
            % Get relative change in FR
            FR1 = mean(predictions{tr_num,nrn_num,iter1});
            FR2 = mean(predictions{tr_num,nrn_num,iter2});
            FR_rel(tr_idx,nrn_idx) = (FR2 - FR1)/FR1;
            FR_rel_abs(tr_idx,nrn_idx) = abs(FR2 - FR1)/FR1;
            gain(tr_idx,nrn_idx) = FR2/FR1;
            
        end
    end

    TW_strength_reshape = reshape(TW_strength,[],1);
    FR_rel_reshape = reshape(FR_rel,[],1);
    FR_rel_abs_reshape = reshape(FR_rel_abs,[],1);
    
    FR_rel_mean = mean(FR_rel_reshape);
    FR_rel_SD = std(FR_rel_reshape);
    TW_FR_corr = corr(TW_strength_reshape,FR_rel_reshape);
    TW_FR_abs_corr = corr(TW_strength_reshape,FR_rel_abs_reshape);
    
    figure
    disp(['Delta_FR mean: ' num2str(FR_rel_mean) ', Delta_FR SD: ' num2str(FR_rel_SD)])
    disp(['Correlation between TW mag and Delta_FR: ' num2str(TW_FR_corr)])
    disp(['Correlation between TW mag and abs(Delta_FR): ' num2str(TW_FR_abs_corr)])
    subplot(1,5,1:4)
    scatter(TW_strength_reshape,FR_rel_reshape)
    
    subplot(1,5,5)
    hist(FR_rel_reshape,100)
    set(gca,'view',[90 -90])
    
    gain_reshape = reshape(gain,[],1);
    gain_mean = mean(gain_reshape);
    gain_SD = std(gain_reshape);
    TW_gain_corr = corr(TW_strength_reshape,gain_reshape);
    
    figure
    disp(['Gain mean: ' num2str(gain_mean) ', Gain SD: ' num2str(gain_SD)])
    disp(['Correlation between TW mag and Gain: ' num2str(TW_gain_corr)])
    subplot(1,5,1:4)
    scatter(TW_strength_reshape,gain_reshape)
    
    subplot(1,5,5)
    hist(gain_reshape,100)
    set(gca,'view',[90 -90])
    
end

%% Visually compare simulated warps with inferred warps

% This tool allows you to compare the simulation parameters with the
% inferred parameters for the DTW portion of the model. 

trials_to_use = [];
iter = 5; % *** automatic iteration chooser
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

        mat_diff = compare_warp_matrices(Results(1).warp_matrices_sim{tr},Results(rep).warp_matrices_inf{tr,wv,iter});
        
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

% This tool computes some comparison statistics using the simulated and
% inferred time warps. 

if ~real_data

last_iter = find(~isnan(LLHD_mean),1,'last');

iter = last_iter;
trials_to_use = 1:num_trials;

plot_this = 0;

% Initialize
warp_matrix_distances = nan(num_trials,num_WV);
warp_matrix_distances_shuffle = nan(num_trials,num_WV);
trials_rand = randperm(num_trials); % For shuffle trials
    
for tr_idx = 1:length(trials_to_use)
    tr = trials_to_use(tr_idx);
       
    for wv = 1:num_WV
        warp_matrix_distances(tr,wv) = compare_warp_matrices(Results(rep).warp_matrices_inf{tr,wv,iter},Results(1).warp_matrices_sim{tr});
        tr_rand = trials_rand(tr);
        warp_matrix_distances_shuffle(tr,wv) = compare_warp_matrices(Results(rep).warp_matrices_inf{tr,wv,iter},Results(1).warp_matrices_sim{tr_rand});
    end
end

% Averaging across WV errors instead of using each one individually (they're not independent)
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

end

%% Statistically compare simulated beta with inferred beta

% This tool computes some comparison statistics using the simulated and
% inferred time LNP parameters (aka beta, or receptive field)

if ~real_data
    
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
        
        beta_corr(nrn_num,cv) = corr(vec1',vec2'); % Linear correlation: simulated and recovered
        beta_corr_noDTW(nrn_num,cv) = corr(vec1',vec3'); % Linear correlation: simulated and recovered without DTW
        beta_corr_shuffle(nrn_num,cv) = corr(vec1',vec4'); % Linear correlation: simulated and shuffled recovered
        
    end
end

beta_distances_reshape = mean(beta_distances,2); % Average over cv folds, so as not to treat them as indep data points
beta_distances_noDTW_reshape = mean(beta_distances_noDTW,2); % Average over cv folds, so as not to treat them as indep data points
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

end

%% Save additional results

if ~quest
    save([Results_fpath Results_fname],'Results','-v7.3')
else
    save(Results_fname,'Results','-v7.3') % Save only Results

    % Break things up for Quest cluster
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

% end % For "repetition"