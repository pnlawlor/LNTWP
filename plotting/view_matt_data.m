%% View everything

trials = 11:20;

reach_tr = cell2mat(Data2.trial_num);

for tr_idx = 1:length(trials)
    tr = trials(tr_idx);

    idx_reach = find(reach_tr==tr);
    
    figure

    subplot(5,5,2:5)

    x_vel = Data.kinematics{tr}(:,3);
    y_vel = Data.kinematics{tr}(:,4);
    
%     plot(x_vel); hold on
%     plot(y_vel); hold on

    thresh = 5;

    spd = sqrt(x_vel.^2 + y_vel.^2); 
%     plot(find(spd<thresh),spd(spd<thresh),'LineWidth',3); hold on
%     plot(find(spd>=thresh),spd(spd>=thresh),'LineWidth',.5); hold on
    plot(spd); hold on
    title(['Trial: ' num2str(tr)])
    
    for rch = 1:length(idx_reach)
        plot((10+2*rch)*Data2.time_window{idx_reach(rch)}); hold on
        idx_window_on = find(Data2.time_window{idx_reach(rch)},1);
        idx_target_on = find(Data2.target_on{idx_reach(rch)});
        idx_line = idx_window_on + idx_target_on;
        line([idx_line idx_line],[0 10],'LineWidth',2); hold on
    end
    
    
    line([1, length(Data.kinematics{tr}(:,4))], [8, 8],'LineStyle','--')
    xlim([1 length(Data.kinematics{tr}(:,4))])
    ylim([0 30])
    
    subplot(5,5,[17 18 19 20 22 23 24 25])
    imagesc(Data.neural_data_M1{tr},[0 2])
    
    subplot(5,5,[7 8 9 10 12 13 14 15])
    imagesc(Data.neural_data_PMd{tr},[0 2])
    
    hold off

end

%% view kinematics

trials = 4;

reach_tr = cell2mat(Data2.trial_num);

for tr_idx = 1:length(trials)
    tr = trials(tr_idx);
    
    idx_reach = find(reach_tr==tr);
    
    x_pos = Data.kinematics{tr}(:,1);
    y_pos = Data.kinematics{tr}(:,2);
    
    figure
    plot(x_pos,y_pos); hold on
    
    for i = 1:length(idx_reach)
        rch = idx_reach(i);
        scatter(Data2.reach_pos_st{rch}(1),Data2.reach_pos_st{rch}(2),'b^')
        scatter(Data2.reach_pos_end{rch}(1),Data2.reach_pos_end{rch}(2),'rv')
        disp(['Reach directon: ' num2str(180/pi*Data2.reach_dir{rch})])
    end
    
    axis([-15 15 -15 15])    
    
end

%% View reach by reach

reaches = 81:100;

for reach_idx = 1:length(reaches)
    reach = reaches(reach_idx);
    
    figure
    
    subplot(5,1,1)
    x_vel = Data2.kinematics{reach}(:,3);
    y_vel = Data2.kinematics{reach}(:,4);
%     plot(x_vel); hold on
%     plot(y_vel); hold on
    plot(smooth(sqrt(x_vel.^2+y_vel.^2))); hold on
    
    subplot(5,1,[2 3])
    imagesc(Data2.neural_data_M1{reach}); hold on
    
    subplot(5,1,[4 5])
    imagesc(Data2.neural_data_PMd{reach})
    hold off   
    
end