function [ trials_good ] = reach_by_st_loc( Data, box ) % box = [x1 x2 y1 y2], x2 > x1, y2 > y1
    
    num_trials = size(Data.kinematics,1);
    trials_good = nan(num_trials,1);
    
    x1 = box(1); x2 = box(2); y1 = box(3); y2 = box(4);
    
    for tr = 1:num_trials
        x_st = Data.reach_pos_st{tr}(1);
        y_st = Data.reach_pos_st{tr}(2);
        
        if (x_st > x1) && (x_st < x2) && ...
                (y_st > y1) && (y_st < y2)
            trials_good(tr) = 1;
        else
            trials_good(tr) = 0;
        end
    end
    
    trials_good = logical(trials_good);

end
