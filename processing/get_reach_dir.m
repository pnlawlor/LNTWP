function [ dir_reach ] = get_reach_dir( kinematics )
    
    x_pos = kinematics(:,1);
    y_pos = kinematics(:,2);    
    x_vel = kinematics(:,3);
    y_vel = kinematics(:,4);
    x_vel_sm = smooth(x_vel);
    y_vel_sm = smooth(y_vel);
    
    distance = sqrt(x_pos.^2 + y_pos.^2);
    distance_sm = smooth(distance);
    [~, bin_farthest_point] = max(distance_sm);
    
    speed = sqrt(x_vel.^2 + y_vel.^2);
    speed_sm = smooth(speed(1:bin_farthest_point));
    [~, bin_max_spd] = max(speed_sm);
    
    dir_reach = cart2pol(x_vel_sm(bin_max_spd),y_vel_sm(bin_max_spd));

end

