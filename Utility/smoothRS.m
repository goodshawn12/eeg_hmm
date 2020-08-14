function smooth_rs = smoothRS(rs, rt_latency, win_len_sec, Fs)    
    % smooth using 90 sec window
    facing_trial_num_list = zeros(size(rt_latency));
    smoothing_range_offset = -win_len_sec * Fs;
    offset_latency = rt_latency + smoothing_range_offset;
    for trial = 2:length(rt_latency)
        trial_facing_threshold = offset_latency(trial);
        facing_trial_num = 0;
        facing_trial_index = trial - facing_trial_num - 1;
        while trial_facing_threshold < rt_latency(facing_trial_index)
            facing_trial_num = facing_trial_num + 1;
            facing_trial_index = facing_trial_index - 1;
            if facing_trial_index < 1
                break;
            end
        end
        facing_trial_num_list(trial) = facing_trial_num;
    end
    
    smooth_rs = zeros(size(rs));
    for trial = 1:length(smooth_rs)
        trials_selection = (trial - facing_trial_num_list(trial)):trial;
        smooth_rs(trial) = median(rs(trials_selection));
    end    
end