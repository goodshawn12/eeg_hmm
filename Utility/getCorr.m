function [correlation] = getCorr(Gamma, K, Fs, rs, rt_latency, win_len_sec, smoothing_range_sec)
    win_offset = -win_len_sec * Fs;
    win_Gamma_mean = zeros(length(rt_latency), K);    
    for row = 1:length(win_Gamma_mean)
        latency = rt_latency(row);
        win_Gamma_mean(row,:) = mean(Gamma(latency+win_offset:latency,:),1);
    end
    
    % smooth using 90 sec window
    facing_trial_num_list = zeros(length(rt_latency), 1);
    smoothing_range_offset = (smoothing_range_sec - win_len_sec) * Fs;
    offset_latency = rt_latency - smoothing_range_offset;
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
    
    Gamma_median = zeros(length(rt_latency), K);
    for row = 1:length(Gamma_median)
        trials_selection = (row - facing_trial_num_list(row)):row;
        Gamma_median(row, :) = median(win_Gamma_mean(trials_selection,:), 1);
    end    
  
    correlation = zeros(1, K);
    for state = 1:K
        correlation(state) = corr(Gamma_median(:, state), rs');
    end
end