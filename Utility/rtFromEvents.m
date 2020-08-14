function [rt, rt_off, rt_latency, rs, event_sparsity] = rtFromEvents(events, Fs, isClean)
    rt = zeros(1,length(events));
    rt_off = zeros(1,length(events));
    rt_latency = zeros(1,length(events));
    rs = zeros(1, length(events));
    event_sparsity = zeros(1, length(events));
    
    prev_offset = nan;
    for i = 1:length(events)
        event = events(i);
        event_type = event.type;
        if isnumeric(event_type)
            event_type = num2str(event_type);
        end      
        if strcmp(event_type, '251') || strcmp(event_type, '252')
            rt_latency(i) = event.latency;
            rt(i) = events(i+1).latency - event.latency;
            rt_off(i) = events(i+2).latency - event.latency;
            rs(i) = 1/(rt(i)/Fs); %rt_speed is in (1/seconds)
            
            if ~isnan(prev_offset)
                event_sparsity(i) = event.latency - prev_offset;
            end
            prev_offset = events(i+2).latency;
        end
    end
    
    nonempty_index = any(rt_latency, 1);
    rt = rt(nonempty_index);
    rt_off = rt_off(nonempty_index);
    rt_latency = rt_latency(nonempty_index);
    rs = rs(nonempty_index);
    event_sparsity = event_sparsity(nonempty_index);
    event_sparsity = event_sparsity(2:end); % The first 251 event has no prior LDT offset
    
    if isClean
        lower_remove_index = rt < 0.1 * Fs;
        upper_remove_index = rt > 10 * Fs;
        rt_off_remove_index = (rt_off - rt) < 0.1 * Fs;
        remove_index = lower_remove_index | upper_remove_index | rt_off_remove_index;
%         rt_removal_percentage_list(idx) = sum(remove_index) / length(remove_index);

        rt = rt(~remove_index);
        rt_off = rt_off(~remove_index);
        rt_latency = rt_latency(~remove_index);    
        rs = rs(~remove_index);

        [rs, speed_remove_index] = rmoutliers(rs);
        rt = rt(~speed_remove_index);
        rt_off = rt_off(~speed_remove_index);
        rt_latency = rt_latency(~speed_remove_index);
    end
end