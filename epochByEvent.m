function epoched_data = epochByEvent(data, event_latency, epoch_start_offset, epoch_end_offset)
    epoched_data = zeros(length(event_latency), epoch_end_offset - epoch_start_offset + 1); 
    for i = 1:length(event_latency)
        latency = event_latency(i);
        sliced_epoch_start = latency + epoch_start_offset;
        sliced_epoch_end = latency + epoch_end_offset;        
        sliced_epoch = data(sliced_epoch_start:sliced_epoch_end);
        epoched_data(i,:) = sliced_epoch';
    end    
    epoched_data(~any(epoched_data,2),:) = []; 
end