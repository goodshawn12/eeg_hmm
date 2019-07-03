cd '/home/ting/Documents/eeg_hmm';
addpath('/home/ting/Documents/eeglab')
eeglab;
addpath(genpath('HMM-MAR'))
addpath('/data/projects/Shawn/2019_HMM/data')
filename = 'session53_ASR20.set';
eegdata = pop_loadset(filename);

%% Colormap by epochs
n_epochs = 1;
if ~isempty(eegdata.epoch)
    n_epochs = size(eegdata.epoch, 2);
end

state_by_epoch = [];
if n_epochs > 1
    figure;
    state_by_epoch = (reshape(vpath, [], size(Gamma, 1)/n_epochs));
    imagesc(state_by_epoch);
    set(gca,'Title',text('String','State arranged by epochs'))
    colormap(parula)
else
    zero_padded_vpath = [zeros(training_data_size(1)-size(vpath,1),1); vpath];
    
    all_events = eegdata.event;
    epoch_start_offset = -2 * Fs;
    epoch_end_offset = 4 * Fs;
    state_by_epoch = zeros(length(all_events), epoch_end_offset - epoch_start_offset);
    
    rt = [];
    rt_off = [];
    for i = 1:length(eegdata.event)
        event = all_events(i);
        sliced_epoch_start = epoch_start_offset + event.latency;
        sliced_epoch_end = epoch_end_offset + event.latency -1;
        if (strcmp(event.type, '251') || strcmp(event.type, '252')) && ...
            sliced_epoch_start >= start_timepoint && ...
            sliced_epoch_end <= end_timepoint
        
            rt = [rt, (all_events(i+1).latency-all_events(i).latency)];
            rt_off = [rt_off, (all_events(i+2).latency-all_events(i).latency)];
            sliced_epoch = zero_padded_vpath(sliced_epoch_start:sliced_epoch_end);
            state_by_epoch(i,:) = sliced_epoch';
        end
    end
    state_by_epoch(~any(state_by_epoch,2),:) = [];
end

[sortedRT,sortIdx] = sort(rt);
[sortedRT_off,sortIdx_off] = sort(rt_off);


figure;
imagesc(state_by_epoch(sortIdx,:));
set(gca,'Title',text('String','State arranged by epochs'));
cmap = [0.1,0.1,1; 0.1,1,0.1; 1,0.1,0.1; 1,1,0.1];
colormap(cmap);

hold on,
plot(sortedRT-epoch_start_offset,1:length(rt),'linewidth',2,'color','w');
line([-epoch_start_offset, -epoch_start_offset], [1, length(rt)],'linewidth',2,'color','k');
plot(sortedRT_off-epoch_start_offset,1:length(rt_off),'linewidth',2,'color',[0.7,0.7,0.7]);

