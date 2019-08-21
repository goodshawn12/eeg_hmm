addpath('/home/ting/Documents/eeglab');
addpath('/home/ting/Documents/eeg_hmm')
addpath(genpath('HMM-MAR'));
addpath('/data/projects/Shawn/2016 JNE/dataset/');
%% Manual setup
cd '/home/ting/Documents/eeg_hmm/results';

filename_list = {'GAU_session16_ASR20_0.mat',
    'GAU_session30_ASR20_0.mat',
    'GAU_session34_ASR20_0.mat',
    'GAU_session36_ASR20_0.mat',
    'GAU_session38_ASR20_0.mat',
    'GAU_session44_ASR20_0.mat',
    'GAU_session48_ASR20_0.mat',
    'GAU_session53_ASR20_0.mat',
    'GAU_session60_ASR20_0.mat',
    'GAU_session62_ASR20_0.mat',
    'MAR_session16_ASR20_0.mat',
    'MAR_session30_ASR20_0.mat',
    'MAR_session34_ASR20_0.mat',
    'MAR_session36_ASR20_0.mat',
    'MAR_session38_ASR20_0.mat',
    'MAR_session44_ASR20_0.mat',
    'MAR_session48_ASR20_0.mat',
    'MAR_session53_ASR20_0.mat',
    'MAR_session60_ASR20_0.mat',
    'MAR_session62_ASR20_0.mat',
    };
filename_list = {'MAR_session16_ASR20_0.mat'};
n_of_files = length(filename_list);
K = 3;

%% Prepare raw data
Fs = 250;
n_epochs = 1;

K_list = cell(n_of_files,1);
epoch_length_list = cell(n_of_files,1);
timepoint_start_list = cell(n_of_files,1);
timepoint_end_list = cell(n_of_files,1);
hmm_train_list = cell(n_of_files,1);
Gamma_list = cell(n_of_files,1);
vpath_list = cell(n_of_files,1);

for i = 1:n_of_files
    load(filename_list{i});  
    
    K_list{i,1} = hmm.K;
    epoch_length_list{i,1} = training_data_size(1);
    if select_start == 0
        timepoint_start_list{i,1} = 1;
    else
        timepoint_start_list{i,1} = floor(training_data_size(1) * select_start);
    end
    timepoint_end_list{i,1} = floor(training_data_size(1) * select_end);
    hmm_train_list{i,1} = hmm.train;
    Gamma_list{i,1} = Gamma;
    vpath_list{i,1} = vpath;
end

results = table(filename_list, K_list, epoch_length_list, timepoint_start_list, ...
    timepoint_end_list, hmm_train_list, Gamma_list, vpath_list);

%% Smoothing
smoothing_window_length = 25;
smoothed_vpath_list = cell(n_of_files,1);

for idx = 1:n_of_files
    smoothed_vpath = results{idx, 'vpath_list'}{1};
    
    last_window_index = ceil(length(smoothed_vpath)/smoothing_window_length);
    for i = 1:last_window_index
        if i < last_window_index
           window_start = (i-1)*smoothing_window_length + 1;
           window_end = i*smoothing_window_length;
           smoothed_vpath(window_start:window_end) = mode(smoothed_vpath(window_start:window_end));
        elseif i == last_window_index
           window_start = (i-1)*smoothing_window_length + 1;
           window_end = length(smoothed_vpath);
           smoothed_vpath(window_start:window_end) = mode(smoothed_vpath(window_start:window_end));
        end
    end
    smoothed_vpath_list{idx,1} = smoothed_vpath;
end

% smoothing_window_length = 10;
% smoothing_stride = 2;
% smoothed_vpath_list = cell(n_of_files,1);
% 
% for idx = 1:n_of_files
%     vpath = results{1, 'vpath_list'}{1};
%     smoothed_vpath = vpath;
%     
%     i = floor(smoothing_window_length/2)+1;
%     while i <= length(vpath)-ceil(smoothing_window_length/2)
%        window_start = i - floor(smoothing_window_length/2);
%        window_end = i + ceil(smoothing_window_length/2) - 1;
%        smoothed_vpath(i) = mode(vpath(window_start:window_end));
%        i = i + smoothing_stride;
%     end
%     smoothed_vpath_list{idx,1} = smoothed_vpath;
% end

results = addvars(results, smoothed_vpath_list);

%% Load events
events_list = cell(n_of_files, 1);
eeglab;
for i = 1:n_of_files
    filename = split(filename_list{i}, '_');
    filename = strcat(filename{2}, '_', filename{3}, '.set');
    eeg_data = pop_loadset(filename);
    
    events_list{i} = eeg_data.event;
end

results = addvars(results, events_list);
%% Prepare for graphing
cd '/home/ting/Documents/eeg_hmm/figures';

if K == 4
    cmap = [0.1,0.1,1; 0.1,1,0.1; 1,0.1,0.1; 1,1,0.1];
elseif K == 3
    cmap = [0.1,0.1,1; 0.1,1,0.1; 1,0.1,0.1];
elseif K == 2
    cmap = [0.1,0.1,1; 0.1,1,0.1];
else
    cmap = parula;
end

visible = 'off';
%% Colormap by epochs timelock at events
rt_list = cell(n_of_files, 1);
rt_off_list = cell(n_of_files, 1);

for idx = 1:n_of_files
    smoothed_vpath = results{idx, 'smoothed_vpath_list'}{1};
    training_data_size = results{idx, 'epoch_length_list'}{1};
    start_timepoint = results{idx, 'timepoint_start_list'}{1};
    end_timepoint = results{idx, 'timepoint_end_list'}{1};
    all_events = results{idx, 'events_list'}{1};
    filename = results{idx, 'filename_list'}{1};
    filename = split(filename, '.');
    filename = filename{1}
    zero_padded_vpath = [zeros(training_data_size(1)-size(smoothed_vpath,1),1); smoothed_vpath];
    
    % Timelock 251/252, lane-departure task introduced
    epoch_start_offset = -2 * Fs;
    epoch_end_offset = 4 * Fs;
    rt = [];
    rt_off = [];
    state_by_epoch = zeros(length(all_events), epoch_end_offset - epoch_start_offset); 
    counter = 0;
    for i = 1:length(all_events)
        event = all_events(i);
        if (strcmp(event.type, '251') || strcmp(event.type, '252'))           
            sliced_epoch_start = event.latency + epoch_start_offset;
            sliced_epoch_end = event.latency + epoch_end_offset - 1;
            if sliced_epoch_start >= start_timepoint && sliced_epoch_end <= end_timepoint
                rt = [rt, (all_events(i+1).latency-all_events(i).latency)];
                rt_off = [rt_off, (all_events(i+2).latency-all_events(i).latency)];
                sliced_epoch = zero_padded_vpath(sliced_epoch_start:sliced_epoch_end);
                state_by_epoch(i,:) = sliced_epoch';
                counter = counter + 1;
            end                   
        end
    end
    
    state_by_epoch(~any(state_by_epoch,2),:) = [];
    [sortedRT, sortIdx] = sort(rt);
    sortedRT_off = rt_off(:, sortIdx);
    state_by_epoch = state_by_epoch(sortIdx,:);
    state_by_epoch = verticalSmoothing(state_by_epoch, 3);

    rt_list{idx} = rt;
    rt_off_list{idx} = rt_off;   
    
    figure('Visible', visible);
    imagesc(state_by_epoch);
    set(gca,'Title',text('String','Time locked to lane-departure event'));
    set(gca,'XLabel',text('String','Timepoints in an epoch'));
    set(gca,'YLabel',text('String','Epochs'));
    colormap(cmap);

    hold on,
    line([-epoch_start_offset, -epoch_start_offset], [1, length(rt)], 'linewidth',2,'color','k');
    plot(sortedRT-epoch_start_offset, 1:length(rt), 'linewidth',2,'color','w');
    plot(sortedRT_off-epoch_start_offset, 1:length(rt), 'linewidth',0.01,'color', [0.8, 0.8, 0.8]);
  
    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf, strcat(filename,'_','251','.fig'))
    print(strcat(filename,'_','251'), '-djpeg')    
    
    % Timelock 253, response onset
    epoch_start_offset = -3 * Fs;
    epoch_end_offset = 3 * Fs;
    rt = [];
    rt_off = [];
    state_by_epoch = zeros(length(all_events), epoch_end_offset - epoch_start_offset);   
    for i = 1:length(all_events)
        event = all_events(i);
        if (strcmp(event.type, '253'))          
            sliced_epoch_start = event.latency + epoch_start_offset;
            sliced_epoch_end = event.latency + epoch_end_offset - 1;           
            if sliced_epoch_start >= start_timepoint && sliced_epoch_end <= end_timepoint
                rt = [rt, (all_events(i).latency-all_events(i-1).latency)];
                rt_off = [rt_off, (all_events(i+1).latency-all_events(i-1).latency)];
                sliced_epoch = zero_padded_vpath(sliced_epoch_start:sliced_epoch_end);
                state_by_epoch(i,:) = sliced_epoch';
            end                   
        end
    end
    state_by_epoch(~any(state_by_epoch,2),:) = [];
    [sortedRT, sortIdx] = sort(rt);
    sortedRT_off = rt_off(:, sortIdx);
    state_by_epoch = state_by_epoch(sortIdx,:);
    state_by_epoch = verticalSmoothing(state_by_epoch, 3);

    figure('Visible', visible);
    imagesc(state_by_epoch);
    set(gca,'Title',text('String','Time locked to car driver response start'));
    set(gca,'XLabel',text('String','Timepoints in an epoch'));
    set(gca,'YLabel',text('String','Epochs'));
    colormap(cmap);

    hold on,
    plot(-epoch_start_offset-sortedRT, 1:length(rt), 'linewidth',2,'color','k')
    line([-epoch_start_offset, -epoch_start_offset], [1, length(rt)],'linewidth',2,'color','w')
    plot(-epoch_start_offset-sortedRT+sortedRT_off, 1:length(rt), 'linewidth',0.01,'color',[0.8,0.8,0.8])

    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf, strcat(filename,'_','253','.fig'))
    print(strcat(filename,'_','253'), '-djpeg')
    
    % Timelock 254, response off-set
    epoch_start_offset = -4 * Fs;
    epoch_end_offset = 2 * Fs;
    rt = [];
    rt_off = [];
    state_by_epoch = zeros(length(all_events), epoch_end_offset - epoch_start_offset);   
    for i = 1:length(all_events)
        event = all_events(i);
        if (strcmp(event.type, '254'))     
            sliced_epoch_start = event.latency + epoch_start_offset;
            sliced_epoch_end = event.latency + epoch_end_offset - 1;           
            if sliced_epoch_start >= start_timepoint && sliced_epoch_end <= end_timepoint
                rt = [rt, (all_events(i-1).latency-all_events(i-2).latency)];
                rt_off = [rt_off, (all_events(i).latency-all_events(i-2).latency)];
                sliced_epoch = zero_padded_vpath(sliced_epoch_start:sliced_epoch_end);
                state_by_epoch(i,:) = sliced_epoch';
            end                   
        end
    end
    state_by_epoch(~any(state_by_epoch,2),:) = [];
    [sortedRT, sortIdx] = sort(rt);
    sortedRT_off = rt_off(:, sortIdx);
    state_by_epoch = state_by_epoch(sortIdx,:);
    state_by_epoch = verticalSmoothing(state_by_epoch, 3);

    figure('Visible', visible);
    imagesc(state_by_epoch);
    set(gca,'Title',text('String','Time locked to driver response end'));
    set(gca,'XLabel',text('String','Timepoints in an epoch'));
    set(gca,'YLabel',text('String','Epochs'));
    colormap(cmap);

    hold on,
    plot(-epoch_start_offset-sortedRT_off, 1:length(rt), 'linewidth',0.01,'color','k')
    plot(-epoch_start_offset-sortedRT_off+sortedRT, 1:length(rt), 'linewidth',0.01,'color','w')
    line([-epoch_start_offset, -epoch_start_offset], [1, length(rt)], 'linewidth',2,'color',[0.8,0.8,0.8])
    
    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf, strcat(filename,'_','254','.fig'))
    print(strcat(filename,'_','254'), '-djpeg')
end

results = addvars(results, rt_list);
%% Plot smoothed state time
for idx = 1:n_of_files
    winLen = 30 * Fs;
    walkLen = 5 * Fs;
    Gamma = results{idx, 'Gamma_list'}{1};
    rt = results{idx, 'rt_list'}{1};
    filename = results{idx, 'filename_list'}{1};
    filename = split(filename, '.');
    filename = filename{1}
    [nTime,nMod] = size(Gamma);

    stateProb = zeros(nMod,floor(nTime/walkLen));
    for it = 1:floor((nTime-winLen+walkLen)/walkLen)
        time_range = (it-1)*walkLen+1:(it-1)*walkLen+winLen;
        stateProb(:,it) = mean(Gamma(time_range,:),1)';
    end
    figure('Visible', visible);
    gamma_axis = subplot(3,1,1); 
    imagesc(stateProb), colorbar,
    colormap(gamma_axis, hot);
    title('State time courses')

    [~, corresponding_vpath] = max(stateProb, [], 1);
    vpath_axis = subplot(3,1,2);
    imagesc(corresponding_vpath), colorbar,
    colormap(vpath_axis, cmap)
    title('Viterbi path')

    subplot(3,1,3), plot(log10(rt))
    title('Reaction Time at Log10')
    
    saveas(gcf, strcat(filename,'_','path','.fig'))
    print(strcat(filename,'_','path'), '-djpeg')
end

%% Plotting transition probability matrix
% figure;
% imagesc(hmm.P)
