cd '/home/ting/Documents/eeg_hmm';
addpath('/home/ting/Documents/eeg_hmm');
addpath(genpath('/home/ting/Documents/eeg_hmm/Utility'))
addpath('/data/projects/Shawn/2019_HMM/data/');
run('/data/common/matlab/eeglab/eeglab');

%% Specify number of inferred states
K = 3;
resultsDir = sprintf('/home/ting/Documents/eeg_hmm/AMICA_results/results_K%d/', K);

%% Manual setup files to process
data_filelist = dir(resultsDir);
filename_list = cell(length(data_filelist), 1);
for i = 1:length(data_filelist)
    file = data_filelist(i);
    if ~file.isdir || strcmp(file.name,'.') || strcmp(file.name,'..')
        continue
    end
    filename_list{i} = file.name;    
end
filename_list = filename_list(~cellfun('isempty', filename_list));
filename_list = unique(filename_list, 'stable');

n_of_files = length(filename_list);

%% Prepare raw data
n_epochs = 1;

session_name_list = cell(n_of_files,1);
K_list = cell(n_of_files,1);
epoch_length_list = cell(n_of_files,1);
Fs_list = cell(n_of_files,1);
v_list = cell(n_of_files,1);

for idx = 1:n_of_files
    fileDir = strcat(resultsDir, filename_list{idx});
    modout = loadmodout15(fileDir);
    
    session_name = filename_list{idx};
    session_name = split(session_name, '.');
    session_name = session_name{1};
    session_name = join(split(session_name, '_'));
    session_name_list{idx} = session_name{1};
    
    K_list{idx} = modout.num_models;
    v_list{idx} = modout.v; % v is K-by-pnts
end

results = table(filename_list, session_name_list, K_list, v_list);

%% Process Gamma and vpath
Gamma_list = cell(n_of_files,1);
vpath_list = cell(n_of_files,1);

for idx = 1:n_of_files
    v = results{idx, 'v_list'}{1};
    rej_index = any(sum(v, 1)==0, 1);
    valid_index = logical(1 - rej_index);
  
    v(:,valid_index) = 10.^v(:, valid_index);
    
    Gamma_list{idx} = v';
    
    [~, vpath] = max(v, [], 1);
    vpath(:,rej_index) = 0;
    vpath_list{idx} = vpath';
end

results = deleteVars(results, {'vpath_list', 'Gamma_list'});
results = addvars(results, vpath_list, Gamma_list);

%% Load EEG related metadata
events_list = cell(n_of_files, 1);
chanlocs_list = cell(n_of_files, 1);
Fs_list = cell(n_of_files,1);

% For all sessions
for idx = 1:n_of_files
    filename = strcat(filename_list{idx}, '.set');
    
    load('-mat', filename);
    events_list{idx} = EEG.event;
    chanlocs_list{idx} = extractfield(EEG.chanlocs, 'labels');
    Fs_list{idx} = EEG.srate;
end

results = deleteVars(results, {'events_list', 'chanlocs_list', 'Fs_list'});
results = addvars(results, events_list, chanlocs_list, Fs_list);

%% Generate rt, rt_off
rt_list = cell(n_of_files, 1);
rt_off_list = cell(n_of_files, 1);
rt_latency_list = cell(n_of_files, 1);
rt_speed_list = cell(n_of_files, 1);

for idx = 1:n_of_files
    Fs = results{idx, 'Fs_list'}{1};
    events = results{idx, 'events_list'}{1};  
    rt = zeros(1,length(events));
    rt_off = zeros(1,length(events));
    rt_latency = zeros(1,length(events));
    rt_speed = zeros(1, length(events));
    
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
            rt_speed(i) = 1./(rt(i)/Fs); %rt_speed is in (1/seconds)
        end
    end
    nonempty_index = any(rt_latency, 1);
    rt_list{idx} = rt(nonempty_index);
    rt_off_list{idx} = rt_off(nonempty_index);
    rt_latency_list{idx} = rt_latency(nonempty_index);
    rt_speed_list{idx} = rt_speed(nonempty_index);
end
results = deleteVars(results, {'rt_list', 'rt_off_list', 'rt_latency_list', 'rt_speed_list'});
results = addvars(results, rt_list, rt_off_list, rt_latency_list, rt_speed_list);

% RT and 1/RT Outlier Removal
rt_clean_list = cell(n_of_files, 1);
rt_off_clean_list = cell(n_of_files, 1);
rt_latency_clean_list = cell(n_of_files, 1);
rt_speed_clean_list = cell(n_of_files, 1); 

for idx = 1:n_of_files
    Fs = results{idx, 'Fs_list'}{1};
    rt = results{idx, 'rt_list'}{1};
    rt_off = results{idx, 'rt_off_list'}{1};
    rt_latency = results{idx, 'rt_latency_list'}{1}; 
    rt_speed = results{idx, 'rt_speed_list'}{1};
    
    lower_remove_index = rt < 0.1 * Fs;
    upper_remove_index = rt > 10 * Fs;
    remove_index = lower_remove_index | upper_remove_index;
    
    rt = rt(~remove_index);
    rt_off = rt_off(~remove_index);
    rt_latency = rt_latency(~remove_index);    
    rt_speed = rt_speed(~remove_index);
    
    [rt_speed, speed_remove_index] = rmoutliers(rt_speed);
    rt = rt(~speed_remove_index);
    rt_off = rt_off(~speed_remove_index);
    rt_latency = rt_latency(~speed_remove_index);
    
    rt_clean_list{idx} = rt;
    rt_off_clean_list{idx} = rt_off;
    rt_latency_clean_list{idx} = rt_latency;
    rt_speed_clean_list{idx} = rt_speed;
end

results = deleteVars(results, {'rt_clean_list', 'rt_off_clean_list', 'rt_latency_clean_list', 'rt_speed_clean_list'});
results = addvars(results, rt_clean_list, rt_off_clean_list, rt_latency_clean_list, rt_speed_clean_list);

%% Calculate state correlation with ERP
state_r_list = zeros(n_of_files, K);
trial_Gamma_mean_list = cell(n_of_files, 1);
win_len_sec = 2;

for idx = 1:n_of_files
    Fs = results{idx, 'Fs_list'}{1};
    Gamma = results{idx, 'Gamma_list'}{1};
    rt_speed = results{idx, 'rt_speed_clean_list'}{1};
    rt_latency = results{idx, 'rt_latency_clean_list'}{1};
    
    win_offset = -win_len_sec * Fs;
    win_Gamma_mean = zeros(length(rt_latency), K);    
    for row = 1:length(win_Gamma_mean)
        latency = rt_latency(row);
        win_Gamma_mean(row,:) = mean(Gamma(latency+win_offset:latency,:),1);
    end
  
    trial_Gamma_mean_list{idx} = win_Gamma_mean;
    for state = 1:K
        state_r_list(idx, state) = round(corr(win_Gamma_mean(:,state), rt_speed'), 3);
    end
end

% The sort is in ascending order so the first index in the permutation_list
% marks the lowest correlation (drowsy model) 
[~, permutation_list] = sort(state_r_list, 2);
inverse_permutation_list = zeros(size(permutation_list));
for row=1:size(permutation_list,1)
    inverse_permutation_list(row, permutation_list(row,:)) = 1:size(permutation_list,2);
end

results = deleteVars(results, {'trial_Gamma_mean_list', 'state_r_list', 'permutation_list', 'inverse_permutation_list'});
results = addvars(results, trial_Gamma_mean_list, state_r_list, permutation_list, inverse_permutation_list);

%% Smoothing vpath
delete(gcp('nocreate')); % shut down any current pool
npar = 27;
parpool(npar);   % request workers from the cluster

smoothing_window_len = 25;
smoothed_vpath_list = cell(n_of_files,1);
vpath_list = results.vpath_list;

parfor (idx = 1:n_of_files, npar)
    vpath = vpath_list{idx};
    smoothed_vpath_list{idx} = movingModeSmoothing(vpath, smoothing_window_len, 1, 1:K);
end

results = deleteVars(results, 'smoothed_vpath_list');
results = addvars(results, smoothed_vpath_list);

delete(gcp('nocreate'));

%% Check for rejection density
% for idx = 4:5
%     vpath = results{idx, 'vpath_list'}{1};
%     vpath(any(vpath~=0,1)) = 1;
%     imagesc(vpath);
% end

%% ######## Graphing Sections ########
%########################################

graphDir = sprintf('/home/ting/Documents/eeg_hmm/AMICA_results/graphs_K%d', K);
cd(graphDir);

% Preparation
if K == 5
    cmap = [1,0.1,0.1; 1,0.8,0.2; 0.8,0.7,0.55; 0.1,1,0.1; 0.1,0.1,1];
    state_description = {'Drowsy', 'Middle1', 'Middle2', 'Middle3', 'Alert'};
elseif K == 4
    cmap = [1,0.1,0.1; 1,0.8,0.2; 0.1,1,0.1; 0.1,0.1,1];
    state_description = {'Drowsy', 'Middle1', 'Middle2', 'Alert'};
elseif K == 3
    cmap = [1,0.1,0.1; 0.1,1,0.1; 0.1,0.1,1];
    state_description = {'Drowsy', 'Middle1', 'Alert'};
elseif K == 2
    cmap = [1,0.1,0.1; 0.1,0.1,1];
    state_description = {'Drowsy', 'Alert'};
else
    cmap = parula;
end

cmap_list = repmat(cmap, 1, 1, n_of_files);
for i = 1:n_of_files
    cmap_list(:,:,i) = cmap_list(results.inverse_permutation_list(i,:),:,i);
end
toSave = 0;
visible = 'off';

%% Colormap by epochs timelock at events
visible = 'on';
toSave = 0;

smooth = 1;
vertical_smoothing_window_len = 5;
isSorted = 1;

if smooth
    thin_line_width = 0.01;
    thick_line_width = 2; 
else
    thin_line_width = 0.01;
    thick_line_width = 0.01; 
end

for idx = 7:7
    cmap = cmap_list(:,:,idx);
    
    if smooth
        vpath = results{idx, 'smoothed_vpath_list'}{1};
    else
        vpath = results{idx, 'vpath_list'}{1};
        cmap = [1 1 1;cmap]; % Pad white for rejected points
    end
    
    num_mods = results{idx, 'K_list'}{1};
    all_events = results{idx, 'events_list'}{1};
    Fs = results{idx, 'Fs_list'}{1};
    rt = results{idx, 'rt_clean_list'}{1};
    rt_off = results{idx, 'rt_off_clean_list'}{1};
    rt_latency = results{idx, 'rt_latency_clean_list'}{1};
    filename = results{idx, 'filename_list'}{1};
    filename = split(filename, '.');
    filename = filename{1}
    
    training_data_size = [length(vpath) num_mods];
    start_timepoint = 1;
    end_timepoint = length(vpath);
    
    zero_padded_vpath = vpath;
    
    % Timelock 251/252, lane-departure task introduced
    epoch_start_offset_251 = -2 * Fs;
    epoch_end_offset_251 = 4 * Fs;

    % Timelock 253, response onset
    epoch_start_offset_253 = -3 * Fs;
    epoch_end_offset_253 = 3 * Fs;

    % Timelock 254, response off-set
    epoch_start_offset_254 = -4 * Fs;
    epoch_end_offset_254 = 2 * Fs;
    
    % Prevent epoch exceeding vpath dimension
    removeIndex = any(rt_latency + epoch_start_offset_251 < start_timepoint, 1)...
        | any(rt_latency + rt_off + epoch_end_offset_254 > end_timepoint, 1);    
    rt(removeIndex) = [];
    rt_off(removeIndex) = [];
    rt_latency(removeIndex) = [];

    
    [sortedRT, sortIdx] = sort(rt);
    sortedRT_off = rt_off(:, sortIdx);
    
    if ~isSorted
        sortedRT = rt;
        sortedRT_off = rt_off;
        sortIdx = 1:length(rt);
    end
    
    % Plot 251/252
    state_by_epoch = epochByEvent(zero_padded_vpath, rt_latency, epoch_start_offset_251, epoch_end_offset_251);    
    state_by_epoch = state_by_epoch(sortIdx,:);
    state_by_epoch = movingModeSmoothing(state_by_epoch, vertical_smoothing_window_len, 1, 1:K); 
    
    figure('Visible', visible);
    imagesc(state_by_epoch);
    set(gca,'Title',text('String','Time locked to lane-departure event'));
    set(gca,'XLabel',text('String','Time offset (sec)'));
    set(gca,'YLabel',text('String','Epochs'));
    xticks(1:Fs:Fs*6+1);
    xticklabels({'-2', '-1', '0', '1', '2', '3', '4'});
    colormap(cmap);

    hold on,
    line([-epoch_start_offset_251, -epoch_start_offset_251], [1, length(rt)], 'linewidth',2,'color','k');
    plot(sortedRT-epoch_start_offset_251, 1:length(rt), 'linewidth',thick_line_width,'color','w');
    plot(sortedRT_off-epoch_start_offset_251, 1:length(rt), 'linewidth',thin_line_width,'color', [0.8, 0.8, 0.8]);
  
    if toSave
        set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(filename,'_','251','.fig'))
        print(strcat(filename,'_','251'), '-djpeg')  
    end
    
    % Plot 253
    state_by_epoch = epochByEvent(zero_padded_vpath, rt_latency+rt, epoch_start_offset_253, epoch_end_offset_253);   
    state_by_epoch = state_by_epoch(sortIdx,:);
    state_by_epoch = movingModeSmoothing(state_by_epoch, vertical_smoothing_window_len, 1, 1:K); 

    figure('Visible', visible);
    imagesc(state_by_epoch);
    set(gca,'Title',text('String','Time locked to car driver response start'));
    set(gca,'XLabel',text('String','Time offset (sec)'));
    set(gca,'YLabel',text('String','Epochs'));
    xticks(1:Fs:Fs*6+1);
    xticklabels({'-3', '-2', '-1', '0', '1', '2', '3'});
    colormap(cmap);

    hold on,
    plot(-epoch_start_offset_253 - sortedRT, 1:length(rt), 'linewidth',thick_line_width,'color','k')
    line([-epoch_start_offset_253, -epoch_start_offset_253], [1, length(rt)],'linewidth',2,'color','w')
    plot(-epoch_start_offset_253 - sortedRT+sortedRT_off, 1:length(rt), 'linewidth',thin_line_width,'color',[0.8,0.8,0.8])

    if toSave
        set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(filename,'_','253','.fig'))
        print(strcat(filename,'_','253'), '-djpeg')
    end
       
    % Plot 254
    state_by_epoch = epochByEvent(zero_padded_vpath, rt_latency+rt_off, epoch_start_offset_254, epoch_end_offset_254);   
    state_by_epoch = state_by_epoch(sortIdx,:);
    state_by_epoch = movingModeSmoothing(state_by_epoch, vertical_smoothing_window_len, 1, 1:K); 

    figure('Visible', visible);
    imagesc(state_by_epoch);
    set(gca,'Title',text('String','Time locked to driver response end'));
    set(gca,'XLabel',text('String','Time offset (sec)'));
    set(gca,'YLabel',text('String','Epochs'));
    xticks(1:Fs:Fs*6+1);
    xticklabels({'-4', '-3', '-2', '-1', '0', '1', '2'});
    colormap(cmap);

    hold on,
    plot(-epoch_start_offset_254 - sortedRT_off, 1:length(rt), 'linewidth',thin_line_width,'color','k')
    plot(-epoch_start_offset_254 - sortedRT_off+sortedRT, 1:length(rt), 'linewidth',thin_line_width,'color','w')
    line([-epoch_start_offset_254, -epoch_start_offset_254], [1, length(rt)], 'linewidth',2,'color',[0.8,0.8,0.8])
    
    if toSave
        set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(filename,'_','254','.fig'))
        print(strcat(filename,'_','254'), '-djpeg')
    end
end

%% Plot smoothed state time
visible = 'off';
toSave = 1;
for idx = 1:n_of_files
    Fs = results{idx, 'Fs_list'}{1};
    Gamma = results{idx, 'Gamma_list'}{1};
    rt = results{idx, 'rt_clean_list'}{1};
    rt_latency = results{idx, 'rt_latency_clean_list'}{1};
    rt_speed = results{idx, 'rt_speed_clean_list'}{1};
    state_r_list = results{idx, 'state_r_list'};
    cmap = cmap_list(:,:,idx);
    filename = results{idx, 'filename_list'}{1};
    filename = split(filename, '.');
    filename = filename{1}
    
    winLen = 30 * Fs;
    walkLen = 5 * Fs;
    [nTime,nStates] = size(Gamma);
    stateProb = zeros(nStates,floor(nTime/walkLen));
    for iter = 1:floor((nTime-winLen+walkLen)/walkLen)
        time_range = (iter-1)*walkLen+1:(iter-1)*walkLen+winLen;
        stateProb(:,iter) = mean(Gamma(time_range,:),1)';
    end
    end_minutes = (size(stateProb,2)-1)*walkLen/Fs/60;
    
    % Check if rejection points are significant
    error_rate = sum(any(stateProb==zeros(K,1), 1))/size(stateProb, 2);
    if  error_rate > 0.02
        warning('Over 2% of state path points are biased due to rejection')
        sprintf('Percentage of rejected points: %d', error_rate)
    end
    
    n_subplots = K + 2;
    figure('Visible', visible);
%     imagesc([0 end_minutes], [1 size(stateProb,1)], stateProb), colorbar;
%     colormap(gamma_axis, hot);
%     title('State Probability');
%     ylabel('State')
        
    [~, corresponding_vpath] = max(stateProb, [], 1);
    vpath_axis = subplot(n_subplots,1,1);
    imagesc([0 end_minutes], [1 1], corresponding_vpath);
    
    % Configure colorbar
    ytick_space = (K-1)/2/K;
    yticks = 1-ytick_space:2*ytick_space:K-0.01;
    yticks(1) = [];
    ytickslabel = cell(1,K);
    for state = 1:K
        ytickslabel{state} = num2str(state);
    end
    colorbar('YTick', yticks, 'YTickLabel', ytickslabel);
    colormap(vpath_axis, cmap)
    title('Viterbi path')
    yticks([]);

    rs_axis = subplot(n_subplots,1,2); plot(rt_latency/Fs/60, movmean(rt_speed, 7), 'k');
    title('Reaction Speed')
    ylabel('Speed (1/sec)')
    xlim([0, vpath_axis.XLim(2)*rs_axis.Position(3)/vpath_axis.Position(3)])
    
    for state = 1:K
        state_axis = subplot(n_subplots,1,state+2);
        plot(1:(end_minutes-1)/(size(stateProb,2)-1):end_minutes, movmean(stateProb(state,:),7), 'Color', cmap(state,:))
        r = state_r_list(state);
        title(strcat('State ', num2str(state)))
        ylabel('Probability')
        txt = ['r = ' num2str(r)];
        text(state_axis.Position(3), state_axis.Position(4), txt, 'FontSize', 10)
        xlim([0, vpath_axis.XLim(2)*state_axis.Position(3)/vpath_axis.Position(3)])
    end
    xlabel('Time (minutes)')
    
    if toSave
        output_filename = strcat(filename,'_','path');
        set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(output_filename,'.fig'))
        print(output_filename, '-djpeg')
    end
end
