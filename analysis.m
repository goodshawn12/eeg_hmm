addpath('/home/ting/Documents/eeglab');
addpath('/home/ting/Documents/eeg_hmm');
addpath(genpath('/home/ting/Documents/eeg_hmm/Utility'));
addpath(genpath('/home/ting/Documents/eeg_hmm/HMM-MAR'));
addpath('/data/projects/Shawn/2016 JNE/dataset/');
addpath('/data/projects/Shawn/2019_HMM/data/');

%% ### This pipeline only supports results with the same K (number of inferred states)
cd '/home/ting/Documents/eeg_hmm/results_K4';
K = 4;

%% Manual setup files to process
% Selected sessions
% data_filelist = [dir('GAU*m*1.mat');dir('GAU*ASR20*0.mat');dir('MAR*ASR20*0.mat')];
% All session data
data_filelist = dir('GAU*.mat');

filename_list = cell(length(data_filelist), 1);
for i = 1:length(data_filelist)
    file = data_filelist(i);
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
timepoint_start_list = cell(n_of_files,1);
timepoint_end_list = cell(n_of_files,1);
training_data_size_list = zeros(n_of_files,2);
hmm_list = cell(n_of_files,1);
Fs_list = cell(n_of_files,1);
Gamma_list = cell(n_of_files,1);
vpath_list = cell(n_of_files,1);
vars_to_load = {'hmm', 'training_data_size', 'select_start', 'select_end', 'Gamma', 'vpath'};

for idx = 1:n_of_files
    load(filename_list{idx}, vars_to_load{:});  
    
    session_name = filename_list{idx};
    session_name = split(session_name, '.');
    session_name = session_name{1};
    session_name = join(split(session_name, '_'));
    
    session_name_list{idx} = session_name{1};
    K_list{idx} = hmm.K;
    epoch_length_list{idx} = training_data_size(1);
    training_data_size_list(idx,:) = training_data_size;
    if select_start == 0
        timepoint_start_list{idx} = 1;
    else
        timepoint_start_list{idx} = floor(training_data_size(1) * select_start);
    end
    timepoint_end_list{idx} = floor(training_data_size(1) * select_end);
    hmm_list{idx} = hmm;
    Fs_list{idx} = hmm.train.Fs;
    Gamma_list{idx} = Gamma;
    vpath_list{idx} = vpath;
end

results = table(filename_list, session_name_list, K_list, epoch_length_list, training_data_size_list, ...
    timepoint_start_list, timepoint_end_list, Fs_list, hmm_list, Gamma_list, vpath_list);

%% Load EEG related metadata
events_list = cell(n_of_files, 1);
chanlocs_list = cell(n_of_files, 1);
% For selected sessions
% for idx = 1:n_of_files
%     filename = split(filename_list{idx}, '_');
%     if ~strcmp(filename{3}, 'ASR20')
%         filename{2} = strcat('s', extractAfter(filename{2}, 'session'));
%     end
%         
%     filename = strcat( filename{2}, '_', filename{3}, '.set')    
%     load('-mat', filename); 
%     events_list{idx} = EEG.event;
% end

% For all sessions
for idx = 1:n_of_files
    filename = split(filename_list{idx}, '_');
    filename = join(filename(2:length(filename)-1), '_');
    filename = strcat(filename, '.set');
    
    load('-mat', filename{1});
    events_list{idx} = EEG.event;
    chanlocs_list{idx} = extractfield(EEG.chanlocs, 'labels');
end

results = deleteVars(results, {'events_list', 'chanlocs_list'});
results = addvars(results, events_list, chanlocs_list);

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

%% RT and 1/RT Outlier Removal
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
%% Get RS statistics
rs_mean_list = zeros(n_of_files, 1);
rs_st_list = zeros(n_of_files, 1);
for idx = 1:n_of_files
    rs = results{idx, 'rt_speed_clean_list'}{1};
    rs_mean_list(idx) = mean(rs);
    rs_st_list(idx) = std(rs);
end

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

%% Dominant model RT 
% window_time_sec = 2;
% ttest_h_list = zeros(n_of_files, 3);
% ttest_p_list = zeros(n_of_files, 3);
% dominant_state_rt_list = cell(n_of_files, 3);
% 
% for idx = 1:n_of_files
%     Fs = results{idx, 'Fs_list'}{1};
%     dominant_model_window_length = Fs * window_time_sec;
%     rt = results{idx, 'rt_clean_list'}{1};
%     rt_latency = results{idx, 'rt_latency_clean_list'}{1};    
%     vpath = results{idx, 'vpath_list'}{1};
%     training_data_size = results{idx, 'epoch_length_list'}{1};
%     zero_padded_vpath = [zeros(training_data_size(1)-size(vpath,1),1); vpath];
%     
%     dominant_state = zeros(1, length(rt_latency));
%     for i = 1:length(rt_latency)
%         window_end = rt_latency(i);
%         window_start = rt_latency(i) - dominant_model_window_length + 1;
%         dominant_state(i) = quickMode(zero_padded_vpath(window_start:window_end,:), [1 2 3]);
%     end
%     
%     dominant1 = rt(any(dominant_state==1, 1));
%     dominant2 = rt(any(dominant_state==2, 1));
%     dominant3 = rt(any(dominant_state==3, 1));
%     dominant_state_rt_list(idx, :) = {dominant1, dominant2, dominant3};
%     
%     [h12, p12] = ttest2(dominant1, dominant2);
%     [h13, p13] = ttest2(dominant1, dominant3);
%     [h23, p23] = ttest2(dominant2, dominant3);
%     ttest_h_list(idx, :) = [h12, h13, h23];
%     ttest_p_list(idx, :) = [p12, p13, p23];
% end
% 
% results = deleteVars(results, {'dominant_state_rt_list', 'ttest_p_list'});
% results = addvars(results, dominant_state_rt_list, ttest_p_list);

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

%% Calculate overall state probability
overall_state_prob_list = zeros(n_of_files, K);
for idx = 1:n_of_files
    vpath = results{idx, 'vpath_list'}{1};
    
    for state = 1:K
        overall_state_prob_list(idx, state) = sum(vpath==state)/length(vpath);
    end
end

results = deleteVars(results, 'overall_state_prob_list');
results = addvars(results, overall_state_prob_list);

%% ######## Graphing Sections ########
%########################################
% Preparation
cd 'graphs';

if K == 4
    cmap = [1,0.1,0.1; 1,1,0.1; 0.1,1,0.1; 0.1,0.1,1];
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

%% Box Plot for Reaction time vs dominant model
% dominatn_state_rt_list = results.dominant_state_rt_list;
% 
% for idx = 1:n_of_files    
%     filename = results{idx, 'filename_list'}{1};
%     filename = split(filename, '.');
%     filename = split(filename{1}, '_');
%     filename = join(filename, ' ');
%     
%     dominant_list = dominant_state_rt_list(idx,:);
%     dominant1 = dominant_list{1}; map1 = repmat({'1'}, length(dominant1), 1);
%     dominant2 = dominant_list{2}; map2 = repmat({'2'}, length(dominant2), 1);
%     dominant3 = dominant_list{3}; map3 = repmat({'3'}, length(dominant3), 1);
%     
%     x = [dominant1'; dominant2'; dominant3']./Fs;
%     map = [map1; map2; map3];
%     
%     subplot(4, 5, idx);
%     boxplot(x, map, 'DataLim', [0, 3]);
%     title(filename);
%     
%     % In order to get figure wise ylabel and xlabel location
%     if idx == 1
%         loc =  get(gca, 'position');
%         yhigh = loc(2) + loc(4);
%     end
%     if idx == 16
%        loc =  get(gca, 'position');
%        xlow = loc(1);
%        ylow = loc(2); 
%     end
%     if idx == n_of_files
%         loc = get(gca, 'position');
%         xhigh = loc(1) + loc(3);
%     end
% end
% 
% ax = axes('position', [xlow, ylow, xhigh-xlow, yhigh-ylow], 'Visible', 'off');
% ax.XLabel.Visible = 'on'; ax.YLabel.Visible = 'on';
% xlabel(strcat('Dominant state', ' ', num2str(window_time_sec), ' seconds before lane-departure event'));
% ylabel('RT (sec)');

% saveas(gcf, strcat('20_RT_vs_State_box_', num2str(window_time_sec), 'sec.fig'));

%% Colormap by epochs timelock at events
visible = 'on';
toSave = 0;
raw = 1;
if raw
    smoothing_window_len = 1;
    thin_line_width = 0.01;
    thick_line_width = 0.01;
else
    smoothing_window_len = 5;
    thin_line_width = 0.01;
    thick_line_width = 2;
end

% Timelock 251/252, lane-departure task introduced
epoch_start_offset_251 = -2 * Fs;
epoch_end_offset_251 = 4 * Fs;

% Timelock 253, response onset
epoch_start_offset_253 = -3 * Fs;
epoch_end_offset_253 = 3 * Fs;

% Timelock 254, response off-set
epoch_start_offset_254 = -4 * Fs;
epoch_end_offset_254 = 2 * Fs;

for idx = 24:24
    if raw
        smoothed_vpath = results{idx, 'vpath_list'}{1};
    else
        smoothed_vpath = results{idx, 'smoothed_vpath_list'}{1};
    end    
    training_data_size = results{idx, 'epoch_length_list'}{1};
    start_timepoint = results{idx, 'timepoint_start_list'}{1};
    end_timepoint = results{idx, 'timepoint_end_list'}{1};
    all_events = results{idx, 'events_list'}{1};
    Fs = results{idx, 'Fs_list'}{1};
    rt = results{idx, 'rt_clean_list'}{1};
    rt_off = results{idx, 'rt_off_clean_list'}{1};
    rt_latency = results{idx, 'rt_latency_clean_list'}{1};
    filename = results{idx, 'filename_list'}{1};
    filename = split(filename, '.');
    filename = filename{1}
    
    zero_padded_vpath = [zeros(training_data_size(1)-size(smoothed_vpath,1),1); smoothed_vpath];
    cmap = cmap_list(:,:,idx);
    
    % Prevent epoch exceeding vpath dimension
    removeIndex = any(rt_latency + epoch_start_offset_251 < start_timepoint, 1)...
        | any(rt_latency + rt_off + epoch_end_offset_254 > end_timepoint, 1);    
    rt(removeIndex) = [];
    rt_off(removeIndex) = [];
    rt_latency(removeIndex) = [];

    [sortedRT, sortIdx] = sort(rt);
    sortedRT_off = rt_off(:, sortIdx);
    
    if raw
        sortedRT = rt;
        sortedRT_off = rt_off;
        sortIdx = 1:length(rt);
    end
    
    % Plot 251/252
    state_by_epoch = epochByEvent(zero_padded_vpath, rt_latency, epoch_start_offset_251, epoch_end_offset_251);    
    state_by_epoch = state_by_epoch(sortIdx,:);
    state_by_epoch = movingModeSmoothing(state_by_epoch, smoothing_window_len, 1, 1:K); 
    
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
        saveas(gcf, strcat(filename,'_','251','.fig'))
        print(strcat(filename,'_','251'), '-djpeg')  
    end
    
    % Plot 253
    state_by_epoch = epochByEvent(zero_padded_vpath, rt_latency+rt, epoch_start_offset_253, epoch_end_offset_253);   
    state_by_epoch = state_by_epoch(sortIdx,:);
    state_by_epoch = movingModeSmoothing(state_by_epoch, smoothing_window_len, 1, 1:K); 

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
        saveas(gcf, strcat(filename,'_','253','.fig'))
        print(strcat(filename,'_','253'), '-djpeg')
    end
       
    % Plot 254
    state_by_epoch = epochByEvent(zero_padded_vpath, rt_latency+rt_off, epoch_start_offset_254, epoch_end_offset_254);   
    state_by_epoch = state_by_epoch(sortIdx,:);
    state_by_epoch = movingModeSmoothing(state_by_epoch, smoothing_window_len, 1, 1:K); 

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
        saveas(gcf, strcat(filename,'_','254','.fig'))
        print(strcat(filename,'_','254'), '-djpeg')
    end
end

%% Prepare to plot all session timelock
all_session_251 = [];
all_session_253 = [];
all_session_254 = [];
all_session_RT = [];
all_session_RT_off = [];
all_session_RT_latency = [];

for idx = 1:n_of_files
    smoothed_vpath = results{idx, 'smoothed_vpath_list'}{1};
    training_data_size = results{idx, 'epoch_length_list'}{1};
    start_timepoint = results{idx, 'timepoint_start_list'}{1};
    end_timepoint = results{idx, 'timepoint_end_list'}{1};
    all_events = results{idx, 'events_list'}{1};
    Fs = results{idx, 'Fs_list'}{1};
    rt = results{idx, 'rt_clean_list'}{1};
    rt_off = results{idx, 'rt_off_clean_list'}{1};
    rt_latency = results{idx, 'rt_latency_clean_list'}{1};
    permutation = results{idx, 'permutation_list'};
    filename = results{idx, 'filename_list'}{1}
    
    zero_padded_vpath = [zeros(training_data_size(1)-size(smoothed_vpath,1),1); smoothed_vpath];
    ordered_vpath = zero_padded_vpath;
    ordered_vpath(zero_padded_vpath == permutation(1)) = 1;
    ordered_vpath(zero_padded_vpath == permutation(2)) = 2;
    ordered_vpath(zero_padded_vpath == permutation(3)) = 3;
        
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
    sortedRT_latency = rt_latency(:, sortIdx);
    
    all_session_RT = cat(2, all_session_RT, rt);
    all_session_RT_latency = cat(2, all_session_RT_latency, sortedRT_latency);
    all_session_RT_off = cat(2, all_session_RT_off, sortedRT_off);
    
    % Add 251/252
    state_by_epoch = epochByEvent(ordered_vpath, rt_latency, epoch_start_offset_251, epoch_end_offset_251);   
    state_by_epoch = state_by_epoch(sortIdx,:);
    all_session_251 = cat(1, all_session_251, state_by_epoch); 
    
    % Add 253
    state_by_epoch = epochByEvent(ordered_vpath, rt_latency+rt, epoch_start_offset_253, epoch_end_offset_253);   
    state_by_epoch = state_by_epoch(sortIdx,:);
    all_session_253 = cat(1, all_session_253, state_by_epoch); 
    
    % Add 254
    state_by_epoch = epochByEvent(ordered_vpath, rt_latency+rt, epoch_start_offset_253, epoch_end_offset_253);   
    state_by_epoch = state_by_epoch(sortIdx,:);
    all_session_254 = cat(1, all_session_254, state_by_epoch);   
end

[sortedRT, sortIdx] = sort(all_session_RT);
sortedRT_off = all_session_RT_off(:, sortIdx);
sortedRT_latency = all_session_RT_latency(:, sortIdx);
all_session_251 = all_session_251(sortIdx, :);
all_session_253 = all_session_253(sortIdx, :);
all_session_254 = all_session_254(sortIdx, :);

%% Plot all session timelock
visible = 'on';
toSave = 0;
cmap = [1,0.1,0.1; 0.1,1,0.1; 0.1,0.1,1];

figure('Visible', visible);
imagesc(all_session_251);
set(gca,'Title',text('String','Time locked to lane-departure event'));
set(gca,'XLabel',text('String','Time offset (sec)'));
set(gca,'YLabel',text('String','Epochs'));
xticks(1:Fs:Fs*6+1);
xticklabels({'-2', '-1', '0', '1', '2', '3', '4'});
colormap(cmap);

hold on,
line([-epoch_start_offset_251, -epoch_start_offset_251], [1, length(sortedRT)], 'linewidth',2,'color','k');
plot(sortedRT-epoch_start_offset_251, 1:length(sortedRT), 'linewidth',2,'color','w');
plot(sortedRT_off-epoch_start_offset_251, 1:length(sortedRT), 'linewidth',0.01,'color', [0.8, 0.8, 0.8]);

if toSave
    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf, strcat('all_session_251','.fig'))
    print('all_session_251', '-djpeg')  
end

% Plot 253
figure('Visible', visible);
imagesc(all_session_253);
set(gca,'Title',text('String','Time locked to car driver response start'));
set(gca,'XLabel',text('String','Time offset (sec)'));
set(gca,'YLabel',text('String','Epochs'));
xticks(1:Fs:Fs*6+1);
xticklabels({'-3', '-2', '-1', '0', '1', '2', '3'});
colormap(cmap);

hold on,
plot(-epoch_start_offset_253 - sortedRT, 1:length(sortedRT), 'linewidth',2,'color','k')
line([-epoch_start_offset_253, -epoch_start_offset_253], [1, length(sortedRT)],'linewidth',2,'color','w')
plot(-epoch_start_offset_253 - sortedRT+sortedRT_off, 1:length(sortedRT), 'linewidth',0.01,'color',[0.8,0.8,0.8])

if toSave
    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf, strcat('all_session_253','.fig'))
    print('all_session_253', '-djpeg')
end

% Plot 254
figure('Visible', visible);
imagesc(all_session_254);
set(gca,'Title',text('String','Time locked to driver response end'));
set(gca,'XLabel',text('String','Time offset (sec)'));
set(gca,'YLabel',text('String','Epochs'));
xticks(1:Fs:Fs*6+1);
xticklabels({'-4', '-3', '-2', '-1', '0', '1', '2'});
colormap(cmap);

hold on,
plot(-epoch_start_offset_254 - sortedRT_off, 1:length(sortedRT), 'linewidth',0.01,'color','k')
plot(-epoch_start_offset_254 - sortedRT_off+sortedRT, 1:length(sortedRT), 'linewidth',0.01,'color','w')
line([-epoch_start_offset_254, -epoch_start_offset_254], [1, length(sortedRT)], 'linewidth',2,'color',[0.8,0.8,0.8])

if toSave
    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf, strcat('all_session_254','.fig'))
    print('all_session_254', '-djpeg')
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
    for it = 1:floor((nTime-winLen+walkLen)/walkLen)
        time_range = (it-1)*walkLen+1:(it-1)*walkLen+winLen;
        stateProb(:,it) = mean(Gamma(time_range,:),1)';
    end
    end_minutes = (size(stateProb,2)-1)*walkLen/Fs/60;
    
    n_subplots = K + 2;
    figure('Visible', visible);
%     imagesc([0 end_minutes], [1 size(stateProb,1)], stateProb), colorbar;
%     colormap(gamma_axis, hot);
%     title('State Probability');
%     ylabel('State')
        
    [~, corresponding_vpath] = max(stateProb, [], 1);
    vpath_axis = subplot(n_subplots,1,1);
    imagesc([0 end_minutes], [1 1], corresponding_vpath); colorbar;
%     colorbar_ytick = 0:(K+1)/K:(K+1);
%     colorbar_ytick(1) = [];
%     colorbar('YTick', [1.33 2 2.66], 'YTickLabel', {'1', '2', '3'});
    colormap(vpath_axis, cmap)
    title('Viterbi path')
    yticks([]);

    rt_axis = subplot(n_subplots,1,2); plot(rt_latency/Fs/60, movmean(rt_speed, 7), 'k');
    title('Reaction Speed')
    ylabel('Speed (1/sec)')
    xlim([0, vpath_axis.XLim(2)*rt_axis.Position(3)/vpath_axis.Position(3)])
    
    for state = 1:K
        state_axis = subplot(n_subplots,1,state+2);
        plot(1:(end_minutes-1)/(size(stateProb,2)-1):end_minutes, movmean(stateProb(state,:),7), 'Color', cmap(state,:))
        r = state_r_list(state);
        title(strcat('State ', num2str(state)))
        ylabel('Probability')
        txt = ['r = ' num2str(r)];
        text(state_axis.Position(3), state_axis.Position(4), txt, 'FontSize', 10)
    end
    
    if toSave
        output_filename = strcat(filename,'_','path');
        set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf, strcat(output_filename,'.fig'))
        print(output_filename, '-djpeg')
    end
end

% %% Table view overall state probability
% figure;
% t = uitable('Data', results.overall_state_prob_list, 'ColumnName', {'State 1', 'State 2', 'State 3'}, 'RowName', results.filename_list);
% t.Position(3:4) = t.Extent(3:4);


%% Heatmap overall state probability
% figure;
% heatmap({'1','2','3'}, results.session_name_list, round(results.overall_state_prob_list, 3));
% title('Overall State Probability');
% xlabel('State');

%% Get top model specific state probability trial
% drowsy_top5_trial_r_list = zeros(n_of_files, 1);
% drowsy_top30_trial_r_list = zeros(n_of_files, 1);
% alert_top5_trial_r_list = zeros(n_of_files, 1);
% alert_top30_trial_r_list = zeros(n_of_files, 1);
% top = 1:5;
% wide_top = 1:100;
% 
% for idx = 1:n_of_files
%     trial_Gamma_mean = results{idx, 'trial_Gamma_mean_list'}{1};
%     state_type_index = results{idx, 'permutation_list'};
%     rs = results{idx, 'rt_speed_clean_list'}{1};
%     
%     drowsy_index = state_type_index(1);
%     alert_index = state_type_index(3);
%     [sorted_drowsy_mean, drowsy_sort_index] = sort(trial_Gamma_mean(:, drowsy_index), 1);
%     [sorted_alert_mean, alert_sort_index] = sort(trial_Gamma_mean(:, alert_index), 1);
%     
%     drowsy_top5_trial_r_list(idx) = corr(sorted_drowsy_mean(top, :), rs(drowsy_sort_index(top))');
%     drowsy_top30_trial_r_list(idx) = corr(sorted_drowsy_mean(wide_top, :), rs(drowsy_sort_index(wide_top))');
%     alert_top5_trial_r_list(idx) = corr(sorted_alert_mean(top, :), rs(alert_sort_index(top))');
%     alert_top30_trial_r_list(idx) = corr(sorted_alert_mean(wide_top, :), rs(alert_sort_index(wide_top))');      
% end

%% Plotting transition probability matrix
visible = 'off';
toSave = 1;
for idx = 1:n_of_files
    hmm = results{idx, 'hmm_list'}{1};
    trans_P = hmm.P;
    state_permutation = results{idx, 'permutation_list'};
    filename = results{idx, 'filename_list'}{1};
    filename = split(filename, '.');
    filename = filename{1}
    
    trans_P = trans_P(state_permutation, state_permutation);
    trans_P = round(trans_P, 3);
    figure('Visible', visible)
    imagesc(trans_P), colorbar
    set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Drowsy', 'Middle', 'Alert'});
    set(gca, 'YTick', [1 2 3], 'YTickLabel', {'Drowsy', 'Middle', 'Alert'});
        
    textStrings = num2str(trans_P(:));       % Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    [x, y] = meshgrid(1:3);  % Create x and y coordinates for the strings
    hStrings = text(x(:), y(:), textStrings(:), 'HorizontalAlignment', 'center');
    textColors = repmat(trans_P(:) < 0.5, 1, 3);  % Choose white or black 
    set(hStrings, {'Color'}, num2cell(textColors, 2), {'FontSize'}, num2cell(repmat(13, 9, 1)));
      
    if toSave
        output_filename = strcat(filename, '_', 'transProb');
        set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf, strcat(output_filename,'.fig'))
        print(output_filename, '-djpeg')
    end
end

%% Plot dendrograms of covariance matrix
isVisible = "off";
toSave = 1;

set(0, 'ShowHiddenHandles', 'on');
for idx = 1:n_of_files
    hmm = results{idx, 'hmm_list'}{1};
    chanlocs = results{idx, 'chanlocs_list'}{1};
    permutation_list = results{idx, 'permutation_list'};
    filename = results{idx, 'filename_list'}{1};
    filename = split(filename, '.');
    filename = filename{1}
        
    for state = permutation_list
        data = getFuncConn(hmm, state);
        data = reshape(normalize(reshape(data, 1, [])), size(data));
        
        cg = clustergram(data, 'Cluster', 2, 'ColumnLabels', chanlocs);
        allhandles = get(0, 'Children');
        cgfigidx = strcmp('Clustergram', get(allhandles, 'Tag'));
        close(allhandles(cgfigidx));
        
        [~, perm] = ismember(get(cg, 'ColumnLabels'), chanlocs);        
        data = data(perm, perm);
        chanlocs = chanlocs(perm);
        cg = clustergram(data, 'Cluster', 2, 'Colormap', redblue, 'Dendrogram', 'default', ...
            'ColumnLabels', chanlocs, 'ColumnLabelsRotate', 60, 'RowLabels', chanlocs);
        allhandles = get(0, 'Children');
        cgfigidx = strcmp('Clustergram', get(allhandles, 'Tag'));
        cgfighandle = allhandles(cgfigidx);
        cgfighandle = cgfighandle(end);
        set(cgfighandle, 'Position', [100 100 500 500]);
        addTitle(cg, strcat('Functional Connectivity of state', {' '}, state_description{state}));
        
        % Create colorbar
        cbButton = findall(cgfighandle, 'tag', 'HMInsertColorbar');
        ccb = get(cbButton, 'ClickedCallback');
        set(cbButton, 'State', 'on');
        ccb{1}(cbButton, [], ccb{2});        
        cgfighandle.Visible = isVisible;
        
        if toSave
            output_filename = strcat(filename, '_', 'FunConn', '_', state_description{state});
            set(cgfighandle, 'PaperPositionMode', 'auto');
            saveas(cgfighandle, strcat(output_filename,'.fig'))
            print(output_filename, '-djpeg')
        end
        
    end
end

set(0, 'ShowHiddenHandles', 'off');

%% Plot dendrograms of covariance matrix
isVisible = "off";
toSave = 1;

for idx = 1:n_of_files
    hmm = results{idx, 'hmm_list'}{1};
    chanlocs = results{idx, 'chanlocs_list'}{1};
    permutation_list = results{idx, 'permutation_list'};
    filename = results{idx, 'filename_list'}{1};
    filename = split(filename, '.');
    filename = filename{1}
        
    for state = permutation_list
        data = getFuncConn(hmm, state);
%         data = reshape(quantilenorm(reshape(data, [], 1)), size(data,1), []);
        figure('Visible', isVisible);
        imagesc(data);
        colormap(redblue); colorbar;
        set(gca, 'YDir', 'normal');
        set(gca, 'XTick', 1:length(chanlocs), 'XTickLabel', chanlocs, 'XTickLabelRotation', 60);
        set(gca, 'YTick', 1:length(chanlocs), 'YTickLabel', chanlocs);
        title(strcat('Functional Connectivity of state', {' '}, state_description{state}));
        
        if toSave
            output_filename = strcat(filename, '_', 'FunConnRaw', '_', state_description{state});
            set(gcf, 'PaperPositionMode', 'auto');
%             saveas(cgfighandle, strcat(output_filename,'.fig'))
            print(output_filename, '-djpeg')
        end
        
    end
end

set(0, 'ShowHiddenHandles', 'off');