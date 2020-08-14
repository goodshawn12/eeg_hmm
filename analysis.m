cd '/home/ting/Documents/eeg_hmm';
addpath('/home/ting/Documents/eeg_hmm');
addpath(genpath('/home/ting/Documents/eeg_hmm/Utility'))
addpath(genpath('/home/ting/Documents/eeg_hmm/HMM-MAR'));
addpath('/data/projects/Shawn/2016 JNE/dataset/');
addpath('/data/projects/Shawn/2019_HMM/data/');

%% Specify number of inferred states
K = 3
method = 'AMICA'
if strcmp(method, 'AMICA')
    isAMICA = 1;
else
    isAMICA = 0;
end
if isAMICA
    resultsDir = sprintf('/home/ting/Documents/eeg_hmm/AMICA_results/results_K%d/', K);
    graphDir = sprintf('/home/ting/Documents/eeg_hmm/AMICA_results/graphs_K%d/', K);
else
    resultsDir = sprintf('/home/ting/Documents/eeg_hmm/HMM_results/results_K%d/', K);
    graphDir = sprintf('/home/ting/Documents/eeg_hmm/HMM_results/graphs_K%d/', K);
end

if ~exist(graphDir, 'dir')
    mkdir(graphDir)
end

%% Manual setup files to process
if isAMICA
    data_filelist = dir(resultsDir);
else
    data_filelist = dir(strcat(resultsDir, method, '*_2.mat'));
end

filename_list = cell(length(data_filelist), 1);
for j = 1:length(data_filelist)
    file = data_filelist(j);
    if contains(file.name, 'ASR20') || strcmp(file.name,'.') || strcmp(file.name,'..') || (isAMICA && ~file.isdir)
        continue;
    end
    filename = split(file.name, '.');
    filename_list{j} = filename{1};
end
filename_list = filename_list(~cellfun('isempty', filename_list));
filename_list = unique(filename_list, 'stable');

n_of_files = length(filename_list)

% Prepare raw data
training_data_size_list = zeros(n_of_files,2);
Gamma_list = cell(n_of_files,1);
vpath_list = cell(n_of_files,1);
if ~isAMICA
    hmm_list = cell(n_of_files,1);
    vars_to_load = {'hmm', 'training_data_size', 'Gamma', 'vpath'};
end

if isAMICA
    run('/data/common/matlab/eeglab/eeglab');
end

for idx = 1:n_of_files
    if strcmp(method, 'AMICA')        
        fileDir = strcat(resultsDir, filename_list{idx});
        AMICA_out = loadmodout15(fileDir);
        training_data_size_list(idx,:) = [size(AMICA_out.v, 2), AMICA_out.data_dim];
        [vpath, Gamma] = processVpathAMICA(AMICA_out.v);
        Gamma_list{idx} = Gamma;
        vpath_list{idx} = vpath;
    else
        HMM_out = load(strcat(resultsDir, filename_list{idx}, '.mat'), vars_to_load{:});
        training_data_size_list(idx,:) = HMM_out.training_data_size;
        hmm_list{idx} = HMM_out.hmm;
        Gamma_list{idx} = HMM_out.Gamma;
        vpath_list{idx} = HMM_out.vpath;
    end    
end

results = table(filename_list, training_data_size_list, Gamma_list, vpath_list);

if ~isAMICA
    results = deleteVars(results, {'hmm_list'});
    results = addvars(results, hmm_list);
end

fprintf('Raw results data loaded.\n')

%% Load EEG related metadata
events_list = cell(n_of_files, 1);
chanlocs_list = cell(n_of_files, 1);
Fs_list = zeros(n_of_files,1);

for idx = 1:n_of_files
    if ~isAMICA
        filename = split(results.filename_list{idx}, '_');
        filename = join(filename(2:length(filename)-1), '_');
        filename = strcat(filename{1}, '.set');
    else
        filename = strcat(results.filename_list{idx}, '.set');
    end
    
    load('-mat', filename);
    events_list{idx} = EEG.event;
    chanlocs_list{idx} = extractfield(EEG.chanlocs, 'labels');
    Fs_list(idx) = EEG.srate;
end

metaInfo = table(events_list, chanlocs_list, Fs_list);

fprintf('EEG meta data loaded.\n')

%% Generate cleaned rt, rt_off
rt_list = cell(n_of_files, 1);
rt_off_list = cell(n_of_files, 1);
rt_latency_list = cell(n_of_files, 1);
rs_list = cell(n_of_files, 1);
event_sparsity_list = cell(n_of_files, 1);

for idx = 1:n_of_files
    Fs = metaInfo{idx, 'Fs_list'};
    events = metaInfo{idx, 'events_list'}{1};
    
    [rt, rt_off, rt_latency, rs, event_sparsity] = rtFromEvents(events, Fs, 1);
    
    rt_list{idx} = rt;
    rt_off_list{idx} = rt_off;
    rt_latency_list{idx} = rt_latency;
    rs_list{idx} = rs;
    event_sparsity_list{idx} = event_sparsity;
    
end

metaInfo = deleteVars(metaInfo, {'rt_list', 'rt_off_list', 'rt_latency_list', 'rs_list', 'event_sparsity_list'});
metaInfo = addvars(metaInfo, rt_list, rt_off_list, rt_latency_list, rs_list, event_sparsity_list);

fprintf('rt generated.\n');

%% Smoothed state correlation with RS
state_r_list = zeros(n_of_files, K);
win_Gamma_median_list = cell(n_of_files, 1);

win_len_sec = 5;
smoothing_range_sec = 90;

for idx = 1:n_of_files
    Fs = metaInfo{idx, 'Fs_list'};
    Gamma = results{idx, 'Gamma_list'}{1};
    rs = metaInfo{idx, 'rs_list'}{1};
    rt_latency = metaInfo{idx, 'rt_latency_list'}{1};

    [correlation, win_Gamma_median] = getCorr(Gamma, K, Fs, rs, rt_latency, win_len_sec, smoothing_range_sec);
    state_r_list(idx, :) = round(correlation, 3);
    win_Gamma_median_list{idx} = win_Gamma_median;
end

% The sort is in ascending order so the first index in the permutation_list
% marks the lowest correlation (drowsy model) 
[~, permutation_list] = sort(state_r_list, 2);
inverse_permutation_list = zeros(size(permutation_list));
for row=1:size(permutation_list,1)
    inverse_permutation_list(row, permutation_list(row,:)) = 1:size(permutation_list,2);
end

results = deleteVars(results, {'state_r_list', 'permutation_list', 'inverse_permutation_list', 'win_Gamma_median_list'});
results = addvars(results, state_r_list, permutation_list, inverse_permutation_list, win_Gamma_median_list);

fprintf('correlation generated.\n');

%% Examine fractional occupancy
fo_list = zeros(n_of_files, K);

for idx = 1:n_of_files
    Gamma = results{idx, 'Gamma_list'}{1};
    Gamma_valid = Gamma(any(Gamma,2),:);
    permutation = results{idx, 'permutation_list'};
    
    fo = mean(Gamma_valid, 1);
    fo_list(idx,:) = fo(permutation);
end
results = deleteVars(results, {'fo_list'});
results = addvars(results, fo_list);

fprintf('fractional occupancy generated.\n');

%% Compute transitional probability matrix
transProb_list = cell(n_of_files, 1);
for idx = 1:n_of_files
    if isAMICA
        vpath = results{idx, 'vpath_list'}{1};
        transProb = getTransProbMatrix(vpath, K);
    else
        hmm = results{idx, 'hmm_list'}{1};
        transProb = hmm.P;
    end
    transProb_list{idx} = transProb;  
end

results = deleteVars(results, {'transProb_list'});
results = addvars(results, transProb_list);

fprintf('transition probability matrix generated.\n')

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

%% ######## Graphing Sections ########
%########################################

% Preparation
[cmap, state_description] = getCmap(K);

cmap_list = repmat(cmap, 1, 1, n_of_files);
for i = 1:n_of_files
    cmap_list(:,:,i) = cmap_list(results.inverse_permutation_list(i,:),:,i);
end
toSave = 0;
visible = 'off';

fprintf('graphing prepared.\n');

%% Plot smoothed state time path
visible = 'off';
toSave = 1;
for idx = 1:n_of_files
    Fs = metaInfo{idx, 'Fs_list'};
    Gamma_median = results{idx, 'win_Gamma_median_list'}{1};
%     vpath = results{idx, 'vpath'}{1};
    rt = metaInfo{idx, 'rt_list'}{1};
    rt_latency = metaInfo{idx, 'rt_latency_list'}{1};
    rs = metaInfo{idx, 'rs_list'}{1};
    state_r_list = results{idx, 'state_r_list'};
    filename = results{idx, 'filename_list'}{1}
    
    smooth_rs = smoothRS(rs, rt_latency, 90, Fs);
    cmap = cmap_list(:,:,idx);
       
    xData = rt_latency;
    unit_time_minutes = 10;
    rt_latency_minutes = rt_latency/Fs/60;
    max_minutes = ceil(rt_latency_minutes(end)/unit_time_minutes)*unit_time_minutes;
    x_ticks_label_data = 0:unit_time_minutes:max_minutes;
    x_ticks_label = cell(1, length(x_ticks_label_data));
    x_ticks = x_ticks_label_data*Fs*60;
    for j = 1:length(x_ticks_label_data)
        x_ticks_label{j} = num2str(x_ticks_label_data(j));
    end
    
    n_subplots = K+1;
    fig1 = figure('Visible', visible);
    p = uipanel('Parent', fig1);   

    rs_axis = subplot(n_subplots,1,1, 'Parent', p);
    plot(xData, rs, '.', 'LineWidth', 0.2, 'Color', [0.8,0.8,0.8]); hold on;
    plot(xData, smooth_rs, '-', 'LineWidth', 2, 'Color', 'k'); 
    title('Reaction Speed by Events')
    ylabel('Speed (1/sec)')
    xticks(x_ticks);
    xticklabels(x_ticks_label);
    
    [~, permutation_list] = sort(state_r_list, 'descend');
    for j = 1:K
        state = permutation_list(j);
        state_axis = subplot(n_subplots,1,j+1);
        yData = Gamma_median(:,state);
                
        plot(xData, yData, 'Color', cmap(state,:))
        r = state_r_list(state);
        title(state_description(K-j+1));
        ylabel('Probability');
        xticks(x_ticks);
        xticklabels(x_ticks_label);
        txt = ['r = ' num2str(r)];
        y_center = mean(state_axis.YLim);
        if mean(yData(rt_latency < unit_time_minutes*60*Fs)) > mean(y_center)
            text_y_position = y_center - range(state_axis.YLim)/4;
        else
            text_y_position = y_center + range(state_axis.YLim)/4;
        end    
        text(2*60*Fs, text_y_position, txt, 'FontSize', 10, 'Parent', state_axis)
    end
    xlabel('Time (minutes)')
    
    if toSave
        output_filename = strcat(graphDir, filename,'_','path');
        set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(output_filename,'.fig'))
        print(output_filename, '-djpeg')
    end
end

fprintf('vpath plot done.\n')

%% Epoching by timelock at events
epoch_251_list = cell(n_of_files,1);
epoch_253_list = cell(n_of_files,1);
epoch_254_list = cell(n_of_files,1);
epoch_rt_remove_list = cell(n_of_files,1);

for idx = 1:n_of_files
    vpath = results{idx, 'vpath_list'}{1}; 
    training_data_size = results{idx, 'training_data_size_list'};
    Fs = metaInfo{idx, 'Fs_list'};
    rt = metaInfo{idx, 'rt_list'}{1};
    rt_off = metaInfo{idx, 'rt_off_list'}{1};
    rt_latency = metaInfo{idx, 'rt_latency_list'}{1};
    
    vpath = [zeros(training_data_size(1)-size(vpath,1),1); vpath];
    
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
    epoch_rt_remove_list{idx} = any(rt_latency + epoch_start_offset_251 < 1, 1)...
        | any(rt_latency + rt_off + epoch_end_offset_254 > size(vpath,1), 1);    
    
    % Epoching to timelock window
    epoch_251_list{idx} = epochByEvent(vpath, rt_latency, epoch_start_offset_251, epoch_end_offset_251);        
    epoch_253_list{idx} = epochByEvent(vpath, rt_latency+rt, epoch_start_offset_253, epoch_end_offset_253);   
    epoch_254_list{idx} = epochByEvent(vpath, rt_latency+rt_off, epoch_start_offset_254, epoch_end_offset_254);   
end

results = deleteVars(results, {'epoch_251_list', 'epoch_253_list', 'epoch_254_list', 'epoch_rt_remove_list'});
results = addvars(results, epoch_251_list, epoch_253_list, epoch_254_list, epoch_rt_remove_list);

%% Colormap epoched data
visible = 'off';
toSave = 1;

vertical_smoothing_window_len = 5;
isSorted = 1;
win_len_sec = 6;

for idx = 1:n_of_files
    cmap = cmap_list(:,:,idx);
    
    Fs = metaInfo{idx, 'Fs_list'};
    rt = metaInfo{idx, 'rt_list'}{1};
    rt_off = metaInfo{idx, 'rt_off_list'}{1};
    rt_latency = metaInfo{idx, 'rt_latency_list'}{1};
    removeIndex = results{idx, 'epoch_rt_remove_list'}{1};
    filename = results{idx, 'filename_list'}{1}
    
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
    rt(removeIndex) = [];
    rt_off(removeIndex) = [];
    rt_latency(removeIndex) = [];
    
    if isSorted
        [rt, sortIdx] = sort(rt, 'descend');
        rt_off = rt_off(:, sortIdx);
        rt_latency = rt_latency(:, sortIdx);
    end
    
    % Plot 251/252
    epoch_251 = results{idx, 'epoch_251_list'}{1};
    epoch_251 = epoch_251(sortIdx,:);
    epoch_251 = movingModeSmoothing(epoch_251, vertical_smoothing_window_len, 1, 1:K);     
    figure('Visible', visible);
    plotTimelock(gca, '251', epoch_251, epoch_start_offset_251, Fs, win_len_sec, rt, rt_off, cmap, 1);

    if toSave
        set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(filename,'_','251','.fig'))
        print(strcat(graphDir, filename,'_','251'), '-djpeg')  
    end
    
    % Plot 253
    epoch_251 = results{idx, 'epoch_253_list'}{1};
    epoch_251 = epoch_251(sortIdx,:);
    epoch_251 = movingModeSmoothing(epoch_251, vertical_smoothing_window_len, 1, 1:K); 

    figure('Visible', visible);
    plotTimelock(gca, '253', epoch_251, epoch_start_offset_253, Fs, win_len_sec, rt, rt_off, cmap, 1);
    
    if toSave
        set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(filename,'_','253','.fig'))
        print(strcat(graphDir, filename,'_','253'), '-djpeg')
    end
       
    % Plot 254
    epoch_251 = results{idx, 'epoch_254_list'}{1};
    epoch_251 = epoch_251(sortIdx,:);
    epoch_251 = movingModeSmoothing(epoch_251, vertical_smoothing_window_len, 1, 1:K); 

    figure('Visible', visible);
    plotTimelock(gca, '254', epoch_251, epoch_start_offset_254, Fs, win_len_sec, rt, rt_off, cmap, 1);
    
    if toSave
        set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(filename,'_','254','.fig'))
        print(strcat(graphDir, filename,'_','254'), '-djpeg')
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
    rt = metaInfo{idx, 'rt_list'}{1};
    rt_off = metaInfo{idx, 'rt_off_list'}{1};
    rt_latency = metaInfo{idx, 'rt_latency_list'}{1};
    
    permutation = results{idx, 'permutation_list'};    
    ordered_vpath = vpath;
    for state = 1:K
        ordered_vpath(vpath == permutation(state)) = state;
    end
    
    % Prevent epoch exceeding vpath dimension
    rt(removeIndex) = [];
    rt_off(removeIndex) = [];
    rt_latency(removeIndex) = [];
    
    all_session_RT = cat(2, all_session_RT, rt);
    all_session_RT_latency = cat(2, all_session_RT_latency, rt_latency);
    all_session_RT_off = cat(2, all_session_RT_off, rt_off);
    
    epoch_251 = results{idx, 'epoch_251_list'}{1};
    epoch_253 = results{idx, 'epoch_253_list'}{1};
    epoch_254 = results{idx, 'epoch_253_list'}{1};
    
    for state = 1:K
        epoch_251(epoch_251 == permutation(state)) = state;
        epoch_253(epoch_253 == permutation(state)) = state;
        epoch_254(epoch_254 == permutation(state)) = state;
    end
    
    all_session_251 = cat(1, all_session_251, epoch_251);  
    all_session_253 = cat(1, all_session_253, epoch_251);     
    all_session_254 = cat(1, all_session_254, epoch_251);   
end

[rt, sortIdx] = sort(all_session_RT, 'descend');
rt_off = all_session_RT_off(:, sortIdx);
all_session_251 = all_session_251(sortIdx, :);
all_session_253 = all_session_253(sortIdx, :);
all_session_254 = all_session_254(sortIdx, :);

%% Plot all session timelock
visible = 'on';
toSave = 1;

vertical_smoothing_window_len = 5;
cmap = getCmap(K);
filename = 'all_session_timelock';

% Plot 251
all_session_251 = movingModeSmoothing(all_session_251, vertical_smoothing_window_len, 1, 0:K);
figure('Visible', visible);
plotTimelock(gca, '251', all_session_251, epoch_start_offset_251, Fs, win_len_sec, sortedRT, sortedRT_off, cmap, 0);

if toSave
    set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(filename,'_','251','.fig'))
    print(strcat(graphDir, filename,'_251'), '-djpeg')  
end

% Plot 253
all_session_253 = movingModeSmoothing(all_session_253, vertical_smoothing_window_len, 1, 0:K);
figure('Visible', visible);
plotTimelock(gca, '253', all_session_253, epoch_start_offset_253, Fs, win_len_sec, sortedRT, sortedRT_off, cmap, 0);

if toSave
    set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(filename,'_','253','.fig'))
    print(strcat(graphDir, filename,'_','253'), '-djpeg')
end

% Plot 254
all_session_254 = movingModeSmoothing(all_session_254, vertical_smoothing_window_len, 1, 0:K); 
figure('Visible', visible);
plotTimelock(gca, '254', all_session_254, epoch_start_offset_254, Fs, win_len_sec, sortedRT, sortedRT_off, cmap, 0);

if toSave
    set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(filename,'_','254','.fig'))
    print(strcat(graphDir, filename,'_','254'), '-djpeg')
end

