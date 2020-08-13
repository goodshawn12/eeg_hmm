% cd '/home/ting/Documents/eeg_hmm';
% addpath('/home/ting/Documents/eeg_hmm');
% addpath(genpath('/home/ting/Documents/eeg_hmm/Utility'))
% addpath('/data/projects/Shawn/2019_HMM/data/');
% run('/data/common/matlab/eeglab/eeglab');

%% Specify number of inferred states
K = 5;
resultsDir = sprintf('/home/ting/Documents/eeg_hmm/AMICA_results/results_K%d/', K);
graphDir = sprintf('/home/ting/Documents/eeg_hmm/AMICA_results/graphs_K%d/', K)
if ~exist(graphDir, 'dir')
    mkdir(graphDir)
end

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

n_of_files = length(filename_list)

%% Prepare raw data
n_epochs = 1;

v_list = cell(n_of_files,1);

for idx = 1:n_of_files
    fileDir = strcat(resultsDir, filename_list{idx});
    modout = loadmodout15(fileDir);
    
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
event_sparsity_list = cell(n_of_files, 1);

for idx = 1:n_of_files
    Fs = results{idx, 'Fs_list'}{1};
    events = results{idx, 'events_list'}{1};  
    rt = zeros(1,length(events));
    rt_off = zeros(1,length(events));
    rt_latency = zeros(1,length(events));
    rt_speed = zeros(1, length(events));
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
            rt_speed(i) = 1./(rt(i)/Fs); %rt_speed is in (1/seconds)
            
            if ~isnan(prev_offset)
                event_sparsity(i) = event.latency - prev_offset;
            end
            prev_offset = events(i+2).latency;
        end
    end
    nonempty_index = any(rt_latency, 1);
    rt_list{idx} = rt(nonempty_index);
    rt_off_list{idx} = rt_off(nonempty_index);
    rt_latency_list{idx} = rt_latency(nonempty_index);
    rt_speed_list{idx} = rt_speed(nonempty_index);
    event_sparsity = event_sparsity(nonempty_index);
    event_sparsity_list{idx} = event_sparsity(2:end); % The first 251 event has no prior LDT offset 
end
results = deleteVars(results, {'rt_list', 'rt_off_list', 'rt_latency_list', 'rt_speed_list', 'event_sparsity_list'});
results = addvars(results, rt_list, rt_off_list, rt_latency_list, rt_speed_list, event_sparsity_list);

% RT and 1/RT Outlier Removal
rt_clean_list = cell(n_of_files, 1);
rt_off_clean_list = cell(n_of_files, 1);
rt_latency_clean_list = cell(n_of_files, 1);
rt_speed_clean_list = cell(n_of_files, 1);
rt_removal_percentage_list = zeros(n_of_files, 1);

for idx = 1:n_of_files
    Fs = results{idx, 'Fs_list'}{1};
    rt = results{idx, 'rt_list'}{1};
    rt_off = results{idx, 'rt_off_list'}{1};
    rt_latency = results{idx, 'rt_latency_list'}{1}; 
    rt_speed = results{idx, 'rt_speed_list'}{1};
    
    lower_remove_index = rt < 0.1 * Fs;
    upper_remove_index = rt > 10 * Fs;
    rt_off_remove_index = (rt_off - rt) < 0.1 * Fs;
    remove_index = lower_remove_index | upper_remove_index | rt_off_remove_index;
    rt_removal_percentage_list(idx) = sum(remove_index) / length(remove_index);
    
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

results = deleteVars(results, {'rt_clean_list', 'rt_off_clean_list', 'rt_latency_clean_list', 'rt_speed_clean_list', 'rt_removal_percentage_list'});
results = addvars(results, rt_clean_list, rt_off_clean_list, rt_latency_clean_list, rt_speed_clean_list, rt_removal_percentage_list);

% Get RS statistics
rs_mean_list = zeros(n_of_files, 1);
rs_st_list = zeros(n_of_files, 1);
for idx = 1:n_of_files
    rs = results{idx, 'rt_speed_clean_list'}{1};
    rs_mean_list(idx) = mean(rs);
    rs_st_list(idx) = std(rs);
end

fprintf('rt processing done.\n')

%% Global state correlation with ERP
state_r_list = zeros(n_of_files, K);
trial_Gamma_mean_list = cell(n_of_files, 1);
win_len_sec = 5;

for idx = 1:n_of_files
    Fs = results{idx, 'Fs_list'}{1};
    Gamma = results{idx, 'Gamma_list'}{1};
    rt_speed = results{idx, 'rt_speed_clean_list'}{1};
    rt_latency = results{idx, 'rt_latency_clean_list'}{1};
    
    win_offset = -win_len_sec * Fs;
    win_Gamma_mean = zeros(length(rt_latency), K);    
    for row = 1:length(win_Gamma_mean)
        latency = rt_latency(row);
        Gamma_window = Gamma(latency+win_offset:latency,:);
        Gamma_window_valid = Gamma_window(any(Gamma_window,2),:);
        win_Gamma_mean(row,:) = mean(Gamma_window_valid ,1);
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

%% Global smoothed state correlation with ERP
state_r_list = zeros(n_of_files, K);
trial_Gamma_mean_list = cell(n_of_files, 1);

smoothed_state_r_list = zeros(n_of_files, K);
smoothed_Gamma_median_list = cell(n_of_files, 1);

win_len_sec = 5;
smoothing_range_sec = 90;

for idx = 1:n_of_files
    Fs = results{idx, 'Fs_list'}{1};
    Gamma = results{idx, 'Gamma_list'}{1};
    rt_speed = results{idx, 'rt_speed_clean_list'}{1};
    rt_latency = results{idx, 'rt_latency_clean_list'}{1};
    
    win_offset = -win_len_sec * Fs;
    win_Gamma_mean = zeros(length(rt_latency), K);    
    for row = 1:length(win_Gamma_mean)
        latency = rt_latency(row);
        Gamma_window = Gamma(latency+win_offset:latency,:);
        Gamma_window_valid = Gamma_window(any(Gamma_window,2),:);
        win_Gamma_mean(row,:) = mean(Gamma_window_valid ,1);
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
    smoothed_Gamma_median = zeros(length(rt_latency), K);
    for row = 1:length(smoothed_Gamma_median)
        trials_selection = (row - facing_trial_num_list(row)):row;
        smoothed_Gamma_median(row, :) = median(win_Gamma_mean(trials_selection,:), 1);
    end    
  
    % Get correlation
    trial_Gamma_mean_list{idx} = win_Gamma_mean;
    for state = 1:K
        state_r_list(idx, state) = round(corr(win_Gamma_mean(:,state), rt_speed'), 3);
    end
    
    smoothed_Gamma_median_list{idx} = smoothed_Gamma_median;
    for state = 1:K
        smoothed_state_r_list(idx, state) = round(corr(smoothed_Gamma_median(:, state), rt_speed'), 3);
    end
end

% The sort is in ascending order so the first index in the permutation_list
% marks the lowest correlation (drowsy model) 
[~, permutation_list] = sort(state_r_list, 2);
inverse_permutation_list = zeros(size(permutation_list));
for row=1:size(permutation_list,1)
    inverse_permutation_list(row, permutation_list(row,:)) = 1:size(permutation_list,2);
end

[~, smoothed_permutation_list] = sort(smoothed_state_r_list, 2);
smoothed_inverse_permutation_list = zeros(size(smoothed_permutation_list));
for row=1:size(smoothed_permutation_list,1)
    smoothed_inverse_permutation_list(row, smoothed_permutation_list(row,:)) = 1:size(smoothed_permutation_list,2);
end

results = deleteVars(results, {'trial_Gamma_mean_list', 'state_r_list', 'permutation_list', 'inverse_permutation_list', ...
    'smoothed_Gamma_median_list', 'smoothed_state_r_list', 'smoothed_permutation_list', 'smoothed_inverse_permutation_list'});
results = addvars(results, trial_Gamma_mean_list, state_r_list, permutation_list, inverse_permutation_list, ...
    smoothed_Gamma_median_list, smoothed_state_r_list, smoothed_permutation_list, smoothed_inverse_permutation_list);

sort(results.smoothed_state_r_list,2)

%% Examine fractional occupancy
fo_list = zeros(n_of_files, K);

for idx = 1:n_of_files
    Gamma = results{idx, 'Gamma_list'}{1};
    Gamma_valid = Gamma(any(Gamma,2),:);
    permutation = results{idx, 'smoothed_permutation_list'};
    
    fo = mean(Gamma_valid, 1);
    fo_list(idx,:) = fo(permutation);
end
results = deleteVars(results, {'fo_list'});
results = addvars(results, fo_list);
fprintf('fractional occupancy generated');

fo_list
keyboard;

%% Examine overall state-state transition occurance frequency
transFreq_list = cell(n_of_files, 1);
for idx = 1:n_of_files
    vpath = results{idx, 'vpath_list'}{1};
    transFreq = getTransFreqMatrix(vpath, K);
    transFreq_list{idx} = transProb;  
end

results = deleteVars(results, {'transFreq_list'});
results = addvars(results, transFreq_list);

%% Compute transitional probability matrix
transProb_list = cell(n_of_files, 1);
for idx = 1:n_of_files
    vpath = results{idx, 'vpath_list'}{1};
    transProb = getTransProbMatrix(vpath, K);
    transProb_list{idx} = transProb;  
end

results = deleteVars(results, {'transProb_list'});
results = addvars(results, transProb_list);

%% ######## Graphing Sections ########
%########################################

graphDir = sprintf('/home/ting/Documents/eeg_hmm/AMICA_results/graphs_K%d', K);
cd(graphDir);

% Preparation
[cmap, state_description] = getCmap(K);

cmap_list = repmat(cmap, 1, 1, n_of_files);
for i = 1:n_of_files
    cmap_list(:,:,i) = cmap_list(results.inverse_permutation_list(i,:),:,i);
end
toSave = 0;
visible = 'off';

%% Plot overall transitional Probability matrix
visible = 'on';
toSave = 1;

overall_transProb = zeros(K);
var_matrix = zeros(K);
for idx = 1:n_of_files
    transProb = results{idx, 'transProb_list'}{1};
    state_permutation = results{idx, 'smoothed_permutation_list'};
    
    transProb = transProb(state_permutation, state_permutation);
    overall_transProb = overall_transProb + transProb;
    var_matrix = var_matrix + transProb.^2;
end
overall_transProb = overall_transProb ./ n_of_files;
var_matrix = var_matrix ./ n_of_files - overall_transProb .^ 2;

figure('Visible', visible)
imagesc(overall_transProb), colorbar

TickLabel = cell(1, K);
TickLabel{1} = 'Drowsy'; TickLabel{K} = 'Alert';
TickLabel{2:K-1} = 'Middle';
set(gca, 'XTick', 1:K, 'XTickLabel', TickLabel, 'YTick', 1:K, 'YTickLabel', TickLabel);

overall_transProb = round(overall_transProb, 3);
textStrings = num2str(overall_transProb(:));       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:K);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:), 'HorizontalAlignment', 'center');
textColors = repmat(overall_transProb(:) < 0.5, 1, K);  % Choose white or black 
set(hStrings, {'Color'}, num2cell(textColors, 2), {'FontSize'}, num2cell(repmat(13, 9, 1)));

if toSave
    output_filename = strcat('mean_transProb', '_K', num2str(K));
    set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(output_filename,'.fig'))
    print(output_filename, '-djpeg')
end

%% Plotting transition probability matrix
visible = 'on';
toSave = 0;
for idx = 9
    transProb = results{idx, 'transProb_list'}{1};
    state_permutation = results{idx, 'smoothed_permutation_list'};
    filename = results{idx, 'filename_list'}{1};
    
    transProb = transProb(state_permutation, state_permutation);
    figure('Visible', visible)
    imagesc(transProb), colorbar
    
    TickLabel = cell(1, K);
    TickLabel{1} = 'Drowsy'; TickLabel{K} = 'Alert';
    TickLabel{2:K-1} = 'Middle';
    set(gca, 'XTick', 1:K, 'XTickLabel', TickLabel, 'YTick', 1:K, 'YTickLabel', TickLabel);
    
    transProb = round(transProb, 3);
    textStrings = num2str(transProb(:));       % Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    [x, y] = meshgrid(1:K);  % Create x and y coordinates for the strings
    hStrings = text(x(:), y(:), textStrings(:), 'HorizontalAlignment', 'center');
    textColors = repmat(transProb(:) < 0.5, 1, K);  % Choose white or black 
    set(hStrings, {'Color'}, num2cell(textColors, 2), {'FontSize'}, num2cell(repmat(13, 9, 1)));
      
    if toSave
        output_filename = strcat(filename, '_', 'transProb');
        set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(output_filename,'.fig'))
        print(output_filename, '-djpeg')
    end
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

%% Check for rejection density
% for idx = 4:5
%     vpath = results{idx, 'vpath_list'}{1};
%     vpath(any(vpath~=0,1)) = 1;
%     imagesc(vpath);
% end

%% Colormap by epochs timelock at events
visible = 'on';
toSave = 0;

smooth = 0;
vertical_smoothing_window_len = 1;
isSorted = 1;
win_len_sec = 6;


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
   
    [sortedRT, sortIdx] = sort(rt, 'descend');
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
    plotTimelock(gca, '251', state_by_epoch, epoch_start_offset_251, Fs, win_len_sec, sortedRT, sortedRT_off, cmap);

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
    plotTimelock(gca, '253', state_by_epoch, epoch_start_offset_253, Fs, win_len_sec, sortedRT, sortedRT_off, cmap);
    
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
    plotTimelock(gca, '254', state_by_epoch, epoch_start_offset_254, Fs, win_len_sec, sortedRT, sortedRT_off, cmap);
    
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

%% Prepare to plot all session timelock
all_session_251 = [];
all_session_253 = [];
all_session_254 = [];
all_session_RT = [];
all_session_RT_off = [];
all_session_RT_latency = [];

for idx = 1:n_of_files
    vpath = results{idx, 'smoothed_vpath_list'}{1};
    num_mods = results{idx, 'K_list'}{1};
    Fs = results{idx, 'Fs_list'}{1};
    rt = results{idx, 'rt_clean_list'}{1};
    rt_off = results{idx, 'rt_off_clean_list'}{1};
    rt_latency = results{idx, 'rt_latency_clean_list'}{1};
    
    start_timepoint = 1;
    end_timepoint = length(vpath);
    
    permutation = results{idx, 'permutation_list'};    
    ordered_vpath = vpath;
    for state = 1:K
        ordered_vpath(vpath == permutation(state)) = state;
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
    
    % Prevent epoch exceeding vpath dimension
    removeIndex = any(rt_latency + epoch_start_offset_251 < start_timepoint, 1)...
        | any(rt_latency + rt_off + epoch_end_offset_254 > end_timepoint, 1);    
    rt(removeIndex) = [];
    rt_off(removeIndex) = [];
    rt_latency(removeIndex) = [];
    
    all_session_RT = cat(2, all_session_RT, rt);
    all_session_RT_latency = cat(2, all_session_RT_latency, rt_latency);
    all_session_RT_off = cat(2, all_session_RT_off, rt_off);
    
    % Add 251/252
    state_by_epoch = epochByEvent(ordered_vpath, rt_latency, epoch_start_offset_251, epoch_end_offset_251);   
    all_session_251 = cat(1, all_session_251, state_by_epoch); 
    
    % Add 253
    state_by_epoch = epochByEvent(ordered_vpath, rt_latency+rt, epoch_start_offset_253, epoch_end_offset_253);   
    all_session_253 = cat(1, all_session_253, state_by_epoch); 
    
    % Add 254
    state_by_epoch = epochByEvent(ordered_vpath, rt_latency+rt_off, epoch_start_offset_254, epoch_end_offset_254);   
    all_session_254 = cat(1, all_session_254, state_by_epoch);   
end

[sortedRT, sortIdx] = sort(all_session_RT, 'descend');
sortedRT_off = all_session_RT_off(:, sortIdx);
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