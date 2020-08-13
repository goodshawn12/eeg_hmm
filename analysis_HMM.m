addpath('/home/ting/Documents/eeglab');
addpath('/home/ting/Documents/eeg_hmm');
addpath(genpath('/home/ting/Documents/eeg_hmm/Utility'));
addpath(genpath('/home/ting/Documents/eeg_hmm/HMM-MAR'));
addpath('/data/projects/Shawn/2016 JNE/dataset/');
addpath('/data/projects/Shawn/2019_HMM/data/');

%% ### This pipeline only supports results with the same K (number of inferred states)
K = 3;
method = 'GAU';
resultsDir = sprintf('/home/ting/Documents/eeg_hmm/HMM_results/results_K%d', K);
cd(resultsDir);

%% Manual setup files to process
% Selected sessions
% data_filelist = [dir('GAU*m*1.mat');dir('GAU*ASR20*0.mat');dir('MAR*ASR20*0.mat')];
% All session data
data_filelist = dir(strcat(method, '*_2.mat'));

filename_list = cell(length(data_filelist), 1);
for i = 1:length(data_filelist)
    file = data_filelist(i);
    if ~contains(file.name, 'ASR20')
        filename_list{i} = file.name;
    end
end
filename_list = filename_list(~cellfun('isempty', filename_list));
filename_list = unique(filename_list, 'stable');

n_of_files = length(filename_list)

%% Prepare raw data
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
fehist_list = cell(n_of_files,1);
vars_to_load = {'hmm', 'training_data_size', 'select_start', 'select_end', 'Gamma', 'vpath', 'fehist'};

for idx = 1:n_of_files
    try
        load(filename_list{idx}, vars_to_load{:});
    catch error
        disp(strcat('error loading', {' '}, filename_list{idx}));
        continue;
    end
    
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
    fehist_list{idx} = fehist;
end

hmm.train

results = table(filename_list, session_name_list, K_list, epoch_length_list, training_data_size_list, ...
    timepoint_start_list, timepoint_end_list, Fs_list, hmm_list, Gamma_list, vpath_list, fehist_list);

fprintf('Raw results data loaded.\n')

%% Examine free energy history
% isVisible = 'on';
% final_fehist = zeros(n_of_files,1);
% 
% figure('Visible', isVisible);
% for idx = 1:n_of_files
%     fehist = real(results{idx, 'fehist_list'}{1});
%     final_fehist(idx) = fehist(end);
%     plot(fehist);
%     hold on;
% end
% legend();
% 
% final_fehist
% 
% keyboard
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
% ends

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

fprintf('EEG meta data loaded.\n')

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

%% Global movmean correlation
% movmean_state_r_list = zeros(n_of_files, K);
% 
% for idx = 1:n_of_files
%     Fs = results{idx, 'Fs_list'}{1};
%     Gamma = results{idx, 'Gamma_list'}{1};
%     rs = transpose(results{idx, 'rt_speed_clean_list'}{1});
%     rt_latency = results{idx, 'rt_latency_clean_list'}{1};    
%     
%     win_len = Fs * 30;
%     smoothed_Gamma = movmean(Gamma, win_len, 1);
%     
%     rs_valid_range = rt_latency(1):rt_latency(end);
%     rs_interpolation = transpose(interp1(rt_latency, rs, rs_valid_range));
%     smoothed_rs = movmean(rs_interpolation, win_len);
%     
%     for state = 1:K
%         movmean_state_r_list(idx, state) = corr(smoothed_Gamma(rs_valid_range,state), smoothed_rs);
%     end    
% end
% 
% % The sort is in ascending order so the first index in the permutation_list
% % marks the lowest correlation (drowsy model) 
% [~, movmean_permutation_list] = sort(movmean_state_r_list, 2);
% movmean_inverse_permutation_list = zeros(size(movmean_permutation_list));
% for row=1:size(movmean_permutation_list,1)
%     movmean_inverse_permutation_list(row, movmean_permutation_list(row,:)) = 1:size(movmean_permutation_list,2);
% end
% 
% results = deleteVars(results, {'movmean_state_r_list', 'movmean_inverse_permutation_list'});
% results = addvars(results, movmean_state_r_list, movmean_inverse_permutation_list);

%% Global trend correlation
% trend_state_r_list = zeros(n_of_files, K);
% trial_Gamma_mean_list = cell(n_of_files, 1);
% 
% poly_deg = 6;
% 
% for idx = 1:n_of_files
%     Gamma = results{idx, 'Gamma_list'}{1};
%     rs = transpose(results{idx, 'rt_speed_clean_list'}{1});
%     rt_latency = results{idx, 'rt_latency_clean_list'}{1};    
%     
%     Gamma_trend = Gamma - detrend(Gamma, poly_deg);
%     
%     rs_valid_range = rt_latency(1):rt_latency(end);
%     rs_interpolation = transpose(interp1(rt_latency, rs, rs_valid_range));
%     rs_trend = rs_interpolation - detrend(rs_interpolation, poly_deg);
%     
%     for state = 1:K
%         trend_state_r_list(idx, state) = corr(Gamma_trend(rs_valid_range,state), rs_trend);
%     end    
% end
% 
% % The sort is in ascending order so the first index in the permutation_list
% % marks the lowest correlation (drowsy model) 
% [~, trend_permutation_list] = sort(trend_state_r_list, 2);
% trend_inverse_permutation_list = zeros(size(trend_permutation_list));
% for row=1:size(trend_permutation_list,1)
%     trend_inverse_permutation_list(row, trend_permutation_list(row,:)) = 1:size(trend_permutation_list,2);
% end
% 
% results = deleteVars(results, {'trend_state_r_list', 'trend_inverse_permutation_list'});
% results = addvars(results, trend_state_r_list, trend_inverse_permutation_list);

%% Examine fractional occupancy
fo_list = zeros(n_of_files, K);

for idx = 1:n_of_files
    Gamma = results{idx, 'Gamma_list'}{1};
    permutation = results{idx, 'smoothed_permutation_list'};
    
    fo = mean(Gamma, 1);
    fo_list(idx,:) = fo(permutation);
end
results = deleteVars(results, {'fo_list'});
results = addvars(results, fo_list);
fprintf('fractional  occupancy generated');

%% Calculate overall state probability
% overall_state_prob_list = zeros(n_of_files, K);
% for idx = 1:n_of_files
%     vpath = results{idx, 'vpath_list'}{1};
%     
%     for state = 1:K
%         overall_state_prob_list(idx, state) = sum(vpath==state)/length(vpath);
%     end
% end
% 
% results = deleteVars(results, 'overall_state_prob_list');
% results = addvars(results, overall_state_prob_list);

%% ######## Graphing Sections ########
%########################################
smooth_90s = 1;

% Preparation
if smooth_90s
    graphDir = sprintf('/home/ting/Documents/eeg_hmm/HMM_results/90s_corr_smoothed_graphs_K%d', K);
else
    graphDir = sprintf('/home/ting/Documents/eeg_hmm/HMM_results/graphs_K%d', K);
end
if ~exist(graphDir, 'dir')
    mkdir(graphDir);
end
cd(graphDir);

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
    if smooth_90s
        cmap_list(:,:,i) = cmap_list(results.smoothed_inverse_permutation_list(i,:),:,i);
    else
        cmap_list(:,:,i) = cmap_list(results.inverse_permutation_list(i,:),:,i);
    end
end
toSave = 0;
visible = 'off';

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

%% Epoch Gamma by timelock
Gamma_251_prior = cell(n_of_files, 1);
Gamma_251 = cell(n_of_files, 1);
Gamma_253 = cell(n_of_files, 1);
Gamma_254 = cell(n_of_files, 1);
epoched_rt = cell(n_of_files, 1);
epoched_rt_off = cell(n_of_files, 1);
epoched_rt_latency = cell(n_of_files, 1);

epoched_table = table();
smoothing_win_len = 1;

for idx = 1:n_of_files 
    Gamma = results{idx, 'Gamma_list'}{1};
    Gamma = movmean(Gamma, smoothing_win_len);
    Fs = results{idx, 'Fs_list'}{1};
    training_data_size = results{idx, 'epoch_length_list'}{1};
    rt = results{idx, 'rt_clean_list'}{1};
    rt_off = results{idx, 'rt_off_clean_list'}{1};
    rt_latency = results{idx, 'rt_latency_clean_list'}{1};
    rt = rt(2:end);
    rt_off = rt_off(2:end);
    rt_latency = rt_latency(2:end);
    
    zero_padded_Gamma = [zeros(training_data_size(1)-size(Gamma,1),K); Gamma];
    
    %Timelock prior to 251
    epoch_start_offset_251_prior = -8 * Fs;
    epoch_end_offset_251_prior = 4 * Fs;
    
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
    removeIndex = any(rt_latency + epoch_start_offset_251_prior < 1, 1)...
        | any(rt_latency + rt_off + epoch_end_offset_254 > size(zero_padded_Gamma, 1));
    rt(removeIndex) = [];
    rt_off(removeIndex) = [];
    rt_latency(removeIndex) = [];
    
    epoched_rt{idx} = rt;
    epoched_rt_off{idx} = rt_off;
    epoched_rt_latency{idx} = rt_latency;
   
    Gamma_251_prior{idx} = permute(epochByEvent(zero_padded_Gamma, rt_latency, epoch_start_offset_251_prior, epoch_end_offset_251_prior), [1 3 2]);
    Gamma_251{idx} = permute(epochByEvent(zero_padded_Gamma, rt_latency, epoch_start_offset_251, epoch_end_offset_251), [1 3 2]);
    Gamma_253{idx} = permute(epochByEvent(zero_padded_Gamma, rt_latency+rt, epoch_start_offset_253, epoch_end_offset_253), [1 3 2]);
    Gamma_254{idx} = permute(epochByEvent(zero_padded_Gamma, rt_latency+rt_off, epoch_start_offset_254, epoch_end_offset_254), [1 3 2]);
end

epoched_table = deleteVars(epoched_table, {'Gamma_251_prior', 'Gamma_251', 'Gamma_253', 'Gamma_254', 'epoched_rt', 'epoched_rt_off', 'epoched_rt_latency'});
epoched_table = addvars(epoched_table, Gamma_251_prior, Gamma_251, Gamma_253, Gamma_254, epoched_rt, epoched_rt_off, epoched_rt_latency);

fprintf('Gamma epoching done.\n')

%% Process Global Gamma epoched
total_epochs = 0;
for idx = 1:n_of_files
    total_epochs = total_epochs + length(epoched_table{idx, 'epoched_rt'}{1});
end
epoched_length = size(epoched_table{1, 'Gamma_251'}{1}, 2);

global_Gamma_251_prior = zeros(total_epochs, epoch_end_offset_251_prior-epoch_start_offset_251_prior+1, K);
global_Gamma_timelock = zeros(3, total_epochs, epoched_length, K);
global_rt = zeros(total_epochs, 1);
global_rt_off = zeros(total_epochs, 1);
global_rt_latency = zeros(total_epochs, 1);

curr_epoch = 1;
Gamma_events = {'Gamma_251', 'Gamma_253', 'Gamma_254'};
for idx = 1:n_of_files
    permutation = results{idx, 'smoothed_permutation_list'};
    rt = epoched_table{idx, 'epoched_rt'}{1};   
    n_of_epochs = length(rt);
    end_epoch = curr_epoch + n_of_epochs - 1;
    
    global_rt(curr_epoch:end_epoch) = rt;
    global_rt_off(curr_epoch:end_epoch) = epoched_table{idx, 'epoched_rt_off'}{1};
    global_rt_latency(curr_epoch:end_epoch) = epoched_table{idx, 'epoched_rt_latency'}{1};
    
    epoched_Gamma = epoched_table{idx, 'Gamma_251_prior'}{1};
    global_Gamma_251_prior(curr_epoch:end_epoch,:,:) = epoched_Gamma;
    for event = 1:3
        epoched_Gamma = epoched_table{idx, Gamma_events{event}}{1};
        global_Gamma_timelock(event, curr_epoch:end_epoch,:,:) = epoched_Gamma(:,:,permutation);
    end
    curr_epoch = end_epoch + 1;
end

% Sort everything by rt
[sorted_global_rt, sortIdx] = sort(global_rt);
sorted_global_rt_off = global_rt_off(sortIdx);
sorted_global_rt_latency = global_rt_latency(sortIdx);

sorted_global_Gamma_timelock = zeros(size(global_Gamma_timelock));
sorted_global_Gamma_251_prior = global_Gamma_251_prior(sortIdx,:,:);
for event = 1:3
    sorted_global_Gamma_timelock(event,:,:,:) = global_Gamma_timelock(event,sortIdx,:,:);
end

fprintf('Global Gamma and rt preparation done.\n')

%% Global colormap Gamma timelock prior to event
visible = 'on';
toSave = 0;

n_of_epochs = size(sorted_global_Gamma_timelock,2);
thick_line_width = 2;
thin_line_width = 0.01;
f = figure('Visible', visible); f.Position(3) = f.Position(4) * 3;
for state = 1:K
    gradient_cmap = ones(256,3);
    gradient_cmap(:,setdiff(1:K,state)) = repmat(((255:-1:0)/255)',1,2);
        
    epoched_Gamma = squeeze(sorted_global_Gamma_251_prior(:,:,state));
    ax = subplot(1,3,state);
    imagesc(epoched_Gamma); hold on; 
    axis('square');
    xtick_labels = cell(1, 12); xtick_labels_num = -8:1:4;
    for i = 1:13
        xtick_labels{i} = num2str(xtick_labels_num(i));
    end        
    xticks(1:Fs:Fs*12+1); xticklabels(xtick_labels); yticks([]); yticklabels([]);
    xlabel('seconds'), ylabel('epochs');
    
    line([-epoch_start_offset_251_prior, -epoch_start_offset_251_prior], [1, n_of_epochs], 'linewidth',2,'color','k');
    plot(sorted_global_rt-epoch_start_offset_251_prior, 1:n_of_epochs, 'linewidth',thick_line_width,'color',[0.5 0.5 0.5]); 
    colormap(ax, gradient_cmap)
end

if toSave
    set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(method, '_K', num2str(K), '_global_Gamma_timelock_251_prior_', num2str(smoothing_win_len), '.fig'))
    print(strcat(method, '_K', num2str(K), '_global_Gamma_timelock_251_prior_', num2str(smoothing_win_len)), '-djpeg')
end

fprintf('Global Gamma prior to 251 timelock plots done.\n')

%% Global colormap of state Gamma timelock to events
visible = 'on';
toSave = 0;

n_of_epochs = size(sorted_global_Gamma_timelock,2);
thick_line_width = 2;
thin_line_width = 0.01;
gray_color = [0.5 0.5 0.5];
f1 = figure('Visible', visible); f1.Position(3) = f1.Position(4) * 3;
for state = 1:K
    gradient_cmap = ones(256,3);
    gradient_cmap(:,setdiff(1:K,state)) = repmat(((255:-1:0)/255)',1,2);
    
    ax_251 = subplot(1, K, state);
    imagesc(squeeze(sorted_global_Gamma_timelock(1,:,:,state))); hold on;
    ax_251.XTick = 1:Fs:Fs*6+1; ax_251.XTickLabel = {'-2', '-1', '0', '1', '2', '3', '4'}; ax_251.YTick = []; ax_251.YTickLabel = [];
    axis('square');
    line([-epoch_start_offset_251, -epoch_start_offset_251], [1, n_of_epochs], 'linewidth',thick_line_width,'color','k');
    plot(sorted_global_rt-epoch_start_offset_251, 1:n_of_epochs, 'linewidth',thick_line_width,'color',gray_color);
%     plot(sorted_global_rt_off-epoch_start_offset_251, 1:n_of_epochs, 'linewidth',thin_line_width,'color', [0.8, 0.8, 0.8]);
    colormap(ax_251, gradient_cmap);
end

if toSave
    set(gcf, 'PaperPositionMode', 'auto');
%     saveas(gcf, strcat(method, '_K', num2str(K), '_global_Gamma_timelock_movmean_', num2str(smoothing_win_len), '.fig'))
    print(strcat(method, '_K', num2str(K), '_global_Gamma_timelock_251_movmean_', num2str(smoothing_win_len)), '-djpeg')
end

f2 = figure('Visible', visible); f2.Position(3) = f2.Position(4) * 3;
for state = 1:K
    gradient_cmap = ones(256,3);
    gradient_cmap(:,setdiff(1:K,state)) = repmat(((255:-1:0)/255)',1,2);
    
    ax_253 = subplot(1, K, state); 
    imagesc(squeeze(sorted_global_Gamma_timelock(2,:,:,state))); hold on;
    ax_253.XTick = 1:Fs:Fs*6+1; ax_253.XTickLabel = {'-3', '-2', '-1', '0', '1', '2', '3'}; ax_253.YTick = []; ax_253.YTickLabel = [];
    axis('square');
    plot(-epoch_start_offset_253 - sorted_global_rt, 1:n_of_epochs, 'linewidth',thick_line_width,'color',gray_color)
    line([-epoch_start_offset_253, -epoch_start_offset_253], [1, n_of_epochs],'linewidth',thick_line_width,'color','k')
%     plot(-epoch_start_offset_253 - sorted_global_rt+sorted_global_rt_off, 1:n_of_epochs, 'linewidth',thin_line_width,'color',[0.8,0.8,0.8])

    colormap(ax_253, gradient_cmap);
end

if toSave
    set(gcf, 'PaperPositionMode', 'auto');
%     saveas(gcf, strcat(method, '_K', num2str(K), '_global_Gamma_timelock_movmean_', num2str(smoothing_win_len), '.fig'))
    print(strcat(method, '_K', num2str(K), '_global_Gamma_timelock_253_movmean_', num2str(smoothing_win_len)), '-djpeg')
end

f3 = figure('Visible', visible); f3.Position(3) = f3.Position(4) * 3;
for state = 1:K
    gradient_cmap = ones(256,3);
    gradient_cmap(:,setdiff(1:K,state)) = repmat(((255:-1:0)/255)',1,2);
    
    ax_254 = subplot(1, K, state);
    imagesc(squeeze(sorted_global_Gamma_timelock(3,:,:,state))); hold on;
    ax_254.XTick = 1:Fs:Fs*6+1; ax_254.XTickLabel = {'-4', '-3', '-2', '-1', '0', '1', '2'}; ax_254.YTick = []; ax_254.YTickLabel = [];
    axis('square'); 
%     plot(-epoch_start_offset_254 - sorted_global_rt_off, 1:n_of_epochs, 'linewidth',thin_line_width,'color','k')
%     plot(-epoch_start_offset_254 - sorted_global_rt_off+sorted_global_rt, 1:n_of_epochs, 'linewidth',thin_line_width,'color','w')
    line([-epoch_start_offset_254, -epoch_start_offset_254], [1, n_of_epochs], 'linewidth',thick_line_width,'color', 'k')
    
    colormap(ax_254, gradient_cmap);
end

if toSave
    set(gcf, 'PaperPositionMode', 'auto');
%     saveas(gcf, strcat(method, '_K', num2str(K), '_global_Gamma_timelock_movmean_', num2str(smoothing_win_len), '.fig'))
    print(strcat(method, '_K', num2str(K), '_global_Gamma_timelock_254_movmean_', num2str(smoothing_win_len)), '-djpeg')
end

fprintf('Global Gamma timelock plots done.\n')

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
visible = 'off';
toSave = 1;
raw = 0;
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

for idx = 1:n_of_files
    if raw
        smoothed_vpath = results{idx, 'vpath_list'}{1};
    else
        smoothed_vpath = results{idx, 'smoothed_vpath_list'}{1};
    end    
    training_data_size = results{idx, 'epoch_length_list'}{1};
    start_timepoint = results{idx, 'timepoint_start_list'}{1};
    end_timepoint = results{idx, 'timepoint_end_list'}{1};
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

fprintf('timelock vpath plots done.\n')

%% Prepare to plot all session timelock
% all_session_251 = [];
% all_session_253 = [];
% all_session_254 = [];
% all_session_RT = [];
% all_session_RT_off = [];
% all_session_RT_latency = [];
% 
% for idx = 1:n_of_files
%     smoothed_vpath = results{idx, 'smoothed_vpath_list'}{1};
%     training_data_size = results{idx, 'epoch_length_list'}{1};
%     start_timepoint = results{idx, 'timepoint_start_list'}{1};
%     end_timepoint = results{idx, 'timepoint_end_list'}{1};
%     all_events = results{idx, 'events_list'}{1};
%     Fs = results{idx, 'Fs_list'}{1};
%     rt = results{idx, 'rt_clean_list'}{1};
%     rt_off = results{idx, 'rt_off_clean_list'}{1};
%     rt_latency = results{idx, 'rt_latency_clean_list'}{1};
%     permutation = results{idx, 'permutation_list'};
%     filename = results{idx, 'filename_list'}{1}
%     
%     zero_padded_vpath = [zeros(training_data_size(1)-size(smoothed_vpath,1),1); smoothed_vpath];
%     ordered_vpath = zero_padded_vpath;
%     ordered_vpath(zero_padded_vpath == permutation(1)) = 1;
%     ordered_vpath(zero_padded_vpath == permutation(2)) = 2;
%     ordered_vpath(zero_padded_vpath == permutation(3)) = 3;
%         
%     % Timelock 251/252, lane-departure task introduced
%     epoch_start_offset_251 = -2 * Fs;
%     epoch_end_offset_251 = 4 * Fs;
%     
%     % Timelock 253, response onset
%     epoch_start_offset_253 = -3 * Fs;
%     epoch_end_offset_253 = 3 * Fs;
%     
%     % Timelock 254, response off-set
%     epoch_start_offset_254 = -4 * Fs;
%     epoch_end_offset_254 = 2 * Fs;
%     
%     % Prevent epoch exceeding vpath dimension
%     removeIndex = any(rt_latency + epoch_start_offset_251 < start_timepoint, 1)...
%         | any(rt_latency + rt_off + epoch_end_offset_254 > end_timepoint, 1);    
%     rt(removeIndex) = [];
%     rt_off(removeIndex) = [];
%     rt_latency(removeIndex) = [];
%     [sortedRT, sortIdx] = sort(rt);
%     sortedRT_off = rt_off(:, sortIdx);
%     sortedRT_latency = rt_latency(:, sortIdx);
%     
%     all_session_RT = cat(2, all_session_RT, rt);
%     all_session_RT_latency = cat(2, all_session_RT_latency, sortedRT_latency);
%     all_session_RT_off = cat(2, all_session_RT_off, sortedRT_off);
%     
%     % Add 251/252
%     state_by_epoch = epochByEvent(ordered_vpath, rt_latency, epoch_start_offset_251, epoch_end_offset_251);   
%     state_by_epoch = state_by_epoch(sortIdx,:);
%     all_session_251 = cat(1, all_session_251, state_by_epoch); 
%     
%     % Add 253
%     state_by_epoch = epochByEvent(ordered_vpath, rt_latency+rt, epoch_start_offset_253, epoch_end_offset_253);   
%     state_by_epoch = state_by_epoch(sortIdx,:);
%     all_session_253 = cat(1, all_session_253, state_by_epoch); 
%     
%     % Add 254
%     state_by_epoch = epochByEvent(ordered_vpath, rt_latency+rt, epoch_start_offset_253, epoch_end_offset_253);   
%     state_by_epoch = state_by_epoch(sortIdx,:);
%     all_session_254 = cat(1, all_session_254, state_by_epoch);   
% end
% 
% [sortedRT, sortIdx] = sort(all_session_RT);
% sortedRT_off = all_session_RT_off(:, sortIdx);
% sortedRT_latency = all_session_RT_latency(:, sortIdx);
% all_session_251 = all_session_251(sortIdx, :);
% all_session_253 = all_session_253(sortIdx, :);
% all_session_254 = all_session_254(sortIdx, :);

%% Plot all session timelock
% visible = 'on';
% toSave = 0;
% cmap = [1,0.1,0.1; 0.1,1,0.1; 0.1,0.1,1];
% 
% figure('Visible', visible);
% imagesc(all_session_251);
% set(gca,'Title',text('String','Time locked to lane-departure event'));
% set(gca,'XLabel',text('String','Time offset (sec)'));
% set(gca,'YLabel',text('String','Epochs'));
% xticks(1:Fs:Fs*6+1);
% xticklabels({'-2', '-1', '0', '1', '2', '3', '4'});
% colormap(cmap);
% 
% hold on,
% line([-epoch_start_offset_251, -epoch_start_offset_251], [1, length(sortedRT)], 'linewidth',2,'color','k');
% plot(sortedRT-epoch_start_offset_251, 1:length(sortedRT), 'linewidth',2,'color','w');
% plot(sortedRT_off-epoch_start_offset_251, 1:length(sortedRT), 'linewidth',0.01,'color', [0.8, 0.8, 0.8]);
% 
% if toSave
%     set(gcf, 'PaperPositionMode', 'auto');
%     saveas(gcf, strcat('all_session_251','.fig'))
%     print('all_session_251', '-djpeg')  
% end
% 
% % Plot 253
% figure('Visible', visible);
% imagesc(all_session_253);
% set(gca,'Title',text('String','Time locked to car driver response start'));
% set(gca,'XLabel',text('String','Time offset (sec)'));
% set(gca,'YLabel',text('String','Epochs'));
% xticks(1:Fs:Fs*6+1);
% xticklabels({'-3', '-2', '-1', '0', '1', '2', '3'});
% colormap(cmap);
% 
% hold on,
% plot(-epoch_start_offset_253 - sortedRT, 1:length(sortedRT), 'linewidth',2,'color','k')
% line([-epoch_start_offset_253, -epoch_start_offset_253], [1, length(sortedRT)],'linewidth',2,'color','w')
% plot(-epoch_start_offset_253 - sortedRT+sortedRT_off, 1:length(sortedRT), 'linewidth',0.01,'color',[0.8,0.8,0.8])
% 
% if toSave
%     set(gcf, 'PaperPositionMode', 'auto');
%     saveas(gcf, strcat('all_session_253','.fig'))
%     print('all_session_253', '-djpeg')
% end
% 
% % Plot 254
% figure('Visible', visible);
% imagesc(all_session_254);
% set(gca,'Title',text('String','Time locked to driver response end'));
% set(gca,'XLabel',text('String','Time offset (sec)'));
% set(gca,'YLabel',text('String','Epochs'));
% xticks(1:Fs:Fs*6+1);
% xticklabels({'-4', '-3', '-2', '-1', '0', '1', '2'});
% colormap(cmap);
% 
% hold on,
% plot(-epoch_start_offset_254 - sortedRT_off, 1:length(sortedRT), 'linewidth',0.01,'color','k')
% plot(-epoch_start_offset_254 - sortedRT_off+sortedRT, 1:length(sortedRT), 'linewidth',0.01,'color','w')
% line([-epoch_start_offset_254, -epoch_start_offset_254], [1, length(sortedRT)], 'linewidth',2,'color',[0.8,0.8,0.8])
% 
% if toSave
%     set(gcf, 'PaperPositionMode', 'auto');
%     saveas(gcf, strcat('all_session_254','.fig'))
%     print('all_session_254', '-djpeg')
% end

%% Plot smoothed state time 
visible = 'on';
toSave = 0;
smooth_90s = 1;
for idx = 9
    Fs = results{idx, 'Fs_list'}{1};
    Gamma = results{idx, 'Gamma_list'}{1};
    rt = results{idx, 'rt_clean_list'}{1};
    rt_latency = results{idx, 'rt_latency_clean_list'}{1};
    rt_speed = results{idx, 'rt_speed_clean_list'}{1};
    
    % Preparation
    if smooth_90s
        state_r_list = results{idx, 'smoothed_state_r_list'};
    else
        state_r_list = results{idx, 'state_r_list'};
    end
    
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
    stateProb(:, any(sum(stateProb, 1)==0, 1)) = [];
    end_minutes = (size(stateProb,2)-1)*walkLen/Fs/60;
    
    n_subplots = K+2;
    fig1 = figure('Visible', visible);
    p = uipanel('Parent', fig1);
    
    [~, corresponding_vpath] = max(stateProb, [], 1);
    vpath_values = unique(corresponding_vpath);
    n_of_vpath_values = length(vpath_values);
    vpath_axis = subplot(n_subplots,1,1, 'Parent', p);
    vpath_axis.Position(4) = vpath_axis.Position(4)/2;
    imagesc([0 end_minutes], [1 1], corresponding_vpath); colormap(vpath_axis, cmap(vpath_values, :));
    title('Viterbi Path')
    vpath_axis.YTick = [];
    vpath_axis.YTickLabel = [];
    vpath_axis.XTick = [];
    vpath_axis.XTickLabel = []; 
    
    % Configure colorbar
%     ytick_space = (n_of_vpath_values-1)/2/n_of_vpath_values;
%     yticks = 1-ytick_space:2*ytick_space:n_of_vpath_values-0.01;
%     yticks(1) = [];
%     ytickslabel = cell(1,n_of_vpath_values);
%     for state = vpath_values
%         ytickslabel{state} = num2str(state);
%     end
%     colorbar('YTick', yticks, 'YTickLabel', ytickslabel);   

    rs_axis = subplot(n_subplots,1,2, 'Parent', p);
    rs_x_values = rt_latency/Fs/60;
    plot(rs_x_values, rt_speed, '-', 'LineWidth', 0.2, 'Color', [0.7,0.7,0.7]); hold on;
    plot(rs_x_values, rt_speed, 'k.', 'LineWidth', 2, 'Color', 'k'); 
    title('Reaction Speed by Events')
    ylabel('Speed (1/sec)')
    xlim([0 end_minutes])
    rs_axis.XTick = [];
    rs_axis.XTickLabel = [];
%     xlim([0, vpath_axis.XLim(2)*rs_axis.Position(3)/vpath_axis.Position(3)])
    
    [~, permutation_list] = sort(state_r_list);
    for i = 1:K
        state = permutation_list(i);
        state_axis = subplot(n_subplots,1,i+2);
        x_values = 1:(end_minutes-1)/(size(stateProb,2)-1):end_minutes;
        y_values = stateProb(state,:);
                
        plot(x_values, y_values, 'Color', cmap(state,:))
        r = state_r_list(state);
        title(strcat('State ', num2str(state)))
        ylabel('Probability')
        txt = ['r = ' num2str(r)];
        xlim([0 end_minutes])
% %         xlim([0, vpath_axis.XLim(2)*state_axis.Position(3)/vpath_axis.Position(3)])
        if mean(y_values(1:floor(size(y_values,2)/10))) > mean(state_axis.YLim)
            text_y_position = state_axis.YLim(1) + range(state_axis.YLim)/5;
        else
            text_y_position = state_axis.YLim(2) - range(state_axis.YLim)/5;
        end    
        text(2, text_y_position, txt, 'FontSize', 10, 'Parent', state_axis)
    end
    xlabel('Time (minutes)')
    
    if toSave
        output_filename = strcat(filename,'_','path');
        set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(output_filename,'.fig'))
        print(output_filename, '-djpeg')
    end
end

fprintf('vpath plot done.\n')

% %% Table view overall state probability
% figure;
% t = uitable('Data', results.overall_state_prob_list, 'ColumnName', {'State 1', 'State 2', 'State 3'}, 'RowName', results.filename_list);
% t.Position(3:4) = t.Extent(3:4);

%% Plot state timepath with trend
visible = 'on';
toSave = 0;
smooth_90s = 1;
for idx = 11 %1:n_of_files
    Fs = results{idx, 'Fs_list'}{1};
    Gamma = results{idx, 'Gamma_list'}{1};
    vpath = results{idx, 'vpath_list'}{1};
    rt = results{idx, 'rt_clean_list'}{1};
    rt_latency = results{idx, 'rt_latency_clean_list'}{1};
    rt_speed = results{idx, 'rt_speed_clean_list'}{1};
    
    % Preparation
    if smooth_90s
        state_r_list = results{idx, 'smoothed_state_r_list'};
    else
        state_r_list = results{idx, 'state_r_list'};
    end
    
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
    stateProb(:, any(sum(stateProb, 1)==0, 1)) = [];
    end_minutes = (size(stateProb,2)-1)*walkLen/Fs/60;
    
    n_subplots = K+3;
    fig1 = figure('Visible', visible); fig1.Position(3) = fig1.Position(3) * 1.5; fig1.Position(4) = fig1.Position(4) * 1.5;
    p = uipanel('Parent', fig1);
    
    [~, corresponding_vpath] = max(stateProb, [], 1);
    vpath_values = unique(corresponding_vpath);
    n_of_vpath_values = length(vpath_values);
    vpath_axis = subplot(n_subplots,1,1, 'Parent', p);    
    vpath_axis.Position(4) = vpath_axis.Position(4)/2;
    imagesc([0 end_minutes], [1 1], corresponding_vpath); colormap(vpath_axis, cmap(vpath_values, :));
    title('Viterbi Path')
    vpath_axis.YTick = [];
    vpath_axis.YTickLabel = [];
    vpath_axis.XTick = [];
    vpath_axis.XTickLabel = [];    
    
    rs_axis = subplot(n_subplots,1,2, 'Parent', p);
    rs_x_values = rt_latency/Fs/60;
    plot(rs_x_values, rt_speed, '-', 'LineWidth', 0.2, 'Color', [0.7,0.7,0.7]); hold on;
    plot(rs_x_values, rt_speed, 'k.', 'LineWidth', 2, 'Color', 'k');
    title('Reaction Speed by Events')
    ylabel('Speed (1/sec)')
    xlim([0 end_minutes])
    rs_axis.XTick = [];
    rs_axis.XTickLabel = [];
%     xlim([0, vpath_axis.XLim(2)*rs_axis.Position(3)/vpath_axis.Position(3)])
    
    [~, permutation_list] = sort(state_r_list);
    for i = 1:K
        state = permutation_list(i);
        state_axis = subplot(n_subplots,1,i+2);
        x_values = 1:(end_minutes-1)/(size(stateProb,2)-1):end_minutes;
        y_values = stateProb(state,:);
        
        plot(x_values, y_values, 'Color', cmap(state,:))
        r = state_r_list(state);       
        txt = ['r = ' num2str(r)];    
        if mean(y_values(1:floor(size(y_values,2)/10))) > mean(state_axis.YLim)
            text_y_position = state_axis.YLim(1) + range(state_axis.YLim)/5;
        else
            text_y_position = state_axis.YLim(2) - range(state_axis.YLim)/5;
        end
        text_label = text(1, text_y_position, txt, 'FontSize', 10, 'Parent', state_axis);
        text_label.Position(2) = state_axis.YLim(1) + range(state_axis.YLim)/5;
        
        title(strcat('State', {' '}, num2str(state)))
        ylabel('Probability')
        xlim([0 end_minutes])
        state_axis.XTick = [];
        state_axis.XTickLabel = [];
    end
    
    %     detrended_Gamma = detrend(Gamma, 1, 1:floor(length(Gamma)/10):length(Gamma));
%     detrended_Gamma = detrend(Gamma, 6);
%     Gamma_trend = Gamma - detrended_Gamma;
%     detrended_RS = detrend(rt_speed, 6);
%     RS_trend = rt_speed - detrended_RS;
    
%     stateTrend = zeros(nStates,floor(nTime/walkLen));
%     for it = 1:floor((nTime-winLen+walkLen)/walkLen)
%         time_range = (it-1)*walkLen+1:(it-1)*walkLen+winLen;
%         stateTrend(:,it) = mean(Gamma_trend(time_range,:),1)';
%     end
%     stateTrend(:, any(sum(stateTrend, 1)==0, 1)) = [];

    smoothing_win_len = Fs * 300;
    smoothed_Gamma = movmean(Gamma, smoothing_win_len, 1);    
    rs_valid_range = rt_latency(1):rt_latency(end);
    rs_interpolation = transpose(interp1(rt_latency, rt_speed, rs_valid_range));
    smoothed_rs = movmean(rs_interpolation, smoothing_win_len);
    
    trend_axis = subplot(n_subplots,1,K+3);
    x_values = rs_valid_range; %1:(end_minutes-1)/(size(stateTrend,2)-1):end_minutes;
    plot(rs_valid_range, normalize(smoothed_rs), 'Color', [0.7,0.7,0.7], 'LineWidth', 3, 'LineStyle', ':'); hold on;
    for i = 1:K
        y_values = smoothed_Gamma(rs_valid_range,i);        
        plot(rs_valid_range, normalize(y_values), 'Color', cmap(i,:), 'LineWidth', 1.5); hold on;          
    end
    xlim([0 end_minutes*60*Fs])
    trend_axis.Position(4) = trend_axis.Position(4) * 3/2;
    trend_axis.Position(2) = state_axis.OuterPosition(2) - trend_axis.OuterPosition(4);
    trend_axis.YTick = [];
    trend_axis.YLabel = [];
    title('Normalized RS and state trends')
    xlabel('Time (minutes)')
%     legend('RS trend')
    
    if toSave
        output_filename = strcat(filename,'_','trend_path');
        set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(output_filename,'.fig'))
        print(output_filename, '-djpeg')
    end
   
end

fprintf('detrend plot done.\n')

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
%         saveas(gcf, strcat(output_filename,'.fig'))
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
%             saveas(cgfighandle, strcat(output_filename,'.fig'))
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

%% ##### Spectra section #####
% #############################

%% Load EEG signal (Only needed for computing Gaussian HMM spectra)
eeglab;
X_list = cell(n_of_files, 1);
T_list = cell(n_of_files, 1);

% For all sessions
for idx = 1:n_of_files
    filename = split(filename_list{idx}, '_');
    filename = join(filename(2:length(filename)-1), '_');
    filename = strcat(filename, '.set');
    
    datafile = pop_loadset(filename{1});    
    X_list{idx} = transpose(datafile.data);
    T_list{idx} = datafile.pnts;
end

results = deleteVars(results, {'X_list', 'T_list'});
results = addvars(results, X_list, T_list);

%% Compute spectra
delete(gcp('nocreate')); % shut down any current pool
npar = 27;
parpool(npar);   % request workers from the cluster

spectra_list = cell(n_of_files, 1);
Gamma_list = results.Gamma_list;
hmm_list = results.hmm_list;
X_list = results.X_list;
T_list = results.T_list;
options = struct();

if strcmp('GAU', method)
    options.order = 0;    
elseif strcmp('MAR', method)
    hmm = hmm_list{1};   
    options.order = hmm.train.order;
end
options.Fs = results{1, 'Fs_list'}{1};

parfor (idx = 1:n_of_files, npar)
    if strcmp('GAU', method)
        spectra_list{idx} = hmmspectramt(X_list{idx}, T_list{idx}, Gamma_list{idx}, options);
    elseif strcmp('MAR', method)
        spectra_list{idx} = hmmspectramar([], [], hmm_list{idx}, [], options);
    end    
end

results = deleteVars(results, 'spectra_list');
results = addvars(results, spectra_list);
delete(gcp('nocreate'));

%% Save spectra to file
destination_dir = strcat('/home/ting/Documents/eeg_hmm/HMM_results/results_K', num2str(K), '/spectra');
if ~exist(destination_dir, 'dir')
    mkdir(destination_dir);
end

for idx = 1:n_of_files
    filename = results{idx, 'filename_list'}{1};
    filename = split(filename, '.');
    filename = filename{1};
    spectra = results{idx, 'spectra_list'}{1};
    
    output_filename = strcat(filename, '_spectra.mat')
    save(strcat('/home/ting/Documents/eeg_hmm/HMM_results/results_K', num2str(K), '/spectra/', output_filename), '-struct', 'spectra');
end

fprintf('Spectra files saved.\n')

%% Plot state spectra
toSave = 1;
isVisible = "off";

for idx = 1:n_of_files
    spectra = results{idx, 'spectra_list'}{1};
    [~, f] = getPSD(spectra, 1);
    filename = results{idx, 'filename_list'}{1};
    filename = split(filename, '.');
    filename = filename{1}
    
    figure('Visible', isVisible);
    for state = 1:K
        psd = getPSD(spectra, state);
        mean_log_psd = log10(mean(psd, 2));
        plot(mean_log_psd, 'color', cmap(state,:)); hold on;
    end
    xlabel('frequency (Hz)'); ylabel('log10 power');
    legend(state_description);
    
    if toSave
        output_filename = strcat(filename, '_', 'StateSpectra');
        set(gcf, 'PaperPositionMode', 'auto');
%             saveas(cgfighandle, strcat(output_filename,'.fig'))
        print(output_filename, '-djpeg')
    end
end

fprintf('state spectra plot done.\n')

%% Colormap by epochs timelock at events
visible = 'off';
toSave = 1;

% Timelock 251/252, lane-departure task introduced
epoch_start_offset_251 = -2 * Fs;
epoch_end_offset_251 = 4 * Fs;

% Timelock 253, response onset
epoch_start_offset_253 = -3 * Fs;
epoch_end_offset_253 = 3 * Fs;

% Timelock 254, response off-set
epoch_start_offset_254 = -4 * Fs;
epoch_end_offset_254 = 2 * Fs;

for idx = 1:n_of_files
    Gamma = results{idx, 'Gamma_list'}{1};   
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
    
    zero_padded_Gamma = [zeros(training_data_size(1)-size(Gamma,1), size(Gamma,2)); Gamma];
    cmap = cmap_list(:,:,idx);
    xticks_value = 0:Fs:6*Fs; xticks_value(1) = 1;
    
    % Prevent epoch exceeding data dimension
    removeIndex = any(rt_latency + epoch_start_offset_251 < start_timepoint, 1)...
        | any(rt_latency + rt_off + epoch_end_offset_254 > end_timepoint, 1);    
    rt(removeIndex) = [];
    rt_off(removeIndex) = [];
    rt_latency(removeIndex) = [];    
    
    % Plot 251/252
    state_by_epoch = epochByEvent(Gamma, rt_latency, epoch_start_offset_251, epoch_end_offset_251 - 1);    
    
    figure('Visible', visible);
    line_list = cell(state, 1);
    for state = 1:K
        state_course = squeeze(state_by_epoch(:,state,:));
%         sem_value = std(state_course)' ./ sqrt(size(state_course,1));
        mean_state_course = mean(state_course, 1);        
%         x_value = (1:length(std_value))';
        
%         error_band = fill([x_value; flipud(x_value)], [mean_state_course' + sem_value; flipud(mean_state_course' - sem_value)], ...
%             cmap(state,:), 'LineStyle', 'none'); hold on;
%         alpha(error_band, .1); 
        spectra_line = plot(mean(state_course, 1), 'color', cmap(state,:)); hold on;
        line_list{state} = spectra_line;
    end
%     for state = 1:K
%         uistack(line_list{state}, 'top');
%     end
    hold on;
    line([-epoch_start_offset_251, -epoch_start_offset_251], [0, 1], 'linewidth', 2, 'color','k');
    
    set(gca,'Title',text('String','Time locked to lane-departure event'));
    set(gca,'XLabel',text('String','Time offset (sec)'));
    set(gca,'YLabel',text('String','Mean Probability over Trials'));
    set(gca, 'XTick', xticks_value, 'XTickLabel', {'-2', '-1', '0', '1', '2', '3', '4'});  
  
    if toSave
        set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(filename,'_Gamma_mean_','251','.fig'))
        print(strcat(filename,'_Gamma_mean_','251'), '-djpeg')  
    end
    
    % Plot 253
    state_by_epoch = epochByEvent(Gamma, rt_latency + rt, epoch_start_offset_253, epoch_end_offset_253 - 1);    
    
    figure('Visible', visible);
    for state = 1:K
        mean_state_course = mean(squeeze(state_by_epoch(:,state,:)), 1);
        plot(mean_state_course, 'color', cmap(state,:)); hold on;
    end
    hold on;
    line([-epoch_start_offset_253, -epoch_start_offset_253], [0, 1], 'linewidth', 2, 'color', 'k');

    set(gca,'Title',text('String','Time locked to car driver response start'));
    set(gca,'XLabel',text('String','Time offset (sec)'));
    set(gca,'YLabel',text('String','Mean Probability over Trials'));
    set(gca, 'XTick', xticks_value, 'XTickLabel', {'-3', '-2', '-1', '0', '1', '2', '3'}) 

    if toSave
        set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(filename,'_Gamma_mean_','253','.fig'))
        print(strcat(filename,'_Gamma_mean_','253'), '-djpeg')
    end

    % Plot 254
    state_by_epoch = epochByEvent(Gamma, rt_latency + rt_off, epoch_start_offset_254, epoch_end_offset_254 - 1);    
    
    figure('Visible', visible);
    for state = 1:K
        mean_state_course = mean(squeeze(state_by_epoch(:,state,:)), 1);
        plot(mean_state_course, 'color', cmap(state,:)); hold on;
    end
    hold on;
    line([-epoch_start_offset_254, -epoch_start_offset_254], [0, 1], 'linewidth', 2, 'color', 'k');

    set(gca,'Title',text('String','Time locked to car driver response end'));
    set(gca,'XLabel',text('String','Time offset (sec)'));
    set(gca,'YLabel',text('String','Mean Probability over Trials'));
    set(gca, 'XTick', xticks_value, 'XTickLabel', {'-4', '-3', '-2', '-1', '0', '1', '2'}) 

    if toSave
        set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(filename,'_Gamma_mean_','254','.fig'))
        print(strcat(filename,'_Gamma_mean_','254'), '-djpeg')
    end

end

fprintf('Gamma timelock plots done.\n')

%% Calculate global averaged functional connectivity
% set(0, 'ShowHiddenHandles', 'on');
common_chanlocs = results{1, 'chanlocs_list'}{1};
for idx = 2:n_of_files
    common_chanlocs = intersect(common_chanlocs, results{idx, 'chanlocs_list'}{1}, 'stable');
end

mean_covmat = zeros(length(common_chanlocs), length(common_chanlocs), K);
mean_corrmat = zeros(length(common_chanlocs), length(common_chanlocs), K);

for idx = 1:n_of_files
    hmm = results{idx, 'hmm_list'}{1};
    chanlocs = results{idx, 'chanlocs_list'}{1};
    [~, permutation] = sort(results{idx, 'smoothed_state_r_list'}, 2);

    [~, chan_index] = ismember(common_chanlocs, chanlocs);
    for i = 1:K
        state = permutation(i);
        [covmat, corrmat] = getFuncConn(hmm, state);
        mean_covmat(:,:,i) = covmat(chan_index, chan_index) + mean_covmat(:,:,i);
        mean_corrmat(:,:,i) = corrmat(chan_index, chan_index) + mean_corrmat(:,:,i);
    end
end
mean_covmat = mean_covmat ./ n_of_files;
mean_corrmat = mean_corrmat ./ n_of_files;

%% Plot covmat
isVisible = "on";
toSave = 1;

if K == 2
    n_row = 1; n_col = 2;
elseif K == 3
    n_row = 1; n_col = 3;
elseif K == 4
    n_row = 2; n_col = 2;
else
    n_row = ceil(K/3); n_col = 3;
end

% Plot covmat
global_min = min(mean_covmat, [], 'all');
global_max = max(mean_covmat, [], 'all');
f = figure('Visible', isVisible); colormap(copper);
for i = 1:K
    caxis manual; caxis([global_min, global_max]);
    state_axes = subplot(n_row, n_col, i);
    imagesc(flip(mean_covmat(:,:,i), 1));
    xticks(1:length(common_chanlocs)); xticklabels(common_chanlocs); xtickangle(90);
    yticks(1:length(common_chanlocs)); yticklabels(flip(common_chanlocs));
    title(state_description{i});
    axis('square');
end        
f.Position(3) = f.Position(4) * 3;
colorbar_x = state_axes.OuterPosition(1) + state_axes.OuterPosition(3);
colorbar_y = state_axes.Position(2);
colorbar_w = 0.01;
colorbar_h = state_axes.Position(4);
colorbar('Location', 'manual', 'Position', [colorbar_x, colorbar_y, colorbar_w, colorbar_h]);

if toSave
    set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(filename,'_Gamma_mean_','254','.fig'))
    print(strcat(method,'_K',num2str(K),'_meanCovmat_smoothed'), '-djpeg')
end

%% Plot corrmat
isVisible = "on";
toSave = 1;

if K == 2
    n_row = 1; n_col = 2;
elseif K == 3
    n_row = 1; n_col = 3;
elseif K == 4
    n_row = 2; n_col = 2;
else
    n_row = ceil(K/3); n_col = 3;
end

% Plot covmat
global_min = min(mean_corrmat, [], 'all');
global_max = max(mean_corrmat, [], 'all');
f = figure('Visible', isVisible); colormap(copper);
for i = 1:K
    caxis manual; caxis([global_min, global_max]);
    state_axes = subplot(n_row, n_col, i);
    mean_corrmat(1:5,1:5,i)
    imagesc(flip(mean_corrmat(:,:,i), 1));
    xticks(1:length(common_chanlocs)); xticklabels(common_chanlocs); xtickangle(90);
    yticks(1:length(common_chanlocs)); yticklabels(flip(common_chanlocs));
    title(state_description{i});
    axis('square');
end        
f.Position(3) = f.Position(4) * 3;
colorbar_x = state_axes.OuterPosition(1) + state_axes.OuterPosition(3);
colorbar_y = state_axes.Position(2);
colorbar_w = 0.01;
colorbar_h = state_axes.Position(4);
colorbar('Location', 'manual', 'Position', [colorbar_x, colorbar_y, colorbar_w, colorbar_h]);

if toSave
    set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, strcat(filename,'_Gamma_mean_','254','.fig'))
    print(strcat(method,'_K',num2str(K),'_meanCorrmat_smoothed'), '-djpeg')
end

%% Examine Gaussian distribution means
figure;
for idx = 1:n_of_files
    hmm = results{idx, 'hmm_list'}{1};
    permutation = results{idx, 'permutation_list'};
    n_chan = length(results{idx, 'chanlocs_list'}{1});
    
    mean_arr = zeros(K, n_chan);
    
    for i = 1:K
        state = permutation(i);
        mean_arr(i,:) = getMean(hmm, state);
    end
    
    subplot(6, 5, idx);
    imagesc(mean_arr);
end

%% Examine fractional occupancy
for idx = 9    
    Gamma = results{idx, 'Gamma_list'}{1};
    T = results{idx, 'epoch_length_list'}{1};
    hmm = results{idx, 'hmm_list'}{1};
    
    fracOccup = getFractionalOccupancy(Gamma, T, hmm.train) % simply Gamma mean
    switchingRate = getSwitchingRate(Gamma, T, hmm.train) % total sum of diff(Gamma)
    tic
    stateLifeTimes = getStateLifeTimes(Gamma, T, hmm.train)
    stateIntervalTimes = getStateIntervalTimes(Gamma, T, hmm.train) 
    toc
    maxFracOccup = getMaxFractionalOccupancy(Gamma, T, hmm.train)
%     figure;
%     imagesc(fracOccup')    
    
end

%% Examine 













        
        
        
        
        
        
        
        
        