cd '/home/ting/Documents/eeg_hmm';
addpath('/home/ting/Documents/eeg_hmm');
addpath(genpath('/home/ting/Documents/eeg_hmm/Utility'))
addpath(genpath('/home/ting/Documents/eeg_hmm/HMM-MAR'));
addpath('/data/projects/Shawn/2016 JNE/dataset/');
addpath('/data/projects/Shawn/2019_HMM/data/');

%% Specify number of inferred states
K = 3
method = 'GAU'
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
for i = 1:length(data_filelist)
    file = data_filelist(i);
    if contains(file.name, 'ASR20') || ~file.isdir || strcmp(file.name,'.') || strcmp(file.name,'..')
        continue;
    end
    filename_list{i} = file.name;
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
        HMM_out = load(filename_list{idx}, vars_to_load{:});
        training_data_size_list(idx,:) = HMM_out.training_data_size;
        hmm_list{idx} = HMM_out.hmm;
        Gamma_list{idx} = HMM_out.Gamma;
        vpath_list{idx} = HMM_out.vpath;
    end    
end

results = table(filename_list, training_data_size_list, Gamma_list, vpath_list);

if ~AMICA
    results = deleteVars(results, {'hmm_list'});
    results = addvars(results, hmm_list);
end

fprintf('Raw results data loaded.\n')

%% Load EEG related metadata
events_list = cell(n_of_files, 1);
chanlocs_list = cell(n_of_files, 1);
Fs_list = cell(n_of_files,1);

for idx = 1:n_of_files
    filename = split(filename_list{idx}, '_');
    filename = join(filename(2:length(filename)-1), '_');
    filename = strcat(filename, '.set');
    
    load('-mat', filename{1});
    events_list{idx} = EEG.event;
    chanlocs_list{idx} = extractfield(EEG.chanlocs, 'labels');
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
    Fs = metaInfo{idx, 'Fs_list'}{1};
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

fprintf('rt generated');

%% Smoothed state correlation with RS
smoothed_state_r_list = zeros(n_of_files, K);

win_len_sec = 5;
smoothing_range_sec = 90;

for idx = 1:n_of_files
    Fs = metaInfo{idx, 'Fs_list'}{1};
    Gamma = results{idx, 'Gamma_list'}{1};
    rs = metaInfo{idx, 'rt_speed_clean_list'}{1};
    rt_latency = metaInfo{idx, 'rt_latency_clean_list'}{1};

    state_r_list(idx, :) = round(getCorr(Gamma, K, Fs, rs, rt_latency, win_len_sec, smoothing_range_sec), 3);
end

% The sort is in ascending order so the first index in the permutation_list
% marks the lowest correlation (drowsy model) 
[~, permutation_list] = sort(state_r_list, 2);
inverse_permutation_list = zeros(size(permutation_list));
for row=1:size(permutation_list,1)
    inverse_permutation_list(row, permutation_list(row,:)) = 1:size(permutation_list,2);
end

results = deleteVars(results, {'state_r_list', 'permutation_list', 'inverse_permutation_list'});
results = addvars(results, state_r_list, permutation_list, inverse_permutation_list);

fprintf('correlation generated');

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

fprintf('fractional occupancy generated');

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

fprintf('transition probability matrix generated')

%%
