addpath('/home/ting/Documents/eeglab');
addpath('/home/ting/Documents/eeg_hmm')
addpath(genpath('HMM-MAR'));
addpath('/data/projects/Shawn/2016 JNE/dataset/');
addpath('/data/projects/Shawn/2019_HMM/data/');

%% ### This pipeline only supports results with the same K (number of inferred states)
cd '/home/ting/Documents/eeg_hmm/results_K4';
K = 4;

%% Manual setup files to process
data_filelist = dir('MAR*.mat');

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

results = table(filename_list, session_name_list, K_list, epoch_length_list, timepoint_start_list, ...
    timepoint_end_list, Fs_list, hmm_list, Gamma_list, vpath_list);

%% Load events
events_list = cell(n_of_files, 1);

% For all sessions
for idx = 1:n_of_files
    filename = split(filename_list{idx}, '_');
    filename = join(filename(2:length(filename)-1), '_');
    filename = strcat(filename, '.set');
    
    load('-mat', filename{1});
    events_list{idx} = EEG.event;
end

results = deleteVars(results, 'events_list');
results = addvars(results, events_list);

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

%%
mean_list = zeros(2+3+4, 1);
std_list = zeros(2+3+4, 1);

mean_list(1:2) = mean(compare.Var3, 1);
std_list(1:2) = std(compare.Var3, 1);
mean_list(3:5) = mean(compare.Var2, 1);
std_list(3:5) = std(compare.Var2, 1);
mean_list(6:9) = mean(compare.Var1, 1);
std_list(6:9) = std(compare.Var1, 1);

%%
r = [1, 0.1, 0.1];
b = [0.1, 0.1, 1];
y = [1, 1, 0.1];
g = [0.1, 1, 0.1];

yData = [mean_list(1:2); 0; mean_list(3:5); 0; mean_list(6:9)];
errData = [std_list(1:2); 0; std_list(3:5); 0; std_list(6:9)]

figure;
baraxis = bar(1:11, yData);
baraxis.FaceColor = 'flat'
baraxis.CData = [r;b;r;r;g;b;r;r;y;g;b]
set(gca, 'XTick', [1.5,5,9.5])
set(gca, 'XTickLabel', {'K2', 'K3', 'K4'})
hold on
er = errorbar(1:11, yData, errData, errData);
er.LineStyle = 'none';

title('Correlation trends against number of inferred model')
ylabel('Correlation with RS')
xlabel('Models learnt by HMM-Gaussian')

%% ######## Graphing Sections ########
%########################################
% Preparation
cd 'graphs';

if K == 4
    cmap = [1,0.1,0.1; 1,1,0.1; 0.1,1,0.1; 0.1,0.1,1];
elseif K == 3
    cmap = [1,0.1,0.1; 0.1,1,0.1; 0.1,0.1,1];
elseif K == 2
    cmap = [1,0.1,0.1; 0.1,0.1,1];
else
    cmap = parula;
end

cmap_list = repmat(cmap, 1, 1, n_of_files);
for i = 1:n_of_files
    cmap_list(:,:,i) = cmap_list(results.inverse_permutation_list(i,:),:,i);
end
toSave = 0;
visible = 'off';

%%