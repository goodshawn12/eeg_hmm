cd '/home/ting/Documents/eeg_hmm';
addpath('/home/ting/Documents/eeg_hmm');
addpath(genpath('/home/ting/Documents/eeg_hmm/Utility'))
addpath('/data/projects/Shawn/2019_HMM/data/');
eeglab;

%% Load all trained data

K_set = 2:5;
n_of_states = length(K_set);
n_of_files = 0;
results_list = cell(n_of_states, 1);

for state_idx = 1:n_of_states
    resultsDir = sprintf('/home/ting/Documents/eeg_hmm/AMICA_results/results_K%d/', K_set(state_idx));
    
    % Manual setup files to process
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

    % Prepare raw data
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

    % Process Gamma and vpath
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
    
    results_list{state_idx} = results;
end

all_results = table(K_set', results_list);
all_results.Properties.VariableNames{1} = 'K';
meta_info = table(filename_list);

%% Cross session data process
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

meta_info = deleteVars(meta_info, {'events_list', 'chanlocs_list', 'Fs_list'});
meta_info = addvars(meta_info, events_list, chanlocs_list, Fs_list);

%% Generate rt, rt_off
rt_list = cell(n_of_files, 1);
rt_off_list = cell(n_of_files, 1);
rt_latency_list = cell(n_of_files, 1);
rt_speed_list = cell(n_of_files, 1);

for idx = 1:n_of_files
    Fs = meta_info{idx, 'Fs_list'}{1};
    events = meta_info{idx, 'events_list'}{1};  
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
meta_info = deleteVars(meta_info, {'rt_list', 'rt_off_list', 'rt_latency_list', 'rt_speed_list'});
meta_info = addvars(meta_info, rt_list, rt_off_list, rt_latency_list, rt_speed_list);

% RT and 1/RT Outlier Removal
rt_clean_list = cell(n_of_files, 1);
rt_off_clean_list = cell(n_of_files, 1);
rt_latency_clean_list = cell(n_of_files, 1);
rt_speed_clean_list = cell(n_of_files, 1); 

for idx = 1:n_of_files
    Fs = meta_info{idx, 'Fs_list'}{1};
    rt = meta_info{idx, 'rt_list'}{1};
    rt_off = meta_info{idx, 'rt_off_list'}{1};
    rt_latency = meta_info{idx, 'rt_latency_list'}{1}; 
    rt_speed = meta_info{idx, 'rt_speed_list'}{1};
    
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

meta_info = deleteVars(meta_info, {'rt_clean_list', 'rt_off_clean_list', 'rt_latency_clean_list', 'rt_speed_clean_list'});
meta_info = addvars(meta_info, rt_clean_list, rt_off_clean_list, rt_latency_clean_list, rt_speed_clean_list);

%% Calculate state correlation with ERP
win_len_sec = 2;

for state_idx = 1:n_of_states
    K = K_set(state_idx);
    state_r_list = zeros(n_of_files, K);
    trial_Gamma_mean_list = cell(n_of_files, 1);   
    results = all_results{state_idx, 'results_list'}{1};
    
    for idx = 1:n_of_files
        Fs = meta_info{idx, 'Fs_list'}{1};
        Gamma = results{idx, 'Gamma_list'}{1};
        rt_speed = meta_info{idx, 'rt_speed_clean_list'}{1};
        rt_latency = meta_info{idx, 'rt_latency_clean_list'}{1};

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

    all_results{state_idx, 'results_list'}{1} = results;
end

%%
% cd '~/Documents/eeg_hmm/'
% writetable(compare, 'cross_AMICA_statistics.dat')

%%
corr_mean_list = cell(n_of_states, 1);
corr_std_list = cell(n_of_states, 1);

for state_idx = 1:n_of_states
    results = all_results{state_idx, 'results_list'}{1};
    
    sorted_corr = sort(results.state_r_list, 2);
    corr_mean_list{state_idx} = mean(sorted_corr, 1);
    corr_std_list{state_idx} = std(sorted_corr, 1);
end

all_results = deleteVars(all_results, {'corr_mean_list', 'corr_std_list'});
all_results = addvars(all_results, corr_mean_list, corr_std_list);

%% paired ttest
pairwise_count = sum(1:length(K_set)-1);
p_list = repmat(-1, 1, pairwise_count);
neg_p_list = repmat(-1, 1, pairwise_count);

for firstGroup = 1:length(K_set)-1
    for secondGroup = firstGroup+1:length(K_set)
        firstResult = all_results{firstGroup, 'results_list'}{1};
        firstList = sort(firstResult.state_r_list, 2);
        secondResult = all_results{secondGroup, 'results_list'}{1};
        secondList = sort(secondResult.state_r_list, 2);
        
        [~,p] = ttest(firstList(:,end), secondList(:,end));
        [~,neg_p] = ttest(firstList(:,1), secondList(:,1));
        index = find(p_list== -1, 1);
        p_list(index) = p;
        neg_p_list(index) = neg_p;       
    end
end

%% Mean Correlation Trend vs K
r = [1, 0.1, 0.1];
y = [1, 1, 0.1];
yg = [0.8, 1, 0.1];
g = [0.1, 1, 0.1];
b = [0.1, 0.1, 1];

mean_list = all_results.corr_mean_list;
std_list = all_results.corr_std_list;

yData = [];
errData = [];
for state_idx = 1:n_of_states
    yData = [yData mean_list{state_idx}];
    errData = [errData std_list{state_idx}];
    if state_idx ~= n_of_states
        yData = [yData 0];
        errData = [errData 0];
    end
end
yData = yData';
errData = errData';

figure;
baraxis = bar(1:length(yData), yData);
baraxis.FaceColor = 'flat';
baraxis.CData = [r;b; r; r;g;b; r; r;y;g;b; r; r;y;yg;g;b];
set(gca, 'XTick', [1.5, 5, 9.5, 15])
set(gca, 'XTickLabel', {'2', '3', '4', '5'})

hold on
er = errorbar(1:length(yData), yData, errData, errData);
er.LineStyle = 'none';
er.Color = 'k';

hold on

%2 6 11 17
% rej_list = p_list > 0.05;
% position_list = {[2,6],[2,11],[2,17],[6,11],[6,17],[11,17]};
% sigstar(position_list(~rej_list), p_list(~rej_list));
%1 4 8 13
rej_list = neg_p_list > 0.05;
position_list = {[1,4],[1,8],[1,13],[4,8],[4,13],[8,13]};
sigstar(position_list(~rej_list), neg_p_list(~rej_list));

title('Correlation Trends against Number of Inferred Model')
ylabel('Correlation with RS')
xlabel('Number of Inferred Models learnt by AMICA')


%% K2 - K4 MAR
r = [1, 0.1, 0.1];
y = [1, 1, 0.1];
yg = [0.8, 1, 0.1];
g = [0.1, 1, 0.1];
b = [0.1, 0.1, 1];

yData = [mean_MAR_list(1:2); 0; mean_MAR_list(3:5); 0; mean_MAR_list(6:9)];
errData = [std_MAR_list(1:2); 0; std_MAR_list(3:5); 0; std_MAR_list(6:9)];

figure;
baraxis = bar(1:length(yData), yData);
baraxis.FaceColor = 'flat';
baraxis.CData = [r;b; r; r;g;b; r; r;y;g;b;];
set(gca, 'XTick', [1.5,5,9.5])
set(gca, 'XTickLabel', {'2', '3', '4'})

hold on
er = errorbar(1:length(yData), yData, errData, errData);
er.LineStyle = 'none';
er.Color = 'k';

hold on
%2 6 11 17
sigstar({[2,6],[2,11],[6,11]},p_MAR_list);
%1 4 8 13
% sigstar({[1,4],[1,8],[4,8]},neg_p_MAR_list);

title('Correlation Trends against Number of Inferred Model')
ylabel('Correlation with RS')
xlabel('Number of Inferred Models learnt by HMM-MAR')

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