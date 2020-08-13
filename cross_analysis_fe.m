addpath('/home/ting/Documents/eeglab');
addpath('/home/ting/Documents/eeg_hmm');
addpath(genpath('/home/ting/Documents/eeg_hmm/Utility'));
addpath(genpath('/home/ting/Documents/eeg_hmm/HMM-MAR'));
addpath('/data/projects/Shawn/2016 JNE/dataset/');
addpath('/data/projects/Shawn/2019_HMM/data/');
%%
graph_dir = '/home/ting/Documents/eeg_hmm/cross_session_graphs/';

%% Load all correlation
fe_list = readtable('/home/ting/Documents/eeg_hmm/cross_session_results/fe_trend_K.csv');
fe_list = transpose(table2cell(fe_list));
fe_list = cell2table(fe_list(2:end,:), 'VariableNames', fe_list(1,:));
n_of_files = size(fe_list, 2) - 1;

%% Plot session specific fe trend vs K
visible = 'on';
toSave = 1;

figure('Visible', visible)
for idx = 1:n_of_files
    subplot(5, 6, idx);
    plot(fe_list.K, fe_list{:, idx + 1});
    xticks(fe_list.K);
    xticklabels(fe_list.K);
end

if toSave
    set(gcf, 'PaperPositionMode', 'auto');
    output_filename = strcat(graph_dir, 'fe_trend_K');
    saveas(gcf, strcat(output_filename,'.fig'))
%     print(output_filename, '-djpeg')
end

