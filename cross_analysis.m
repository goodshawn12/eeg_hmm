addpath('/home/ting/Documents/eeglab');
addpath('/home/ting/Documents/eeg_hmm');
addpath(genpath('/home/ting/Documents/eeg_hmm/Utility'));
addpath(genpath('/home/ting/Documents/eeg_hmm/HMM-MAR'));
addpath('/data/projects/Shawn/2016 JNE/dataset/');
addpath('/data/projects/Shawn/2019_HMM/data/');
%%
all_corr = csvread('/home/ting/Documents/eeg_hmm/cross_session_results/all_90s_smoothed_correlation.csv', 1, 1);

methods = {'AMICA', 'GAU', 'AR'};
Ks = 2:5;
n_sessions = size(all_corr, 1)/length(methods);
drowsy_col_index = cumsum(1:4); alert_col_index = cumsum(2:5);
drowsy_all = all_corr(:, drowsy_col_index); alert_all = all_corr(:, alert_col_index);

%% Wilcoxon test on treatment K within the same method
[x_coor, y_coor] = meshgrid(1:length(Ks), 1:length(Ks));
x_coor = tril(x_coor, -1);
y_coor = tril(y_coor, -1);
K_pairs = [x_coor(find(x_coor)), y_coor(find(y_coor))];

[p,~,~] = mult_pair_signrank(drowsy_all, n_sessions, 1);
drowsy_K_signrank_p = [K_pairs, p];

[p,~,~] = mult_pair_signrank(alert_all, n_sessions, 1);
alert_K_signrank_p = [K_pairs, p];


%% Wilcoxon test on treatment method within the same K
[x_coor, y_coor] = meshgrid(1:length(methods), 1:length(methods));
x_coor = tril(x_coor, -1);
y_coor = tril(y_coor, -1);
method_pairs = [x_coor(find(x_coor)), y_coor(find(y_coor))];

[p,~,~] = mult_pair_signrank(drowsy_all, n_sessions, 2);
drowsy_method_signrank_p = [method_pairs, p];

[p,~,~] = mult_pair_signrank(alert_all, n_sessions, 2);
alert_method_signrank_p = [method_pairs, p];


