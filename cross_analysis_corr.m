addpath('/home/ting/Documents/eeglab');
addpath('/home/ting/Documents/eeg_hmm');
addpath(genpath('/home/ting/Documents/eeg_hmm/Utility'));
addpath(genpath('/home/ting/Documents/eeg_hmm/HMM-MAR'));
addpath('/data/projects/Shawn/2016 JNE/dataset/');
addpath('/data/projects/Shawn/2019_HMM/data/');
%%
graph_dir = '/home/ting/Documents/eeg_hmm/cross_session_graphs/';

%% Load all correlation
all_corr = csvread('/home/ting/Documents/eeg_hmm/cross_session_results/3_method_all_corr.csv', 1, 1);
all_corr = all_corr(1:54,:);
methods = {'HMM-GAU', 'AMICA'};
n_methods = length(methods);
Ks = 2:5;
n_Ks = length(Ks);
n_sessions = size(all_corr, 1)/length(methods);
drowsy_col_index = cumsum(1:n_Ks); alert_col_index = cumsum(2:(n_Ks+1));
drowsy_all = all_corr(:, drowsy_col_index); 
alert_all = all_corr(:, alert_col_index);

%% Test normality
p_drowsy_list = zeros(n_methods, n_Ks);
p_drowsy_method_list = zeros(1, n_methods);
p_drowsy_K_list = zeros(1, n_Ks);
p_alert_list = zeros(n_methods, n_Ks);
p_alert_method_list = zeros(1, n_methods);
p_alert_K_list = zeros(1, n_Ks);

for row = 1:n_methods
    row_selection = (1:n_sessions) + (row-1)*n_sessions;
    
    [~, p_drowsy_method] = swtest(reshape(drowsy_all(row_selection, :), 1, []));
    p_drowsy_method_list(row) = p_drowsy_method;
    [~, p_alert_method] = swtest(reshape(alert_all(row_selection, :), 1, []));
    p_alert_method_list(row) = p_alert_method;
    
    for col = 1:n_Ks
        if row == 1
            [~, p_drowsy_K] = swtest(drowsy_all(:, col));
            p_drowsy_K_list(col) = p_drowsy_K;
            [~, p_alert_K] = swtest(alert_all(:, col));
            p_alert_K_list(col) = p_alert_K;
        end
        
        [~, p_drowsy] = swtest(drowsy_all(row_selection, col));
        p_drowsy_list(row, col) = p_drowsy;
        [~, p_alert] = swtest(alert_all(row_selection, col));  
        p_alert_list(row, col) = p_alert;
    end
end

p_drowsy_list
p_drowsy_method_list
p_drowsy_K_list

p_alert_list
p_alert_method_list
p_alert_K_list

%% paired two-way ANOVA
drowsy_cell = cell(n_methods, n_Ks);
alert_cell = cell(n_methods, n_Ks);
for row = 1:n_methods
    row_selection = (1:n_sessions) + (row-1) * n_sessions;
    size(row_selection)
    for col = 1:n_Ks        
        drowsy_cell{row, col} = reshape(drowsy_all(row_selection, col), 1, []);
        alert_cell{row, col} = reshape(alert_all(row_selection, col), 1, []);
    end
end

[stats, df, pvals, surrog] = statcond(drowsy_cell, 'mode', 'bootstrap')
[stats, df, pvals, surrog] = statcond(alert_cell, 'mode', 'bootstrap')

%% Wilcoxon test on treatment K within the same method
[x_coor, y_coor] = meshgrid(1:length(Ks), 1:length(Ks));
x_coor = tril(x_coor, -1);
y_coor = tril(y_coor, -1);
K_pairs = [x_coor(find(x_coor)), y_coor(find(y_coor))];

[p,~,~] = mult_pair_signrank(drowsy_all, n_sessions, 1);
drowsy_K_signrank_p = [K_pairs, p]

[p,~,~] = mult_pair_signrank(alert_all, n_sessions, 1);
alert_K_signrank_p = [K_pairs, p]


%% Wilcoxon test on treatment method within the same K
[x_coor, y_coor] = meshgrid(1:length(methods), 1:length(methods));
x_coor = tril(x_coor, -1);
y_coor = tril(y_coor, -1);
method_pairs = [x_coor(find(x_coor)), y_coor(find(y_coor))];

[p,~,~] = mult_pair_signrank(drowsy_all, n_sessions, 2);
drowsy_method_signrank_p = [method_pairs, p]

[p,~,~] = mult_pair_signrank(alert_all, n_sessions, 2);
alert_method_signrank_p = [method_pairs, p]

%% Wilcoxon test on treatment K for all methods
[x_coor, y_coor] = meshgrid(1:length(Ks), 1:length(Ks));
x_coor = tril(x_coor, -1);
y_coor = tril(y_coor, -1);s
K_pairs = [x_coor(find(x_coor)), y_coor(find(y_coor))];

[p,~,~] = mult_pair_signrank(drowsy_all, length(methods)*n_sessions, 1);
drowsy_K_all_methods_signrank_p = [K_pairs, p]

[p,~,~] = mult_pair_signrank(alert_all, length(methods)*n_sessions, 1);
alert_K_all_methods_signrank_p = [K_pairs, p]

%% Wilcoxon test on treatment method for all Ks
drowsy_arranged = zeros(length(Ks)*length(methods)*n_sessions,1);
alert_arranged = zeros(length(Ks)*length(methods)*n_sessions,1);
rept = length(Ks)*n_sessions;
for i = 1:length(methods)
    drowsy_arranged((i-1)*rept+(1:rept),1) = reshape(drowsy_all((i-1)*n_sessions+(1:n_sessions),:), [], 1);
    alert_arranged((i-1)*rept+(1:rept),1) = reshape(alert_all((i-1)*n_sessions+(1:n_sessions),:), [], 1);
end

[x_coor, y_coor] = meshgrid(1:length(methods), 1:length(methods));
x_coor = tril(x_coor, -1);
y_coor = tril(y_coor, -1);
K_pairs = [x_coor(find(x_coor)), y_coor(find(y_coor))];

[p,~,~] = mult_pair_signrank(drowsy_arranged, rept, 2)
drowsy_method_all_K_signrank_p = [K_pairs, p]

[p,~,~] = mult_pair_signrank(alert_arranged, rept, 2)
alert_method_all_K_signrank_p = [K_pairs, p]


%% Two-way ANOVA all 3 methods
toSave = 1;

[~, tbl, drowsy_stats] = anova2(drowsy_all, n_sessions);
tbl
if toSave
    output_filename = strcat(graph_dir, 'all_methods_drowsy_anova');
    set(gcf, 'PaperPositionMode', 'auto');
    print(output_filename, '-djpeg')
end

[~, tbl, alert_stats] = anova2(alert_all, n_sessions);
tbl
if toSave
    output_filename = strcat(graph_dir, 'all_methods_alert_anova');
    set(gcf, 'PaperPositionMode', 'auto');
    print(output_filename, '-djpeg')
end

%% Two-way ANOVA excluding HMM-AR
toSave = 1;

[~, tbl, drowsy_stats] = anova2(drowsy_all(1:2*n_sessions,:), n_sessions);
if toSave
    output_filename = strcat(graph_dir, 'AMICA_GAU_drowsy_anova');
    set(gcf, 'PaperPositionMode', 'auto');
    print(output_filename, '-djpeg')
end

[~, tbl, alert_stats] = anova2(alert_all(1:2*n_sessions,:), n_sessions);
if toSave
    output_filename = strcat(graph_dir, 'AMICA_GAU_alert_anova');
    set(gcf, 'PaperPositionMode', 'auto');
    print(output_filename, '-djpeg')
end

%% Load all correlation overall statistics
corr_stats = csvread('/home/ting/Documents/eeg_hmm/cross_session_results/3_method_summary.csv', 1, 1);

methods = {'AMICA', 'HMM-GAU'};
n_methods = length(methods);
Ks = 2:5;
n_Ks = length(Ks);
corr_stats(3:3:9,:) = [];
drowsy_stats = corr_stats(:, (2:2:2*n_Ks) - 1);
alert_stats = corr_stats(:, 2:2:2*n_Ks);

%% Plot line graph
mean_list = drowsy_stats(1:n_methods,:);
sem_list = drowsy_stats(2*n_methods+(1:n_methods),:);

figure;
errorbar(Ks, -mean_list(1, :), sem_list(1, :), 'LineWidth', 2);
hold on;
errorbar(Ks, -mean_list(2, :), sem_list(2, :), 'LineWidth', 2);
hold on;
ylim([0.1 0.35])
xlabel('Number of states');
ylabel('Mean correlation');

plot(2:K, [0.15 0.15 0.15 0.15], 'k*', 'LineStyle', 'none'); hold on;
plot(2:3, [0.14 0.14], 'k*', 'LineStyle', 'none'); hold on;
sigstar({[1 2], [2 3], [2, 3], [2, 4]}, [0.05 0.01 0.06 0.05])
legend(methods, 'Location', 'northeast')

%% Plot all methods drowsy and alert mean and error;
toSave = 1;

r = [1, 0.2, 0.2];
g = [0.2, 1, 0.2];
b = [0.2, 0.2, 1];

mean_list = corr_stats(1:n_methods,:);
sem_list = corr_stats(2*n_methods+(1:n_methods),:);

xData = repmat(Ks', 1, 2*n_methods);
yData = transpose(reshape(mean_list, 2*n_methods, []));
errData = transpose(reshape(sem_list, 2*n_methods, []));

figure;
baraxis = bar(xData, yData);
colorData = [r;g;b];
for i = 1:2*n_methods
    mod_val = mod(i, n_methods);
    if mod_val == 0
        mod_val = n_methods;
    end
    baraxis(i).FaceColor = colorData(mod_val,:);
end

hold on
bar_x_endpoints = [];
for i = 1:size(xData,2)
    bar_x_endpoints = [bar_x_endpoints; baraxis(i).XEndPoints];
end
er = errorbar(bar_x_endpoints', yData, errData, 'LineStyle', 'none', 'Color', 'k');

legend([baraxis(1:n_methods), er(1)], [methods, 'SEM'], 'Location', 'northwest');
xlabel('Number of states');
ylabel('Mean correlation');

if toSave
    output_filename = strcat(graph_dir, 'AMICA_GAU_K_alert_drowsy_correlation');
    set(gcf, 'PaperPositionMode', 'auto');
%     saveas(cgfighandle, strcat(output_filename,'.fig'))
    print(output_filename, '-djpeg')
end
%%
hold on

if 1%sigstarAlert   
    %2 6 11 17
    rej_list = p_list > 0.05;
    position_list = {[2,6],[2,11],[2,17],[6,11],[6,17],[11,17]};
    sigstar(position_list(~rej_list), p_list(~rej_list));
else
    %1 4 8 13
    rej_list = neg_p_list > 0.05;
    position_list = {[1,4],[1,8],[1,13],[4,8],[4,13],[8,13]};
    sigstar(position_list(~rej_list), neg_p_list(~rej_list));
end

title('Correlation Trends against Number of Inferred Model')
ylabel('Correlation with RS')
xlabel('Number of Inferred Models of HMM-MAR')
toSave = 0;
if toSave
    if sigstarAlert
        output_filename = strcat(method, '_Correlation_K_trend_', 'AlertSigTest_trend');
    else
        output_filename = strcat(method, '_Correlation_K_trend_', 'DrowsySigTest_trend');
    end
    set(gcf, 'PaperPositionMode', 'auto');
%     saveas(cgfighandle, strcat(output_filename,'.fig'))
    print(output_filename, '-djpeg')
end

