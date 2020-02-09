cd '/home/ting/Documents/eeg_hmm';
addpath('/home/ting/Documents/eeglab')
eeglab;
addpath(genpath('HMM-MAR'))
addpath('/data/projects/Shawn/2019_HMM/data')
filename = 'session53_ASR20.set';
eegdata = pop_loadset(filename);

%/data/projects/Shawn/2016 JNE/dataset

%% Prepare raw data
Fs = eegdata.srate;
epoch_length = eegdata.pnts;
n_epochs = 1;
if ~isempty(eegdata.epoch)
    n_epochs = size(eegdata.epoch, 2);
end
timepoints = eegdata.data;
    
%% Prepare training data
X = [];
T = [];
select_start = 0;
select_end = 1;

if n_epochs > 1
    start_epoch = floor(n_epochs * select_start);
    if start_epoch == 0
        start_epoch = 1;
    end
    end_epoch = floor(n_epochs * select_end); 
    X = timepoints(:,:,start_epoch:end_epoch); 
    X = reshape(X, eegdata.nbchan, []);
    X = X';
    
    select_epochs = end_epoch - start_epoch + 1;
    T = epoch_length * ones(select_epochs, 1);
else
    start_timepoint = floor(epoch_length * select_start);
    if start_timepoint == 0
        start_timepoint = 1;
    end
    end_timepoint = floor(epoch_length * select_end);    
    X = timepoints(:,start_timepoint:end_timepoint);
    X = X';
    
    select_timepoints = end_timepoint - start_timepoint + 1;
    T = select_timepoints;
end

%% Prepare global training parameters
% Regarding T, it can be either a (N X 1) vector (where N is the total number of trials 
% or segments for all subjects) containing the length of each trial/segment/subject, or 
% a (no. of subjects X 1) cell with each element containing a vector (no. of trials X 1) 
% reflecting the length of the trials for that particular subject.
K = 3; 
use_stochastic = 0;
method = 'MIX';

cyc = 100;
tol = 1e-5;
initrep = 3;
initcycle = 50;

options = struct();
options.K = K; % number of states 
options.Fs = Fs; 
options.verbose = 1;
options.useParallel = 0;
options.standardise = 0;
options.onpower = 0;
options.DirichletDiag = 100;

options.cyc = cyc;
options.tol = tol;
options.initrep = initrep;
options.initcyc = initcycle;

if strcmp(method,'MAR')
    options.order = 1;
    options.zeromean = 1;
    options.covtype = 'diag';
elseif strcmp(method, 'TDE')
    % For TDE: order = 0, zeromean = 1, covtype = 'full'
    options.embeddedlags = -7:7;
    options.order = 0; % no autoregressive components
    options.zeromean = 0; % model the mean
    options.covtype = 'full'; % full covariance matrix
elseif strcmp(method, "GAUSSIAN")
    options.order = 0;
    options.zeromean = 0;
    options.covtype = 'full';     
    options.onpower = 1; 
elseif strcmp(method, "MIX")
    % default
end

if use_stochastic && n_epochs > 1
    options.BIGNinitbatch = 15;
    options.BIGNbatch = 15;
    options.BIGtol = tol;
    options.BIGcyc = cyc;
    options.BIGinitrep = initrep;
    options.BIGundertol_tostop = 5;
    options.BIGforgetrate = 0.7;
    options.BIGbase_weights = 0.9;
end

%% Prepare training parameters for parallel environment
options_header = {'K', 'order', 'covtype', 'zeromean', 'onpower', 'DirichletDiag'};

method_list =   {'MAR', 'MAR',  'MAR',  'GAU',  'GAU',  'GAU',  'GAU'};
K_list =        {3,     3,      3,      3,      3,      3,      4};
order_list =    {2,     1,      1,      0,      0,      0,      0};
covtype_list =  {'diag','full', 'diag', 'full', 'full', 'diag', 'full'};
zeromean_list = {1,     1,   	1,      0,      0,      0,      0};
onpower_list =  {0,     0,      0,      1,      0,      1,      1};
dirichlet_list ={100,   100,    10000,  100,    100,    100,    100};
options_cellarray = [K_list; order_list; covtype_list; zeromean_list; onpower_list; dirichlet_list];

options_selection = 1;
method_list = method_list(options_selection);
options_cellarray = options_cellarray(:, options_selection);

options_list = cell2struct(options_cellarray, options_header, 1);

%% Train
delete(gcp('nocreate')); % shut down any current pool
npar = size(options_cellarray, 2);
parpool(npar);   % request workers from the cluster

results_list = {};
results_list_header = {'hmm', 'Gamma', 'Xi', 'vpath', 'fehist', ...
    'time_elapsed', 'select_start', 'select_end', 'training_data_size'};

parfor par_id = 1:npar
    opts = mergeStructs(options, options_list(par_id))
    start_time = tic;
    [hmm, Gamma, Xi, vpath, ~, ~, fehist] = hmmmar(X,T,opts);
    time_elapsed = toc(start_time)

    % save output and meta data
    result = struct();
    result.hmm = hmm;
    result.Gamma = Gamma;
    result.Xi = Xi;
    result.vpath = vpath;
    result.fehist = fehist;
    result.time_elapsed = time_elapsed;
    result.select_start = select_start;
    result.select_end = select_end;
    result.training_data_size = size(X);
    
    results_list{par_id} = result;   
end

for i = 1:length(results_list)
    % Find the first non-duplicate filename
    fileindex = 0;
    result_filename_prefix = strcat(method_list{i},'_',extractBefore(filename,'.set'),'_');
    result_filename = strcat(result_filename_prefix,num2str(fileindex),'.mat');
    while isfile(strcat('results/',result_filename))
        fileindex = fileindex + 1;
        result_filename = strcat(result_filename_prefix,num2str(fileindex),'.mat');
    end
    
    result = results_list{i};
    save(strcat('results/', result_filename), '-struct', 'result');
end

%% Plotting paths
figure;
subplot(2,1,1)
plot(Gamma(:,:)), set(gca,'Title',text('String','Estimated state path'))
set(gca,'ylim',[-0.2 1.2]); ylabel('state #')

subplot(2,1,2)
plot(vpath(:)), set(gca,'Title',text('String','Viterbi path'))
set(gca,'ylim',[0 hmm.K+1]); ylabel('state #')

%% Colormap by epochs
state_by_epoch = [];
if n_epochs > 1
    figure;
    state_by_epoch = (reshape(vpath, [], size(Gamma, 1)/n_epochs));
    imagesc(state_by_epoch);
    set(gca,'Title',text('String','State arranged by epochs'))
    colormap(parula)
else
    zero_padded_vpath = [zeros(size(X,1)-size(vpath,1),1); vpath];
    
    all_events = eegdata.event;
    epoch_start_offset = -2 * Fs;
    epoch_end_offset = 4 * Fs;
    state_by_epoch = zeros(length(all_events), epoch_end_offset - epoch_start_offset);
    
    rt = [];
    rt_off = [];
    for i = 1:length(eegdata.event)
        event = all_events(i);
        sliced_epoch_start = epoch_start_offset + event.latency;
        sliced_epoch_end = epoch_end_offset + event.latency -1;
        if (strcmp(event.type, '251') || strcmp(event.type, '252')) && ...
            sliced_epoch_start >= start_timepoint && ...
            sliced_epoch_end <= end_timepoint
        
            rt = [rt, (all_events(i+1).latency-all_events(i).latency)];
            rt_off = [rt_off, (all_events(i+2).latency-all_events(i).latency)];
            sliced_epoch = zero_padded_vpath(sliced_epoch_start:sliced_epoch_end);
            state_by_epoch(i,:) = sliced_epoch';
        end
    end
    state_by_epoch(~any(state_by_epoch,2),:) = [];
end

[sortedRT,sortIdx] = sort(rt);
[sortedRT_off,sortIdx_off] = sort(rt_off);


figure;
imagesc(state_by_epoch(sortIdx,:));
set(gca,'Title',text('String','State arranged by epochs'))
cmap = [0.1,0.1,1; 0.1,1,0.1; 1,0.1,0.1; 1,1,0.1];
colormap(cmap(:,:))

hold on,
plot(sortedRT-epoch_start_offset,1:length(rt),'linewidth',2,'color','w')
line([-epoch_start_offset, -epoch_start_offset], [1, length(rt)],'linewidth',2,'color','k')
plot(sortedRT_off-epoch_start_offset,1:length(rt_off),'linewidth',2,'color',[0.7,0.7,0.7])

%% Colormap by epochs timelock 253
state_by_epoch = [];
if n_epochs > 1
    figure;
    state_by_epoch = (reshape(vpath, [], size(Gamma, 1)/n_epochs));
    imagesc(state_by_epoch);
    set(gca,'Title',text('String','State arranged by epochs'))
    colormap(parula)
else
    zero_padded_vpath = [zeros(size(X,1)-size(vpath,1),1); vpath];
    
    all_events = eegdata.event;
    epoch_start_offset = -2 * Fs;
    epoch_end_offset = 4 * Fs;
    state_by_epoch = zeros(length(all_events), epoch_end_offset - epoch_start_offset);
    
    rt = [];
    rt_off = [];
    for i = 1:length(eegdata.event)
        event = all_events(i);
        sliced_epoch_start = epoch_start_offset + event.latency;
        sliced_epoch_end = epoch_end_offset + event.latency -1;
        if (strcmp(event.type, '253') || strcmp(event.type, '253')) && ...
            sliced_epoch_start >= start_timepoint && ...
            sliced_epoch_end <= end_timepoint
        
            rt = [rt, (all_events(i-1).latency-all_events(i).latency)];
            rt_off = [rt_off, (all_events(i+1).latency-all_events(i).latency)];
            sliced_epoch = zero_padded_vpath(sliced_epoch_start:sliced_epoch_end);
            state_by_epoch(i,:) = sliced_epoch';
        end
    end
    state_by_epoch(~any(state_by_epoch,2),:) = [];
end

[sortedRT,sortIdx] = sort(rt);
[sortedRT_off,sortIdx_off] = sort(rt_off);


figure;
imagesc(state_by_epoch(sortIdx,:));
set(gca,'Title',text('String','State arranged by epochs'))
cmap = [0.1,0.1,1; 0.1,1,0.1; 1,0.1,0.1; 1,1,0.1];
colormap(cmap(1:4,:))

hold on,
plot(sortedRT-epoch_start_offset,1:length(rt),'linewidth',2,'color','w')
line([-epoch_start_offset, -epoch_start_offset], [1, length(rt)],'linewidth',2,'color','k')
plot(sortedRT_off-epoch_start_offset,1:length(rt_off),'linewidth',2,'color',[0.7,0.7,0.7])


%% Plot smoothed state time

winLen = 30 * Fs;
walkLen = 5 * Fs;
[nTime,nMod] = size(Gamma);

stateProb = zeros(nMod,floor(nTime/walkLen));
for it = 1:floor((nTime-winLen+walkLen)/walkLen)
    time_range = (it-1)*walkLen+1:(it-1)*walkLen+winLen;
    stateProb(:,it) = mean(Gamma(time_range,:),1)';
end

figure, imagesc(stateProb)
