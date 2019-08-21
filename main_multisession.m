cd '/home/ting/Documents/eeg_hmm';
addpath('/home/ting/Documents/eeglab')
eeglab;
addpath(genpath('HMM-MAR'))
addpath('/data/projects/Shawn/2016 JNE/dataset')
%/data/projects/Shawn/2016 JNE/dataset

%% Indexing datafile
data_base_dir = '/data/projects/Shawn/2016 JNE/dataset/';
data_filelist = dir(strcat(data_base_dir, '*.set'));

data_filenames = {};
for i = 1:length(data_filelist)
    data_filenames{i} = data_filelist(i).name;
end

%% Prepare raw data
Fs = 250;
n_epoch = 1;
select_start = 0;
select_end = 1;

eegdata_list = {};
T_list = {};
for i = 1:length(data_filenames)
    filename = data_filenames{i};
    eegdata_list{i} = pop_loadset(filename);
    T_list{i} = eegdata_list{i}.pnts;
end

%% Remove session 53
% eegdata_list{8} = [];
% T_list{8} = [];
% eegdata_list = eegdata_list(~cellfun('isempty', eegdata_list));
% T_list = T_list(~cellfun('isempty', T_list));

%%
data_filenames{8} = [];
data_filenames = data_filenames(~cellfun('isempty', data_filenames));
%% Prepare global training parameters
% Regarding T, it can be either a (N X 1) vector (where N is the total number of trials 
% or segments for all subjects) containing the length of each trial/segment/subject, or 
% a (no. of subjects X 1) cell with each element containing a vector (no. of trials X 1) 
% reflecting the length of the trials for that particular subject.
K = 3; 
use_stochastic = 0;
method = 'TDE';

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
elseif strcmp(method, 'GAU')
    options.order = 0;
    options.zeromean = 0;
    options.covtype = 'full';     
    options.onpower = 1; 
elseif strcmp(method, 'MIX')
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

%% Train
delete(gcp('nocreate')); % shut down any current pool
npar = size(eegdata_list, 2);
parpool(npar);   % request workers from the cluster

results_list = {};
results_list_header = {'hmm', 'Gamma', 'Xi', 'vpath', 'fehist', ...
    'time_elapsed', 'select_start', 'select_end', 'training_data_size'};

parfor par_id = 1:npar
    start_time = tic;
    [hmm, Gamma, Xi, vpath, ~, ~, fehist] = ...
        hmmmar(transpose(eegdata_list{par_id}.data), T_list{par_id}, options);
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
    result.training_data_size = [eegdata_list{par_id}.pnts eegdata_list{par_id}.nbchan];
    
    results_list{par_id} = result;   
end

fprintf('Training done');

for i = 1:length(results_list)
    % Find the first non-duplicate filename
    fileindex = 0;
    result_filename_prefix = strcat(method,'_',extractBefore(data_filenames{i},'.set'),'_');
    result_filename = strcat(result_filename_prefix,num2str(fileindex),'.mat');
    while isfile(strcat('results/',result_filename))
        fileindex = fileindex + 1;
        result_filename = strcat(result_filename_prefix,num2str(fileindex),'.mat');
    end
    
    result = results_list{i};
    save(strcat('results/', result_filename), '-struct', 'result');
end

