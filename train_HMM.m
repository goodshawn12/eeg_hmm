cd '/home/ting/Documents/eeg_hmm';
addpath('/home/ting/Documents/eeglab')
addpath(genpath('HMM-MAR'))
addpath('/data/projects/Shawn/2019_HMM/data')
%/data/projects/Shawn/2016 JNE/dataset

%% Indexing datafile
data_base_dir = '/data/projects/Shawn/2019_HMM/data/';
% data_nase_dir = '/data/projects/Shawn/2016 JNE/dataset/';

data_filelist = dir(strcat(data_base_dir, '*.set'));
data_filenames = {};
output_filenames = {};
for i = 1:length(data_filelist)
    filename = data_filelist(i).name;
    output_filename = split(filename, '.');
    output_filename = output_filename{1};
      
    data_filenames{i} = filename;
    output_filenames{i} = output_filename;
end

data_filenames = data_filenames(~cellfun('isempty', data_filenames));
output_filenames = output_filenames(~cellfun('isempty', output_filenames));
n_of_files = length(data_filenames);
%% Prepare raw data
eegdata_list = cell(1, n_of_files);
eeglab;
for i = 1:n_of_files
    filename = data_filenames{i};
    eegdata_list{i} = pop_loadset(filename);
end

%% Prepare global training parameters
% Regarding T, it can be either a (N X 1) vector (where N is the total number of trials 
% or segments for all subjects) containing the length of each trial/segment/subject, or 
% a (no. of subjects X 1) cell with each element containing a vector (no. of trials X 1) 
% reflecting the length of the trials for that particular subject.
n_epoch = 1;
select_start = 0;
select_end = 1;
Fs = 250;

K = 5; 
use_stochastic = 0;
method = 'MAR';

% cyc = 1000;
% tol = 1e-5;
% initrep = 3;
% initcycle = 50;

options = struct();
options.K = K; % number of states
options.Fs = Fs;
options.verbose = 1;
options.useParallel = 0;
options.standardise = 0;
options.onpower = 0;
% options.DirichletDiag = 100;

% options.cyc = cyc;
% options.tol = tol;
% options.initrep = initrep;
% options.initcyc = initcycle;

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
npar = 27;
parpool(npar);   % request workers from the cluster

options_list = repmat(options, 1, n_of_files);
results_list = cell(1, n_of_files);
results_list_header = {'hmm', 'Gamma', 'Xi', 'vpath', 'fehist', ...
    'time_elapsed', 'select_start', 'select_end', 'training_data_size'};

parfor (idx = 1:n_of_files, npar)
    X = transpose(eegdata_list{idx}.data);
    T = eegdata_list{idx}.pnts;
    Fs = eegdata_list{idx}.srate;
    options_list(idx).Fs = Fs;
    
    start_time = tic;
    [hmm, Gamma, Xi, vpath, ~, ~, fehist] = hmmmar(X, T, options_list(idx));
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
    
    results_list{idx} = result;   
end

fprintf('Training done');

for idx = 1:n_of_files
    % Find the first non-duplicate filename
    fileindex = 0;
    result_filename_prefix = strcat(method, '_', output_filenames{idx});
    result_filename = strcat(result_filename_prefix, '_', num2str(fileindex), '.mat');
    while isfile(strcat('results_K5/',result_filename))
        fileindex = fileindex + 1;
        result_filename = strcat(result_filename_prefix, '_', num2str(fileindex), '.mat');
    end
    
    result = results_list{idx};
    save(strcat('results_K5/', result_filename), '-struct', 'result');
end

delete(gcp('nocreate'));
