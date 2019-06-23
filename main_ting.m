cd '/home/ting/Documents/eeg_hmm';
addpath('/home/ting/Documents/eeglab')
eeglab;
addpath(genpath('HMM-MAR'))
addpath('/data/projects/Shawn/2019_HMM/data')
filename = 'session53_ASR20.set';
eegdata = pop_loadset(filename);

%% Prepare data
Fs = eegdata.srate;
epoch_length = eegdata.pnts;
n_epochs = 1;
if ~isempty(eegdata.epoch)
    n_epochs = size(eegdata.epoch, 2);
end
timepoints = eegdata.data;
    
%% training data setup
X = [];
T = [];
select_start = 0.25;
select_end = 0.75;

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

%% training parameters setup
% Regarding T, it can be either a (N X 1) vector (where N is the total number of trials 
% or segments for all subjects) containing the length of each trial/segment/subject, or 
% a (no. of subjects X 1) cell with each element containing a vector (no. of trials X 1) 
% reflecting the length of the trials for that particular subject.
K = 3; 
use_stochastic = 1;
method = 'MAR';

options = struct();
options.K = K; % number of states 
options.Fs = Fs; 
options.verbose = 1;
options.useParallel = 1;
options.standardise = 0;
options.onpower = 0;

if strcmp(method,'MAR')
    options.order = 3;
    options.zeromean = 1;
    options.covtype = 'diag';
    options.DirichletDiag = 250;
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
end

if use_stochastic && n_epochs > 1
    options.BIGNinitbatch = 15;
    options.BIGNbatch = 15;
    options.BIGtol = 1e-5;
    options.BIGcyc = 100;
    options.BIGundertol_tostop = 5;
    %options.BIGdelay = 5;
    options.BIGforgetrate = 0.7;
    options.BIGbase_weights = 0.9;
end

%% 
options
tic
[hmm, Gamma, Xi, vpath] = hmmmar(X,T,options);
% [hmm, Gamma,~,~,~,~,fehist] = hmmmar(data,T,options);
toc
%% plotting results
figure;
subplot(2,1,1)
plot(Gamma(:,:)), set(gca,'Title',text('String','Estimated state path'))
set(gca,'ylim',[-0.2 1.2]); ylabel('state #')

subplot(2,1,2)
plot(vpath(:)), set(gca,'Title',text('String','Viterbi path'))
set(gca,'ylim',[0 hmm.K+1]); ylabel('state #')

% saveas(gcf, 'TDE001_path.jpg')
%% Colormap
figure;
state_by_epochs = (reshape(vpath, [], size(Gamma, 1)/n_epochs));
imagesc(state_by_epochs);
set(gca,'Title',text('String','State arranged by epochs'))
colormap(parula)

% saveas(gcf, 'TDE001_epochs.jpg')