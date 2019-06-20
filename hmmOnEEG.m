eeglab;
eegdata = pop_loadset();

%% Prepare data
timepoints = eegdata.data;

n_epochs = size(timepoints, 3);
Fs = eegdata.srate; 
epoch_length = 6 * Fs;

timepoints = reshape(timepoints, 32, []);
timepoints = timepoints';
size(timepoints)
%% Examine data
% for i = 1:size(timepoints, 2)
%     subplot(5,7,i);
%     plot(timepoints(:,i))
% end

%% PCA run
% [pca_coeff, pca_score, pca_latent] = pca(timepoints);

% variance_percentage = cumsum(pca_latent)/sum(pca_latent);
% plot(variance_percentage)
 
% pca_timepoints = pca_score * pca_coeff(:,1:18);
% 
% for i = 1:18
%     subplot(3,6,i);
%     plot(pca_timepoints(:,i))
% end
    
%% 
select_epochs = floor(n_epochs);
select_time = 1:(select_epochs * epoch_length);

% Regarding T, it can be either a (N X 1) vector (where N is the total number of trials 
% or segments for all subjects) containing the length of each trial/segment/subject, or 
% a (no. of subjects X 1) cell with each element containing a vector (no. of trials X 1) 
% reflecting the length of the trials for that particular subject.
X = timepoints(select_time,:);
T = epoch_length * ones(select_epochs, 1);
K = 3; 
use_stochastic = 1;
options = struct();

options.K = K; % number of states 
options.onpower = 0;
options.standardise = 0;
options.Fs = Fs; 
options.verbose = 1;
options.useParallel = 0;

% For TDE: order = 0, zeromean = 1, covtype = 'full'
options.embeddedlags = -7:7;
options.order = 0; % no autoregressive components
options.zeromean = 1; % model the mean
options.covtype = 'full'; % full covariance matrix

% options.cyc = 10;
% options.initcyc = 10;
% options.initrep = 1;

if use_stochastic
    options.BIGNinitbatch = 15;
    options.BIGNbatch = 15;
    options.BIGtol = 1e-5;
    options.BIGcyc = 500;
    options.BIGundertol_tostop = 5;
    options.BIGdelay = 5;
    options.BIGforgetrate = 0.9;
    options.BIGbase_weights = 0.9;
end

options

%% 
[hmm, Gamma, Xi, vpath] = hmmmar(X,T,options);
% [hmm, Gamma,~,~,~,~,fehist] = hmmmar(data,T,options);

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
zero_padded_vpath = [zeros(size(X,1) - size(vpath,1), 1); vpath];
state_by_epochs = reshape(zero_padded_vpath, [], 1500);
imagesc(state_by_epochs);

saveas(gcf, 'TDE001_epochs.jpg')