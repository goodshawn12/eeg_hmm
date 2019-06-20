%%
% addpath(genpath('.')) % assuming we are in the HMM-MAR directory

%%
K = 4; % number of states
ndim = 3; % number of channels
N = 10; % number of trials
Fs = 200;
T = 10000 * ones(N,1); % number of data points

hmmtrue = struct();
hmmtrue.K = K;
hmmtrue.state = struct();
hmmtrue.train.covtype = 'full';
hmmtrue.train.zeromean = 0;
hmmtrue.train.order = 0;

for k = 1:K
    hmmtrue.state(k).W.Mu_W = rand(1,ndim);
    r = randn(ndim);
    hmmtrue.state(k).Omega.Gam_rate = 0.1 * 1000 * (r') * r + eye(ndim);
    hmmtrue.state(k).Omega.Gam_shape = 1000;
end

hmmtrue.P = rand(K) + 100 * eye(K);  
for j=1:K
    hmmtrue.P(j,:) = hmmtrue.P(j,:) ./ sum(hmmtrue.P(j,:));
end
hmmtrue.Pi = ones(1,K); %rand(1,K);
hmmtrue.Pi = hmmtrue.Pi./sum(hmmtrue.Pi);

%% Generate some simulation data
[X,T,Gammatrue] = simhmmmar(T,hmmtrue,[]);

Flips = binornd(1,0.2,N,3);
unflipped_X = X; 
for n = 1:N
    if mean(Flips(n,:))>0.5 % we force to be more unflipped than flipped
        Flips(n,:) = 1 - Flips(n,:);
    end
    for d = 1:ndim
        ind = (1:T(n)) + sum(T(1:n-1));
        if Flips(n,d) 
            X(ind,d) = -X(ind,d);
        end
    end
end

options_sf = struct(); 
options_sf.maxlag = 10; 
options_sf.noruns = 20; 
options_sf.nbatch = 3;
options_sf.verbose = 1;
flips = findflip(X,T,options_sf);
X = flipdata(X,T,flips);

%% training a hmm model
options = struct();
options.K = K; 
options.Fs = Fs; 
options.covtype = 'full';
options.order = 0;
options.DirichletDiag = 2; 
options.zeromean = 0;
options.verbose = 1;
options.useParallel = 0;

[hmm, Gamma, Xi, vpath] = hmmmar(X,T,options);
%% plotting results
figure; subplot(3,1,1)
plot(Gammatrue(1:1000,:)), set(gca,'Title',text('String','True state path'))
set(gca,'ylim',[-0.2 1.2]); ylabel('state #')

subplot(3,1,2)
plot(Gamma(1:1000,:)), set(gca,'Title',text('String','Estimated state path'))
set(gca,'ylim',[-0.2 1.2]); ylabel('state #')

subplot(3,1,3)
plot(vpath(1:1000)), set(gca,'Title',text('String','Viterbi path'))
set(gca,'ylim',[0 hmm.K+1]); ylabel('state #')

figure
subplot(2,4,1), imagesc(getFuncConn(hmmtrue,1)), colormap('gray'), set(gca,'Title',text('String','Simulated covariance'))
subplot(2,4,2), imagesc(getFuncConn(hmmtrue,2)), colormap('gray'), set(gca,'Title',text('String','Simulated covariance'))
subplot(2,4,3), imagesc(getFuncConn(hmmtrue,3)), colormap('gray'), set(gca,'Title',text('String','Simulated covariance'))
subplot(2,4,4), imagesc(getFuncConn(hmmtrue,4)), colormap('gray'), set(gca,'Title',text('String','Simulated covariance'))

subplot(2,4,5), imagesc(getFuncConn(hmm,1)), colormap('gray'), set(gca,'Title',text('String','Inferred covariance'))
subplot(2,4,6), imagesc(getFuncConn(hmm,2)), colormap('gray'), set(gca,'Title',text('String','Inferred covariance'))
subplot(2,4,7), imagesc(getFuncConn(hmm,3)), colormap('gray'), set(gca,'Title',text('String','Inferred covariance'))
subplot(2,4,8), imagesc(getFuncConn(hmm,4)), colormap('gray'), set(gca,'Title',text('String','Inferred covariance'))

%% training a hmmmar model

options = struct();
options.K = K; 
options.Fs = Fs; 
options.covtype = 'full';
options.order = 2;
options.DirichletDiag = 2; 
options.zeromean = 1;
options.verbose = 1;
options.useParallel = 0;

[hmm, Gamma,Xi] = hmmmar(X,T,options);

%% cv
[mcv,cv] = cvhmmmar(X,T,options);

%% reculculate free energy
fe = hmmfe(X,T,hmm,Gamma,Xi);
sum(fe)

%% getting the MAR spectra
options.Fs = Fs; 
options.completelags = 1;
options.MLestimation = 1; 
options.order = 20; % increase the order
options.p = 0.01;
spectral_info = hmmspectramar(X,T,[],Gamma,options);

figure
for k=1:2
    subplot(1,2,k)
    plot(spectral_info.state(k).f,spectral_info.state(k).psd(:,1,1),'k', 'LineWidth',3)
end

%% getting the non-parametric spectra:
options_mt = struct('Fs',Fs);
options_mt.fpass = [1 48];
options_mt.tapers = [4 7];
options_mt.p = 0;
options_mt.win = 500;
options_mt.Gamma = Gamma;
spectral_info = hmmspectramt(X,T,Gamma,options_mt);

FO = getFractionalOccupancy (Gamma,T); % state fractional occupancies per session
LifeTimes = getStateLifeTimes (Gamma,T); % state life times
Intervals = getStateIntervalTimes (Gamma,T); % interval times between state visits
SwitchingRate =  getSwitchingRate(Gamma,T); % rate of switching between stats