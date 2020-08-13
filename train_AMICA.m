cd '/home/ting/Documents/eeg_hmm';
addpath('/data/projects/Shawn/2019_HMM/data')
eeglab;
%% Indexing datafile
data_base_dir = '/data/projects/Shawn/2019_HMM/data/';

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

%% Train
K = 7;

for idx = 1:n_of_files
    filename = data_filenames{idx};
    EEG = pop_loadset(filename);
    output_filename = output_filenames{idx}
    outdir = sprintf('/home/ting/Documents/eeg_hmm/AMICA_results/results_K%d/%s/', K, output_filename);

    % run AMICA
    [weights,sphere,mods] = runamica15c(EEG.data, 'num_models', K, 'outdir', outdir, ...
                'nodeproc',[8 24],'do_reject',1,'numrej',15,'rejsig',3,'rejint',1);
    
end

fprintf('Training done');

%%
% eeglab
% cd /home/goodshawn12/MATLAB/2017_AMICA_sleep
% addpath(genpath('./'));
% addpath('/home/goodshawn12/MATLAB/Utility/plot/export_fig');
% addpath('/data/projects/Shawn/2017_Sleep')
% 
% %%
% for subjID = 1:10
%     subject = sprintf('nfle%d',subjID);
%     
%     if ~exist([subject '.set'])
%         EEG = pop_biosig(sprintf('/home/goodshawn12/MATLAB/AMICA/sleep/%s.edf',subject));
%         EEG.setname = subject;
%         EEG = pop_select( EEG);
%         EEG = eeg_checkset( EEG );
%         pop_saveset(EEG,'filename',subject,'filepath','/home/goodshawn12/MATLAB/AMICA/sleep/');
%         ScoringReader(sprintf('%s.txt',subject), '/home/goodshawn12/MATLAB/AMICA/sleep/');
%     else
%         EEG = pop_loadset([subject '.set']);
%     end
%     disp(EEG.nbchan);
% end
% 
% %%
% numMod = 8;
% for subjID = 1:10
%     subject = sprintf('nfle%d',subjID);
%     
%     try
%         EEG = pop_loadset(sprintf('sleep/%s.set',subject));
%         outdir = sprintf('/home/goodshawn12/MATLAB/AMICA/sleep/dataamicaout/%s/',subject);
%         if EEG.nbchan >= 5 && ~exist(outdir)
%             % run AMICA
%             [weights,sphere,mods] = runamica15c(EEG.data, 'num_models',numMod, 'outdir',outdir, ...
%                 'nodeproc',[8 24],'do_reject',1,'numrej',15,'rejsig',3,'rejint',1);
%             
%             modout = loadmodout15(outdir);
%             
%             load(sprintf('hyp%s.mat',subject));
%             % igure, plot(hyp(:,2),hyp(:,1));
%             
%             v = zeros(size(hyp,1)-1,numMod);
%             for it = 1:size(hyp,1)-1
%                 dataRange = hyp(it,2)*EEG.srate+1 : hyp(it+1,2)*EEG.srate;
%                 v(it,:) = mean(10.^modout.v(:,dataRange),2);
%             end
%             
%             figure,
%             subplot(5,1,1), plot(hyp(2:end,2)/60,hyp(2:end,1));
%             subplot(5,1,2:5), plot(hyp(2:end,2)/60,v);
%             xlabel('Time (min)'); ylabel('Model Probablity');
%             legend('Model 1','Model 2','Model 3','Model 4','Model 5','Model 6','Model 7','Model 8');
%             set(gcf,'Units','centimeters','Position',[5 5 30 20]);
%             
%             export_fig(sprintf('AMICA_8mod_%s',subject),'-png','-transparent');
%         end
%     end
% end
% 
% 
% %% Testing different number of models
% subjID = 3;
% subject = sprintf('n%d',subjID);
% 
% EEG = pop_loadset(sprintf('%s.set',subject));
% 
% for numMod = 20
%     outdir = sprintf('/home/goodshawn12/MATLAB/2017_AMICA_sleep/dataamicaout/%s_%dmod/',subject,numMod);
%     if EEG.nbchan >= 5 && ~exist(outdir)
%         % run AMICA
%         [weights,sphere,mods] = runamica15c(EEG.data, 'num_models',numMod, 'outdir',outdir, 'nodeproc',[8 24]); %,'do_reject',1,'numrej',15,'rejsig',3,'rejint',1);
%     end
% end
% 
% % plot log likelihood
% subject = sprintf('n%d',subjID);
% numModel = 12;
% maxStep = 2000;
% LL = zeros(maxStep,numModel);
% modProb = zeros(numModel,numModel);
% for modelID = 1:numModel
%     outdir = sprintf('/home/goodshawn12/MATLAB/2017_AMICA_sleep/dataamicaout/%s_%dmod/',subject,modelID);
%     modout = loadmodout15(outdir);
%     LL(:,modelID) = modout.LL;
%     modProb(1:modelID,modelID) = modout.mod_prob;
% end
% 
% % compute consistency of convergence
% numIter = 1;
% finalLL = zeros(numIter,numModel);
% finalLLIndex = zeros(numIter,numModel);
% meanLL = zeros(maxStep,numModel);
% stdLL = zeros(maxStep,numModel);
% for modelID = 1:numModel
%     numConvergeIndex = numIter*ones(1,maxStep);
%     for it = 1:numIter
%         divergeIndex = find(LL(:,modelID)==0,1);
%         if isempty(divergeIndex)
%             finalLL(it,modelID) = LL(maxStep,modelID);
%             finalLLIndex(it,modelID) = maxStep;
%         else
%             finalLL(it,modelID) = LL(divergeIndex-1,modelID);
%             numConvergeIndex(divergeIndex:end) = numConvergeIndex(divergeIndex:end) - 1;
%             finalLLIndex(it,modelID) = divergeIndex - 1;
%         end
%     end
%     meanLL(:,modelID) = LL(:,modelID)./numConvergeIndex';
% end
% 
% % plot results
% figure, plot(meanLL,'LineWidth',2)
% xlabel('Iterations'); ylabel('Log Likelihood'); legend('1 model','2 models','3 models','4 models','5 models','6 models','7 models','8 models','location','SouthEast');
% set(gca,'fontsize',12);
% export_fig results/LL_convergence_zoomin -png -transparent;
% 
% modelRange = 1:12;
% figure, plot(finalLL,'LineWidth',2);
% xlabel('Number of Models'); ylabel('Log Likelihood'); set(gca,'fontsize',12);
% export_fig results/LLxNumModel -png -transparent;
% 
% figure, h = errorbar(modelRange,mean(finalLL(:,modelRange)),std(finalLL(:,modelRange)));%/sqrt(numIter));
% xlabel('Number of Models'); ylabel('Converged Log Likelihood');
% set(gca,'fontsize',12); set(h,'Linewidth',2)
% export_fig results/LLxNumModel_std -png -transparent;
% 
% % model probabilities rh
% numSum = 7;
% sumOfTopProbModel = zeros(numModel,numSum);
% sumofBottomProbModel = zeros(numModel,numSum);
% for sumID = 1:numSum
%     for modelID = 1:numModel
%         sumOfTopProbModel(modelID,sumID) = sum(modProb(1:min(modelID,sumID),modelID));
%         sumofBottomProbModel(modelID,sumID) = sum(modProb(modelID-min(modelID,sumID)+1:modelID,modelID));
%     end
% end
% 
% modelRange = 1:numModel;
% figure, plot(modelRange,sumOfTopProbModel,'LineWidth',2);
% xlabel('Number of Models'); ylabel('Sum of Model Probabilities');
% set(gca,'fontsize',12); legend('1 model','2 models','3 models','4 models','5 models','6 models','7 models');
% 
% figure, plot(modelRange,sumofBottomProbModel,'LineWidth',2);
% xlabel('Number of Models'); ylabel('Cumlative Sum of Smallest Model Prob');
% set(gca,'fontsize',12); legend('1 model','2 models','3 models','4 models','5 models','6 models','7 models');
% 
% 
% 
% %% Visualize brain states transition in model likelihood PCA space
% numMod = 8;
% for subjID = 1:10
%     subject = sprintf('nfle%d',subjID);
% %     subject = sprintf('n%d',subjID);
%     
%     try
%         EEG = pop_loadset(sprintf('sleep/%s.set',subject));
%         outdir = sprintf('/home/goodshawn12/MATLAB/AMICA/sleep/dataamicaout/%s/',subject);
%         modout = loadmodout15(outdir);
%         load(sprintf('hyp%s.mat',subject));
%         load(sprintf('micro_str%s',subject));
%         
%         % correct EEG and hyp times
%         offset = start_time.h*60 + start_time.m + start_time.s/60 - EEG.etc.T0(4)*60 - EEG.etc.T0(5) - EEG.etc.T0(6)/60 - 0.5; % min (Hypnogram precedes EEG data)
%         if start_time.h < 6 % sleep time is after midnight
%             offset = offset + 24*60;
%         end
%         disp('Corrected offset time: '); disp(offset);
%         
%         % non-overlapping sliding window average of model probabilities 
%         v = zeros(numMod,size(hyp,1)-1);
%         timePoint = 1;
%         while timePoint < size(hyp,1) && (hyp(timePoint,2)+offset*60+30)*EEG.srate <= size(modout.v,2)
%             dataRange = round(hyp(timePoint,2)+offset*60)*EEG.srate+1 : round(hyp(timePoint+1,2)+offset*60)*EEG.srate;
%             v(:,timePoint) = mean(10.^modout.v(:,dataRange),2);
%             timePoint = timePoint+1;
%         end
%         
%         figure,
%         subplot(5,1,1), plot(hyp(1:timePoint-1,2)/60,hyp(1:timePoint-1,1));
%         subplot(5,1,2:5), plot(hyp(1:timePoint-1,2)/60,v(:,1:timePoint-1));
%         xlabel('Time (min)'); ylabel('Model Probablity');
%         legend('Model 1','Model 2','Model 3','Model 4','Model 5','Model 6','Model 7','Model 8');
%         set(gcf,'Units','centimeters','Position',[5 5 30 20]);
%         export_fig(sprintf('AMICA_8mod_%s',subject),'-png','-transparent');
%         close
%         
%         [V,D] = eig(v(:,1:timePoint-1)*v(:,1:timePoint-1)');
%         sphere = diag(1./diag(D.^0.5)) * V';
%         y = sphere * v(:,1:timePoint-1);
%         
%         % plot first 3 PC
%         figure, scatter3(y(end,1:timePoint-1), y(end-1,1:timePoint-1), y(end-2,1:timePoint-1),15,hyp(1:timePoint-1,1),'filled'); colorbar; view(-65,30)
%         xlabel('PC1'); ylabel('PC2'); zlabel('PC3'); set(gca,'fontsize',14);
%         set(gcf,'Units','centimeters','Position',[5 5 25 20]);
%         export_fig(sprintf('PC3D_%s',subject),'-png','-transparent');
%         close
%         
%         % plot first 2 PC
%         figure, scatter(y(end,1:timePoint-1), y(end-1,1:timePoint-1),15,hyp(1:timePoint-1,1),'filled'); colorbar;
%         xlabel('PC1'); ylabel('PC2');set(gca,'fontsize',14);
%         set(gcf,'Units','centimeters','Position',[5 5 25 20]);
%         export_fig(sprintf('PC2D_%s',subject),'-png','-transparent');
%         close
%         
%     end
% end
% 
% 
% %% movie 3D
% subject = 'n1';
% c =  hyp(1:timePoint-1,1);
% fid = figure;
% writerObj = VideoWriter(sprintf('Movie_PC3D_%s.avi',subject));
% writerObj.FrameRate = 20;
% open(writerObj);
% 
% figure(fid); set(gcf,'Units','centimeters','Position',[5 5 25 20]);
% ax1 = subplot(4,1,1); plot(hyp(1:timePoint-1,2)/60, hyp(1:timePoint-1,1),'lineWidth',1.5); 
% xlabel('Time (min)'); ylabel('Hypnogram (stage ID)'); set(gca,'FontSize', 9, 'FontWeight','bold');
% l = line([0, 0], ylim,'LineWidth',2,'Color','r');
% ax2 = subplot(4,1,2:4);
% caxis([0 5]); colorbar('FontSize',11,'YTick',[0,1,2,3,4,5],'YTickLabel',[0,1,2,3,4,5]);
% xlabel('PC1'); ylabel('PC2'); zlabel('PC3'); xlim([0 0.12]); ylim([-0.06 0.06]); zlim([-0.1 0.2]); set(gca,'fontsize',12,'FontWeight','bold');
% view(-65,30);
% 
% for i = 1:timePoint-1
%     figure(fid);
%     hold(ax1,'on'); set(l,'XData',[hyp(i,2)/60, hyp(i,2)/60]);
%     hold(ax2,'on'); scatter3(y(end,i), y(end-1,i), y(end-2,i),15,hyp(i+1,1),'filled'); 
%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);
% end
% hold off
% close(writerObj);
% 
% 
% %% matching ICA models and sleep stages
% subjID = 3;
% subject = sprintf('n%d',subjID);
% 
% EEG = pop_loadset(sprintf('sleep/%s.set',subject));
% outdir = sprintf('/home/goodshawn12/MATLAB/AMICA/sleep/dataamicaout/%s/',subject);
% modout = loadmodout15(outdir);
% load(sprintf('hyp%s.mat',subject));
% load(sprintf('micro_str%s',subject));
% 
% % correct EEG and hyp times
% offset = start_time.h*60 + start_time.m + start_time.s/60 - EEG.etc.T0(4)*60 - EEG.etc.T0(5) - EEG.etc.T0(6)/60; % min (Hypnogram precedes EEG data)
% if start_time.h < 6 % sleep time is after midnight
%     offset = offset + 24*60;
% end
% disp('Corrected offset time: '); disp(offset);
% 
% % non-overlapping sliding window average of model probabilities
% v = zeros(numMod,size(hyp,1)-1);
% timePoint = 1;
% while timePoint < size(hyp,1) && (hyp(timePoint,2)+offset*60+30)*EEG.srate <= size(modout.v,2)
%     dataRange = round(hyp(timePoint,2)+offset*60)*EEG.srate+1 : round(hyp(timePoint+1,2)+offset*60)*EEG.srate;
%     v(:,timePoint) = mean(10.^modout.v(:,dataRange),2);
%     timePoint = timePoint+1;
% end
% 
% 
% % plot first 3 models
% figure, scatter3(v(1,:), v(2,:), v(3,:),15,hyp(2:end,1),'filled'); colorbar; view(-65,30)
% xlabel('Model 1'); ylabel('Model 2'); zlabel('Model 3'); set(gca,'fontsize',14);
% set(gcf,'Units','centimeters','Position',[5 5 25 20]);
% export_fig(sprintf('Model3D_%s',subject),'-png','-transparent');
% close
% 
% % plot first 2 models
% figure, scatter(v(1,:), v(2,:),15,hyp(2:end,1),'filled'); colorbar;
% xlabel('Model 1'); ylabel('Model 2'); set(gca,'fontsize',14);
% set(gcf,'Units','centimeters','Position',[5 5 25 20]);
% export_fig(sprintf('Model2D_%s',subject),'-png','-transparent');
% close
% 
% % sorting model probabilities based on sleep stages
% stages = unique(hyp(:,1));
% stageProb = zeros(length(stages),size(v,1));
% for it = 1:length(stages)
%     stageTimeIndex = find(hyp(2:end,1) == stages(it));
%     stageProb(it,:) = mean(v(:,stageTimeIndex)');
% end
% 
% figure, hold on; bar(stageProb);
% set(gca,'XTickLabel',{'Awake','N1','N2','N3','N4','REM'},'XTick',1:size(stageProb,1)); xlabel('Sleep Stages'); ylabel('Model Probabilities'); set(gca,'fontsize',14,'fontweight','bold'); set(gcf,'Units','centimeters','Position',[5 5 25 15]);
% legend('Model 1','Model 2','Model 3','Model 4','Model 5','Model 6','Model 7','Model 8','Location','northeastoutside');
% export_fig(sprintf('ModelProbXSleepStages_%s',subject),'-png','-transparent'); close
% 
% figure, hold on; bar(stageProb'); xlabel('ICA Model Index'); ylabel('Model Probabilities'); 
% set(gca,'XTickLabel',1:size(stageProb,2),'XTick',1:size(stageProb,2)); set(gca,'fontsize',14,'fontweight','bold'); set(gcf,'Units','centimeters','Position',[5 5 25 15]);
% legend('Awake','N1','N2','N3','N4','REM','Location','northeastoutside');
% export_fig(sprintf('StageProbXModels_%s',subject),'-png','-transparent'); close
% 
% 
% %% plot scalp maps and source activities of each ICA model
% channelLabels = {};
% for it = 1:EEG.nbchan
%     channelPairs = strsplit(EEG.chanlocs(it).labels,'-');
%     channelLabels = unique({channelLabels{:}, channelPairs{:}});
% end
% 
% transformationMat = zeros(EEG.nbchan);
% for it = 1:EEG.nbchan
%     channelPairs = strsplit(EEG.chanlocs(it).labels,'-');
%     transformationMat(it,ismember(channelLabels,channelPairs(1))) = 1;
%     transformationMat(it,ismember(channelLabels,channelPairs(2))) = -1;
% end
% invTranMat = pinv(transformationMat);
% 
% EEGtemp = eeg_emptyset;
% EEGtemp.data = invTranMat * EEG.data;
% EEGtemp.srate = EEG.srate;
% EEGtemp = eeg_checkset(EEGtemp);
% EEGtemp.chanlocs = struct('labels',channelLabels);
% EEGtemp = pop_chanedit(EEGtemp, 'lookup','/data/common/matlab/eeglab/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');
% 
% EEGtemp.icasphere = modout.S * transformationMat;
% 
% for modelIndex = 1:modout.num_models
%     EEGtemp.icaweights = modout.W(:,:,modelIndex);
%     EEGtemp.icawinv = invTranMat * modout.A(:,:,modelIndex);
%     
%     pop_topoplot(EEGtemp,0,1:size(EEGtemp.icaweights,1))
%     export_fig(sprintf('ICA_%s_M%d',subject,modelIndex),'-png','-transparent'); close
% end

