% Weights the phases with the bump PDF. It uses the results from a HSMM
% computed from the PCs of all subjects concatenated. It computes the mean
% phease at eact trial and each subeject for each bump

clear all
close all
clc

nBump = 5;                          % number of bumps in HSMM
plotHist = 0;                       % plot head histogram per subject 
nBins = 18*2;                         % number of bins in the phase histogram
nCh = 32;                           % number of channels
freqBand = 'beta';                 % frequency band to be analyzed
nRep = 200;                         % number of repetition to build a random PDF

% prepare data paths to load
pathdata = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/';
%varElec = [pathdata 'electrodes.mat'];
namePDF = ['AllSj100HSMM' num2str(nBump) 'Bout235.mat'];
codeHilb = ['*Bads100.mat'];
namesHilb = dir([pathdata codeHilb]);
nSj = length(namesHilb);

% test that the subject are sorted equaly so loaded datasets agree
pathdataTest = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/';
snamesTest = dir([pathdataTest '*epochs235.mat']);
if length(snamesTest) ~= length(namesHilb), error('Unbalancend number of files per subject'), end
for sj = 1:nSj
    idTest = strfind(snamesTest(sj).name ,'_');
    idHil = strfind(namesHilb(sj).name ,'_');
    if namesHilb(sj).name(1:idHil) ~= snamesTest(sj).name(1:idTest)
        error('Error: Subjects are not the same')
    end
end
clear pathdataTest snamesTest idHil idTest codeHilb

% load variables
load([pathdata 'electrodes.mat'])
load([pathdata namePDF], 'eventprobs')

info = ['1) trBPC trial Bump Phase Clustering. Phase is wighted by bumps PDF ' ...
        ' at each trial a mean trBPC. Then the mean of all trBPC is tested by randomization' ...
        ' 2) sjBPC. The same but the mean es calculated over all trials and points in a subject weighted by the PDF' ...
        ' one HSMMcalculated for all subject concatenated.Fs:100Hz,Band: ' freqBand];

% pValBpcTr = zeros(nCh,nBump,nSj);
pValBpcAll = zeros(nCh,nBump,nSj);
% pValZbpcTr = zeros(nCh,nBump,nSj);
pValZbpcAll = zeros(nCh,nBump,nSj);
% pValBpcTrW = zeros(nCh,nBump,nSj);
% pValZbpcTrW = zeros(nCh,nBump,nSj);

bpcAllScMod = zeros(nBump,nSj);
bpcAllScArg = zeros(nBump,nSj);
% bpcTrScMod = zeros(nBump,nSj);
% bpcTrScModX = zeros(nBump,nSj);
% bpcTrScArg = zeros(nBump,nSj);
% bpcTrScArgW = zeros(nBump,nSj);
% bpcTrScModW = zeros(nBump,nSj);

% pValBpcScTr = zeros(nBump,nSj);
pValBpcScAll = zeros(nBump,nSj);
% pValZbpcScTr = zeros(nBump,nSj);
pValZbpcScAll = zeros(nBump,nSj);
% pValBpcScTrW = zeros(nBump,nSj);
% pValZbpcScTrW = zeros(nBump,nSj);
% pValZbpcScTrX = zeros(nBump,nSj);
% pValBpcScTrX = zeros(nBump,nSj);

% bpcTrMod = zeros(nCh,nBump,nSj);
% bpcTrArg = zeros(nCh,nBump,nSj);
% bpcTrModW = zeros(nCh,nBump,nSj);
% bpcTrArgW = zeros(nCh,nBump,nSj);
bpcAllMod = zeros(nCh,nBump,nSj);
bpcAllArg = zeros(nCh,nBump,nSj);

for sj = 1:nSj
    % load datasets
    data = load([pathdata namesHilb(sj).name], freqBand, 'x', 'y');
    x = data.x;
    y = data.y;
    % do pahse angle from analytic signal
    dataPha = angle(data.(freqBand)); % argment of analytic signal
    dataZx = data.(freqBand);
    clear data
    probump = [];
    % Concatenate trials vertically to remove zeros out of trial
    for tr = 1:length(y)
        yTr = y(tr) - x(tr) + 1;
        probump = vertcat(probump, squeeze(eventprobs(1:yTr,tr,:)));
    end
    % Remove trials from current subject
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --------------- ATTENTION!!!
    eventprobs(:,1:length(y),:) = [];
    %clear eventprobs % Only for the case of ONE Subject
    
    % phase clusering at transtion for all trials taken at once
    for ch = 1:nCh
        for bu = 1:nBump
            bpcAllMod(ch, bu,sj) = squeeze(abs(mean(probump(:,bu) .* exp(1i*dataPha(:,ch)))));
            bpcAllArg(ch, bu,sj) = squeeze(angle(mean(probump(:,bu) .* exp(1i*dataPha(:,ch)))));
        end
    end
    % phase consistency across the scalp
    bpcAllScMod(:,sj) = abs(mean(exp(1i*bpcAllArg(:,:,sj)),1));
    bpcAllScArg(:,sj) = angle(mean(exp(1i*bpcAllArg(:,:,sj)),1));
    
%     % Phase clustering at transitions with a mean parameter per trial
    trMod = zeros(length(y),nCh,nBump);
    trArg = zeros(length(y),nCh,nBump);
%     trScModW = zeros(length(y),nBump);
%     trScArgW = zeros(length(y),nBump);
%     trScMod = zeros(length(y),nBump);
% %     trScArg = zeros(length(y),nBump);
    for tr = 1:length(y)
        for ch = 1:nCh
            for bu = 1:nBump
                trMod(tr,ch, bu) = abs(mean(probump(x(tr):y(tr),bu) .* exp(1i*dataPha(x(tr):y(tr),ch))));
                trArg(tr,ch, bu) = angle(mean(probump(x(tr):y(tr),bu) .* exp(1i*dataPha(x(tr):y(tr),ch))));
            end
        end
%         % phase clustering accross the scap per trial
%         trScModW(tr,:) = abs(mean(trMod(tr,:,:) .* exp(1i*trArg(tr,:,:))));
%         trScArgW(tr,:) = angle(mean(trMod(tr,:,:) .* exp(1i*trArg(tr,:,:))));
%         trScMod(tr,:) = abs(mean(exp(1i*trArg(tr,:,:))));
% %         trScArg(tr,:) = angle(mean(exp(1i*trArg(tr,:,:))));
    end
%     % transition phase clustering
%     bpcTrModW(:,:,sj) = abs(mean(trMod .* exp(1i*trArg), 1));
%     bpcTrArgW(:,:,sj) = angle(mean(trMod .* exp(1i*trArg), 1));
%     bpcTrMod(:,:,sj) = abs(mean(exp(1i*trArg), 1));
%     bpcTrArg(:,:,sj) = angle(mean(exp(1i*trArg), 1));
%     % trial's mean phase clusterin accros the scalp 
%     bpcTrScModW(:,sj) = abs(mean(trScModW .* exp(1i*trScArgW),1));
%     bpcTrScArgW(:,sj) = angle(mean(trScModW .* exp(1i*trScArgW),1));
% %     bpcTrScMod(:,sj) = abs(mean(exp(1i*trScArg),1));
%     bpcTrScModX(:,sj) = mean(trScMod,1);
% %     bpcTrScArg(:,sj) = angle(mean(exp(1i*trScArg),1));
%     clear trScArg 
    % %%%%%%%%% -------- Randomization ---------- %%%%%%%%%%%%%%%%%%%%
    
    % Shuffle time points at random locations for each trial
    bpcAllModRa = zeros(nCh,nBump,nRep);
    bpcAllArgRa = zeros(nCh,nBump,nRep);
    bpcAllScModRa = zeros(nBump,nRep);
%     bpcTrModRa = zeros(nCh,nBump,nRep);
%     bpcTrModWra = zeros(nCh,nBump,nRep);
%     bpcTrScModWra = zeros(nBump,nRep);
%     bpcTrScModRa = zeros(nBump,nRep);
    bpcTrScModXra = zeros(nBump,nRep);
    for r =1:nRep,
        dataRand = [];   
        for tr = 1:length(y),
            idr = randi([y(tr)-x(tr)]);
            dataRand = vertcat(dataRand, dataPha(x(tr)+idr:y(tr),:));
            dataRand = vertcat(dataRand, dataPha(x(tr):x(tr)+idr-1,:));
        end
         % phase clusering at transtion for all trials taken at once
        for ch = 1:nCh
            for bu = 1:nBump
                bpcAllModRa(ch, bu,r) = squeeze(abs(mean(probump(:,bu) .* exp(1i*dataRand(:,ch)))));
                bpcAllArgRa(ch, bu,r) = squeeze(angle(mean(probump(:,bu) .* exp(1i*dataRand(:,ch)))));
            end
        end
            % phase consistency across the scalp
        bpcAllScModRa(:,r) = abs(mean(exp(1i*bpcAllArgRa(:,:,r)),1));
%         % Phase clustering at transitions with a mean parameter per trial
%         trMod = zeros(length(y),nCh,nBump);
%         trArg = zeros(length(y),nCh,nBump);
%         trScModW = zeros(length(y),nBump);
%         trScArgW = zeros(length(y),nBump);
%         trScMod = zeros(length(y),nBump);
%         for tr = 1:length(y)
%             for ch = 1:nCh
%                 for bu = 1:nBump
%                     trMod(tr,ch, bu) = abs(mean(probump(x(tr):y(tr),bu) .* exp(1i*dataRand(x(tr):y(tr),ch))));
%                     trArg(tr,ch, bu) = angle(mean(probump(x(tr):y(tr),bu) .* exp(1i*dataRand(x(tr):y(tr),ch))));
%                 end
%             end
%                 % phase clustering accross the scap per trial
%             trScModW(tr,:) = abs(mean(trMod(tr,:,:) .* exp(1i*trArg(tr,:,:))));
%             trScArgW(tr,:) = angle(mean(trMod(tr,:,:) .* exp(1i*trArg(tr,:,:))));
%             trScMod(tr,:) = abs(mean(exp(1i*trArg(tr,:,:))));
% %             trScArg(tr,:) = angle(mean(exp(1i*trArg(tr,:,:))));
%         end
%         % transition phase clustering
%         bpcTrModWra(:,:,r) = abs(mean(trMod .* exp(1i*trArg), 1));
%         bpcTrModRa(:,:,r) = abs(mean(exp(1i*trArg), 1));
%         % trial's mean phase clusterin accros the scalp 
%         bpcTrScModWra(:,r) = abs(mean(trScModW .* exp(1i*trScArgW),1));
% %         bpcTrScModRa(:,r) = abs(mean(exp(1i*trScArg),1));
%         bpcTrScModXra(:,r) = mean(trScMod,1);
    end
%     clear dataPha dataRand trMod trArg trScModW trScArgW trScMod probump bpcAllArgRa
    clear dataPha dataRand probump bpcAllArgRa

     % Test of random values (null hypo) respect to empirical    
    for bu = 1:nBump,
        for ch = 1:nCh,
%             pValBpcTr(ch,bu,sj) = length(find(bpcTrModRa(ch,bu,:) >= bpcTrMod(ch,bu,sj))) / nRep;
%             pValBpcTrW(ch,bu,sj) = length(find(bpcTrModWra(ch,bu,:) >= bpcTrModW(ch,bu,sj))) / nRep;
            pValBpcAll(ch,bu,sj) = length(find(bpcAllModRa(ch,bu,:) >= bpcAllMod(ch,bu,sj))) / nRep;
            % z-score p-value
%             pValZbpcTrW(ch,bu,sj) = 1 - normcdf((bpcTrModW(ch,bu,sj) - mean(bpcTrModWra(ch,bu,:))) ./ std(bpcTrModWra(ch,bu,:)));
%             pValZbpcTr(ch,bu,sj) = 1 - normcdf((bpcTrMod(ch,bu,sj) - mean(bpcTrModRa(ch,bu,:))) ./ std(bpcTrModRa(ch,bu,:)));
            pValZbpcAll(ch,bu,sj) = 1 - normcdf((bpcAllMod(ch,bu,sj) - mean(bpcAllModRa(ch,bu,:))) ./ std(bpcAllModRa(ch,bu,:)));  
        end
        % p-value of scalp consistency
%         pValBpcScTr(bu,sj) = length(find(bpcTrScModRa(bu,:) >= bpcTrScMod(bu,sj))) / nRep;
%         pValBpcScTrX(bu,sj) = length(find(bpcTrScModXra(bu,:) >= bpcTrScModX(bu,sj))) / nRep;
%         pValBpcScTrW(bu,sj) = length(find(bpcTrScModWra(bu,:) >= bpcTrScModW(bu,sj))) / nRep;
        pValBpcScAll(bu,sj) = length(find(bpcAllScModRa(bu,:) >= bpcAllScMod(bu,sj))) / nRep;
        % z-score 
        pValZbpcScAll(bu,sj) = 1 - normcdf((bpcAllScMod(bu,sj) - mean(bpcAllScModRa(bu,:))) ./ std(bpcAllScModRa(bu,:)));
%         pValZbpcScTr(bu,sj) = 1 - normcdf((bpcTrScMod(bu,sj) - mean(bpcTrScModRa(bu,:))) ./ std(bpcTrScModRa(bu,:)));
%         pValZbpcScTrX(bu,sj) = 1 - normcdf((bpcTrScModX(bu,sj) - mean(bpcTrScModXra(bu,:))) ./ std(bpcTrScModXra(bu,:)));
%         pValZbpcScTrW(bu,sj) = 1 - normcdf((bpcTrScModW(bu,sj) - mean(bpcTrScModWra(bu,:))) ./ std(bpcTrScModWra(bu,:)));
    end
    clear bpcAllModRa bpcAllScModRa
end

savepath = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/histPhase/';
nameSave = [savepath 'wPhaseClust' num2str(nBump) 'Bu' freqBand '.mat'];
% save(nameSave,'info','pValBpcTr','pValBpcAll','pValZbpcTr','pValZbpcAll','pValBpcTrW', ...
%     'pValZbpcTrW','bpcAllScMod','bpcAllScArg','bpcTrScArgW', 'pValBpcScTrX',...
%     'bpcTrScModW','pValBpcScAll','pValZbpcScAll','pValZbpcScTrX', ...
%     'pValBpcScTrW','pValZbpcScTrW','bpcTrMod','bpcTrArg','bpcTrModW','bpcTrArgW','bpcAllMod','bpcAllArg');

save(nameSave,'info','pValBpcAll','pValZbpcAll', ...
    'bpcAllScMod','bpcAllScArg',...
    'pValBpcScAll','pValZbpcScAll', ...
    'bpcAllMod','bpcAllArg');