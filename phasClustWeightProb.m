% makes a histogram of phases with the bump PDF and the pahse of the
% analytic signal at a frequency band. It uses the results from a HSMM
% computed from the PCs of all subjects concatenated

clear all
close all
clc

nBump = 5;                          % number of bumps in HSMM
plotHist = 0;                       % plot head histogram per subject 
nBins = 18*2;                         % number of bins in the phase histogram
edges = linspace(-pi, pi, nBins+1);  % Does nBins from -pi to pi Bin edges
width = abs(edges(1) - edges(2));        % width of bins
binVal = edges(1:end-1) - diff(edges)/2;    % value of each bin
nCh = 32;                           % number of channels
freqBand = 'theta';                 % frequency band to be analyzed
nRep = 400;                         % number of repetition to build a random PDF

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

pValBpcTr = zeros(nCh,nBump,nSj);
pValBpcAll = zeros(nCh,nBump,nSj);
pValZbpcTr = zeros(nCh,nBump,nSj);
pValZbpcAll = zeros(nCh,nBump,nSj);

bpcAllScMod = zeros(nBump,nSj);
bpcAllScArg = zeros(nBump,nSj);
bpcTrScMod = zeros(nBump,nSj);
bpcTrScArg = zeros(nBump,nSj);

pValBpcScTr = zeros(nBump,nSj);
pValBpcScAll = zeros(nBump,nSj);
pValZbpcScTr = zeros(nBump,nSj);
pValZbpcScAll = zeros(nBump,nSj);

bpcTrMod = zeros(nCh,nBump,nSj);
bpcTrArg = zeros(nCh,nBump,nSj);
bpcAllMod = zeros(nCh,nBump,nSj);
bpcAllArg = zeros(nCh,nBump,nSj);

for sj = 1:nSj
    % load datasets
    data = load([pathdata namesHilb(sj).name], freqBand, 'x', 'y');
    x = data.x;
    y = data.y;
    % do pahse angle from analytic signal
    dataPha = angle(data.(freqBand)); % argment of analytic signal
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
    % Phase clustering at transitions with a mean parameter per trial
    trMod = zeros(length(y),nCh,nBump);
    trArg = zeros(length(y),nCh,nBump);
    for tr = 1:length(y)
        for ch = 1:nCh
            for bu = 1:nBump
                trMod(tr,ch, bu) = abs(mean(probump(x(tr):y(tr),bu) .* exp(1i*dataPha(x(tr):y(tr),ch))));
                trArg(tr,ch, bu) = angle(mean(probump(x(tr):y(tr),bu) .* exp(1i*dataPha(x(tr):y(tr),ch))));
            end
        end
    end
    for ch = 1:nCh
        for bu = 1:nBump
            bpcTrMod(ch, bu,sj) = abs(mean(trMod(:,ch,bu) .* exp(1i*trArg(:,ch,bu))));
            bpcTrArg(ch, bu,sj) = angle(mean(trMod(:,ch,bu) .* exp(1i*trArg(:,ch,bu))));
        end
    end
    % phase consistency across the scalp
    bpcAllScMod(:,sj) = abs(mean(exp(1i*bpcAllArg(:,:,sj)),1));
    bpcAllScArg(:,sj) = angle(mean(exp(1i*bpcAllArg(:,:,sj)),1));
    bpcTrScMod(:,sj) = abs(mean(exp(1i*bpcTrArg(:,:,sj)),1));
    bpcTrScArg(:,sj) = abs(mean(exp(1i*bpcTrArg(:,:,sj)),1));
    
    

    
    % %%%%%%%%% -------- Randomization ---------- %%%%%%%%%%%%%%%%%%%%
    
    % Shuffle time points at random locations for each trial
    bpcAllModRa = zeros(nCh,nBump,nRep);
    bpcTrModRa = zeros(nCh,nBump,nRep);
    bpcAllScModRa = zeros(nBump,nRep);
    bpcTrScModRa = zeros(nBump,nRep);
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
            end
        end
        % Phase clustering at transitions with a mean parameter per trial
        trMod = zeros(length(y),nCh,nBump);
        trArg = zeros(length(y),nCh,nBump);
        for tr = 1:length(y)
            for ch = 1:nCh
                for bu = 1:nBump
                    trMod(tr,ch, bu) = abs(mean(probump(x(tr):y(tr),bu) .* exp(1i*dataRand(x(tr):y(tr),ch))));
                    trArg(tr,ch, bu) = angle(mean(probump(x(tr):y(tr),bu) .* exp(1i*dataRand(x(tr):y(tr),ch))));
                end
            end
        end
        for ch = 1:nCh
            for bu = 1:nBump
                bpcTrModRa(ch, bu,r) = abs(mean(trMod(:,ch,bu) .* exp(1i*trArg(:,ch,bu))));
            end
        end
        clear dataRand trMod trArg
        % phase consistency across the scalp
        bpcAllScModRa(:,r) = abs(mean(exp(1i*bpcAllModRa(:,:,r)),1));
        bpcTrScModRa(:,r) = abs(mean(exp(1i*bpcTrModRa(:,:,r)),1));

    end
     % Test of random values (null hypo) respect to empirical    
    for bu = 1:nBump,
        for ch = 1:nCh,
            pValBpcTr(ch,bu,sj) = length(find(bpcTrModRa(ch,bu,:) >= bpcTrMod(ch,bu,sj))) / nRep;
            pValBpcAll(ch,bu,sj) = length(find(bpcAllModRa(ch,bu,:) >= bpcAllMod(ch,bu,sj))) / nRep;
            % z-score p-value
            pValZbpcTr(ch,bu,sj) = 1 - normcdf((bpcTrMod(ch,bu,sj) - mean(bpcTrModRa(ch,bu,:))) ./ std(bpcTrModRa(ch,bu,:)));
            pValZbpcAll(ch,bu,sj) = 1 - normcdf((bpcAllMod(ch,bu,sj) - mean(bpcAllModRa(ch,bu,:))) ./ std(bpcAllModRa(ch,bu,:)));  
        end
        % p-value of scalp consistency
        pValBpcScTr(bu,sj) = length(find(bpcTrScModRa(bu,:) >= bpcTrScMod(bu,sj))) / nRep;
        pValBpcScAll(bu,sj) = length(find(bpcAllScModRa(bu,:) >= bpcAllScMod(bu,sj))) / nRep;
        % z-score 
        pValZbpcScAll(bu,sj) = 1 - normcdf((bpcAllScMod(bu,sj) - mean(bpcAllScModRa(bu,:))) ./ std(bpcAllScModRa(bu,:)));
        pValZbpcScTr(bu,sj) = 1 - normcdf((bpcTrScMod(bu,sj) - mean(bpcTrScModRa(bu,:))) ./ std(bpcTrScModRa(bu,:)));
    end
    clear dataPha bpcTrModRa bpcAllModRa bpcTrScModRa bpcAllScModRa
%     %%%%%%%%%%%%%% ---- visual representation -----
  
end

savepath = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/histPhase/';
nameSave = [savepath 'wPhaseClust' num2str(nBump) 'Bu.mat'];
save(nameSave,'info','pValBpcModAll','pValBpcModTr','pValBpcTr','pValBpcAll','bpcAllScMod','bpcAllScArg','bpcTrScMod','bpcTrScArg', ...
              'bpcTrMod','bpcTrArg','bpcAllMod','bpcAllArg');