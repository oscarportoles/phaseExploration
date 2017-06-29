% Load band pass filtered data and band pass filterd prestimulus data. 
% concatenates both and do the hilbert transform per trial. Then it computes ITPC and Power 
% for each of the 50-msec time windows of anlysis. It requires "eventprob" variable from the 
% HSMM-MVPA output. "eventprob" contains the probability distribution of the location of a 
% bump. "eventprob" is used to estimate the location of the 50-msec time windows.

clear all
close all
clc

disp(['Starts at: ' datestr(now)])

nBump = 5;                          % number of bumps in HSMM
nCh = 32;                           % number of channels
fs = 100;                           % sampling frequency [Hz]
aPnts = 2;                        % number of points around the middle points of a stage to compute connectivty
nBand = 3;                          % number of frequency bands to analyze
win = 2*aPnts+1;                    % number of points in a window of interest

fBand = {'Theta','Alpha','Beta'};   % frequency band to be analyzed
varPos = {{'x','y','thetaEEG'}, ...
          {'x','y','alphaEEG'}, ...
          {'x','y','betaEEG'}};
varPre = {{'x','y','thetaEEGbl'}, ...
          {'x','y','alphaEEGbl'}, ...
          {'x','y','betaEEGbl'}};
      
% limits of each stage [ini, end]
stagLim = zeros(nBump+1,2);
stagLim(:,1) = 1:win:(win)*nBump+1;
stagLim(:,2) = win:win:(win)*(nBump+1);
            %/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20
pathdata = '/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/';
pathporb = '/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/';
savepath = '/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/PhaseConectiv/';
posnames = dir([pathdata '*_EpochBands100.mat']);
prenames = dir([pathdata '*PreStimEpochBands100.mat']);
namePDF = ['AllSj100HSMM' num2str(nBump) 'Bout235.mat'];
nSj = length(posnames);
info = {'Concatenates prestimulus data and poststimulus data. It does' ... 
        'Hilber transformed data (analytic signal) by EEG frequency bands,' ...
        'Computes inter trial phase clustering and power at pre-stimulus and middle of stage window' ...
        ' Fs = 100Hz, . Theta, Alpha and Beta bands'};

% test that the subjects are sorted equaly when loading so loaded subjects agree
pathdataTest = '/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/';
snamesTest = dir([pathdataTest '*epochs235.mat']);
if length(snamesTest) ~= length(posnames), error('Unbalancend number of files per subject'), end
for sj = 1:nSj
    idTest = strfind(snamesTest(sj).name ,'_');
    idpos = strfind(posnames(sj).name ,'_');
    idpre = strfind(prenames(sj).name ,'_');
    if posnames(sj).name(1:idpos) ~= snamesTest(sj).name(1:idTest) | ...
        prenames(sj).name(1:idpre) ~= snamesTest(sj).name(1:idTest) | ...
        posnames(sj).name(1:idpos) ~= prenames(sj).name(1:idpre)
        error('Error: Subjects are not the same')
    end
end
clear pathdataTest snamesTest idpos idpre idTest


for fb = 1:nBand
    % load variables global variables
    eventprobs = load([pathporb namePDF], 'eventprobs');
    eventprobs = eventprobs.eventprobs;

    for sj = 1:nSj
        % load data sets
        dataout = [];                           % put output data
        dataPos = load([pathdata posnames(sj).name], varPos{fb}{:});     % band specific Post-stimulus data
        xpos    = dataPos.(varPos{fb}{1});
        ypos    = dataPos.(varPos{fb}{2});
        dataPos = dataPos.(varPos{fb}{3});
        dataPre = load([pathdata prenames(sj).name], varPre{fb}{:});     % band specific Pre-stimulus data
        xpre    = dataPre.(varPos{fb}{1});
        ypre    = dataPre.(varPos{fb}{2});
        dataPre = dataPre.(varPre{fb}{3});
        nTr     = length(ypos);                     % number of trial in subject

        % pre-stimulus length or baselin
        yBl     = ypre(1) - xpre(1) + 1;        
        blIdx   = [1:yBl];                        % baseline indexes

        %powPre = zeros(yBl, nCh, nTr);
        powTra = zeros(nCh, nBump, nTr);
        powStg = zeros(nCh, nBump+1, nTr);
        %powInPre = zeros(yBl, nCh, nTr);
        powInTra = zeros(nCh, nBump, nTr);
        powInStg = zeros(nCh, nBump+1, nTr);
        itpcPre = zeros(yBl, nCh, nTr);
        itpcTra = zeros(win, nCh, nBump, nTr);
        itpcStg = zeros(win, nCh, nBump+1, nTr);
        
        % limit of random sample for time shifting
        raLim = min(ypos - xpos) - 2 + yBl;

        for tr = 1:nTr
            yTr = ypos(tr) - xpos(tr) + 1;
            % hilber Transform of pre and post stimulus data concatenated, one trial
            dataH = hilbert(vertcat(dataPre(xpre(tr):ypre(tr),:), dataPos(xpos(tr):ypos(tr),:)));
            % define stage bondaries:
            probs   = squeeze(eventprobs(1:yTr,tr,:));
            expLoS  = [1:yTr] * probs;   % bump expected locations
            expLo   = [1 expLoS yTr];     % trial edges
            expLoS  = round(expLoS);
            for bu=1:nBump+1
                % indexes of each stage
                meanSt  = round(mean(expLo(bu:bu+1))); 
                stagIdx = [meanSt-aPnts:meanSt+aPnts];
                if bu ~= nBump+1
                    tranIdx = [expLoS(bu)-aPnts:expLoS(bu)+aPnts];
                    tranIdx = tranIdx + yBl;
                end
                stagIdx = stagIdx + yBl;  % correct by samples in pre-stimulus
                % baseline
                if bu == 1
                    powPre = squeeze(mean(real(dataH(blIdx,:)).^2));
                    powPreStd = squeeze(std(real(dataH(blIdx,:)).^2));
                    powInPre = squeeze(mean(abs(dataH(blIdx,:)).^2));
                    powInPreStd = squeeze(std(abs(dataH(blIdx,:)).^2));
                    itpcPre(:,:,tr) = angle(dataH(blIdx,:));
                end
                % Power and phase per electrode (z-scored from baseline)
                
                powStg(:,bu,tr) = (median(real(dataH(stagIdx,:)).^2) - powPre) ./ powPreStd;
                powInStg(:,bu,tr) = (median(abs(dataH(stagIdx,:)).^2) - powInPre) ./ powInPreStd;
                itpcStg(:,:,bu,tr) = angle(dataH(stagIdx,:));
                if bu ~= nBump+1
                    powTra(:,bu,tr) = (median(real(dataH(tranIdx,:)).^2) - powPre) ./ powPreStd;
                    powInTra(:,bu,tr) = (median(abs(dataH(tranIdx,:)).^2) - powInPre) ./ powInPreStd;
                    itpcTra(:,:,bu,tr) = angle(dataH(tranIdx,:));
                end
            end
        end
        % Remove trials from current subject
        eventprobs(:,1:nTr,:) = [];
        % Inter-trial phase clustering  
        itpcPre = squeeze(abs(mean(exp(1i*itpcPre),3)));
        itpcPreStd = squeeze(std(itpcPre));
        itpcPre = squeeze(median(itpcPre));
        itpcStg = (squeeze(median(squeeze(abs(mean(exp(1i*itpcStg),4))))) - repmat(itpcPre',1,nBump+1));% ./ repmat(itpcPreStd',1,nBump+1);
        itpcTra = (squeeze(median(squeeze(abs(mean(exp(1i*itpcTra),4))))) - repmat(itpcPre',1,nBump));% ./ repmat(itpcPreStd',1,nBump);
        
        % mean power per subject
        %powPre = squeeze(mean(powPre,3));
        powStg = squeeze(mean(powStg,3));
        powTra = squeeze(mean(powTra,3));
        %powInPre = squeeze(mean(powInPre,3));
        powInStg = squeeze(mean(powInStg,3));
        powInTra = squeeze(mean(powInTra,3));

        idSj = regexp(posnames(sj).name, '_');
        namesave = [posnames(sj).name(1:idSj) 'itpcPowFre' fBand{fb} '.mat'];
        save([savepath namesave], 'itpcStg','itpcTra','info','powTra','powStg','powInTra','powInStg','nTr')
        disp(['-- Subject ', num2str(sj),' has finished. --'])
    end
    disp(['-- Done with frquency ' fBand{fb} ' at: ' datestr(now) ' --'])
end
