% Load band pass filtered data and band pass filterd prestimulus data. 
% concatenates both and do the hilbert transform per trial

clear all
close all
clc

nBump = 5;                          % number of bumps in HSMM
nCh = 32;                           % number of channels
freqBand = 'theta';                 % frequency band to be analyzed
fs = 100;                           % sampling frequency [Hz]
pair = combnk([1:nCh],2);        % all pairs of electrodes combinations [nComb, 2]
aPnts = 2;                        % number of points around the middle points of a stage to compute connectivty
doSurro = 1;                        % flag to do analysis of surrogate data
nRep = 300;                         % number of repetitions to generate surrogated data

% limits of each stage [ini, end]
stagLim = zeros(nBump+1,2);
stagLim(:,1) = 1:2*aPnts+1:(2*aPnts+1)*nBump+1;
stagLim(:,2) = 2*aPnts+1:2*aPnts+1:(2*aPnts+1)*(nBump+1);

% pathdata = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/';
% pathporb = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/';
% savepath = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/PhaseConectiv/';
posnames = dir('*_EpochBands100.mat');
prenames = dir('*PreStimEpochBands100.mat');
namePDF = ['AllSj100HSMM' num2str(nBump) 'Bout235.mat'];
nSj = length(posnames);
varPos = {'x','y','thetaEEG'};
varPre = {'x','y','thetaEEGbl'};
info = {'Concatenates prestimulus data and poststimulus data. It does' ... 
        'Hilber transformed data (analytic signal) by EEG frequency bands,' ...
        'Computes connectivity at pre-stimulus and middle of stage window, Fs = 100Hz,' ...
        ' baseline median and absolute median deviation'};
% 
% % test that the subjects are sorted equaly when loading so loaded subjects agree
% pathdataTest = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/';
% snamesTest = dir([pathdataTest '*epochs235.mat']);
% if length(snamesTest) ~= length(posnames), error('Unbalancend number of files per subject'), end
% for sj = 1:nSj
%     idTest = strfind(snamesTest(sj).name ,'_');
%     idpos = strfind(posnames(sj).name ,'_');
%     idpre = strfind(prenames(sj).name ,'_');
%     if posnames(sj).name(1:idpos) ~= snamesTest(sj).name(1:idTest) | ...
%         prenames(sj).name(1:idpre) ~= snamesTest(sj).name(1:idTest) | ...
%         posnames(sj).name(1:idpos) ~= prenames(sj).name(1:idpre)
%         error('Error: Subjects are not the same')
%     end
% end
% clear pathdataTest snamesTest idpos idpre idTest

% load variables global variables
load(namePDF, 'eventprobs')

for sj = 1:nSj,
    % load data sets
    dataout = [];                           % put output data
    dataPos = load(posnames(sj).name, varPos{:});     % band specific Post-stimulus data
    xpos    = dataPos.(varPos{1});
    ypos    = dataPos.(varPos{2});
    dataPos = dataPos.(varPos{3});
    dataPre = load(prenames(sj).name, varPre{:});     % band specific Pre-stimulus data
    xpre    = dataPre.(varPos{1});
    ypre    = dataPre.(varPos{2});
    dataPre = dataPre.(varPre{3});
    nTr     = length(ypos);                     % number of trial in subject
    
    % pre-stimulus length or baselin
    yBl     = ypre(1) - xpre(1) + 1;        
    blIdx   = [1:yBl];                        % baseline indexes

    csdRaS  = zeros(nBump+1,size(pair,1),nRep,nTr);
    csdRaT  = zeros(nBump,size(pair,1),nRep,nTr);
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
        for bu=1:nBump+1,
            % indexes of each stage
            meanSt  = round(mean(expLo(bu:bu+1))); 
            stagIdx = [meanSt-aPnts:meanSt+aPnts];
            % weights Transition PDF at stage points
            if bu ~= nBump+1,
                tranIdx = [expLoS(bu)-aPnts:expLoS(bu)+aPnts];
                % weights Transition PDF at stage points
                tranIdx = tranIdx + yBl;
            end
            stagIdx = stagIdx + yBl;  % correct by samples in pre-stimulus
            % Do surrogated data connectivity matrix 
            if doSurro,
                for r = 1:nRep
                    % genrate surrogate data

                    % takes a random point for each channel
%                     for ch = 1:nCh,
                        % generate random time shift for all trials
                        idrS = randi([yBl,yTr+yBl],nCh,1);
                        idrT = randi([yBl,yTr+yBl],nCh,1);

                        dataHraS = dataH(idrS);
                        % random shifting for each trial
                        dataHraT = dataH(idrT);
%                     end
                    %for pr = 1:size(pair,1)
                        % cross-spectrall density between random points of each electrode
                        csdRaS(bu,:,r,tr)   = dataH(stagIdx(aPnts+1),pair(:,1)) .* conj(dataHraS(pair(:,2))');
                        if bu ~= nBump+1,
                            csdRaT(bu,:,r,tr)	= dataH(tranIdx(aPnts+1),pair(:,1)) .* conj(dataHraT(pair(:,2))');
                        end
                    %end
                end
            end
        end
    end
    % Remove trials from current subject
    eventprobs(:,1:nTr,:) = [];
    % Connectivity
    
    % find connectivity thresholds from surrogate data 
    if doSurro,
        % connectivity
%         ispcRa  = squeeze(abs(mean(exp(1i*angle(csdRa)),3)));
        
        imagsum      = sum(imag(csdRaS),3);
        imagsumW     = sum(abs(imag(csdRaS)),3);
        debiasfactor = sum(imag(csdRaS).^2,3);
        dwpliRaS      = squeeze((imagsum.^2 - debiasfactor)./(imagsumW.^2 - debiasfactor));
        
        imagsum      = sum(imag(csdRaT),3);
        imagsumW     = sum(abs(imag(csdRaT)),3);
        debiasfactor = sum(imag(csdRaT).^2,3);
        dwpliRaT    = squeeze((imagsum.^2 - debiasfactor)./(imagsumW.^2 - debiasfactor));
        
%         wpliRa     = abs(mean(abs(imag(csdRa)).*sign(imag(csdRa)) ,3) )./mean(abs(imag(csdRa)),3);
        % Set threshold as the 95% conficence interval
%         sem = std(ispcRa,[],2)./sqrt(nRep);               % Standard Error
%         ts  = tinv(0.95, nRep-1);                          % T-Score
%         ispcCI  = mean(ispcRa,2) + ts .* sem;             % Confidence Intervals
        
        sem = std(dwpliRaS,[],3)./sqrt(nRep);               % Standard Error
        ts  = tinv(1-0.05/(nCh*(nCh-1)/2), nRep-1);                          % T-Score
        dwpliCIs  = squeeze(mean(dwpliRaS,3) + ts .* sem);            % Confidence Intervals
        
        sem = std(dwpliRaT,[],3)./sqrt(nRep);               % Standard Error
        ts  = tinv(1-0.05/((2*nBump+1)*nCh*(nCh-1)/2), nRep-1);                          % T-Score
        dwpliCIt = squeeze(mean(dwpliRaT,3) + ts .* sem);            % Confidence Intervals
        
%         sem = std(wpliRa,[],2)./sqrt(nRep);               % Standard Error
%         ts  = tinv(0.95, nRep-1);                          % T-Score
%         wpliCI  = mean(wpliRa,2) + ts .* sem;            % Confidence Intervals
    end    
        
    idSj = regexp(posnames(sj).name, '_');
    namesave = [posnames(sj).name(1:idSj) 'dwpliRandTheta.mat'];
    save(namesave,'dwpliCIt','dwpliCIs','info','freqBand')
      
%     save([savepath namesave], 'ispcPre','ispcStg','ispcTra','dwpliPre','dwpliStg','ispcGVpistg','ispcGVpitra', ...
%           'dwpliTra','ispcGV0stg','ispcGV0tra','ispcCI','dwpliCI','info','freqBand','wpdTra','wpdStg', ...
%           'ispcStgArg','ispcTraArg','wpliPre','wpliTra','wpliStg','wpliCI','dwpliPreStd','wpliPreStd','ispcPreStd')  
    clear dwpliCIs dwpliCIt sem ts dataPos dataPre debiasfactor imagsumW imagsum dataH dwpliRaS dwpliRaT csdRaS csdRaT
%     clear ispcCI ispcRa ispcRa ispcStgArg ispcTraArg ispcGVpitra ispcGVpistg expArg ispcPre ispcStg ispcTra ispcGV0stg ispcGV0tra 
%     clear wpliPre wpliTra wpliStg wpliCI wpliPreStd ispcPreStd wpliRa
%     disp(['-- Subject ', num2str(sj),' has finished. --'])
end