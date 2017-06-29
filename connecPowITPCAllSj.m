clear all
close all
clc

nTrm = 0;                           % mean trials across subjects
plotNtw = 1;                         % whether to plot network topology
nSj = 20;                           % number of subjects
nBump = 5;                          % number of bumps in HSMM
nStg = nBump + 1;                   % number of stages in HSMM
nCh = 32;                           % number of channels
fs = 100;                           % sampling frequency [Hz]
aPnts = 2;                        % number of points around the middle points of a stage to compute connectivty
nBand = 3;                          % number of frequency bands to analyze
pairs = combnk2([1:nCh],2);        % all pairs of electrodes combinations [nComb, 2]
fBand = {'Theta','Alpha','Beta'};   % frequency band to be analyzed
evenS = {'S1','S2','S3','s4','S5','S6'}; % cognitive events
evenT = {'onS2','onS3','onS4','onS5','onS6'}; % cognitive events
maxZo = 0;
minZo = 0;

% limits of each stage [ini, end]
stagLim = zeros(nBump+1,2);
stagLim(:,1) = 1:2*aPnts+1:(2*aPnts+1)*nBump+1;
stagLim(:,2) = 2*aPnts+1:2*aPnts+1:(2*aPnts+1)*(nBump+1);

pathdata = '/home/oscar/Documents/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/PhaseConectiv/';
pathporb = '/home/oscar/Documents/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/';
savepath = '/home/oscar/Documents/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/PhaseConectiv/';
namePDF = ['AllSj100HSMM' num2str(nBump) 'Bout235.mat'];
info = {'Computes thresholded connectivity matrices (dwPLI), power and ITPC',...
        'By subjects and globally, Fs = 100Hz, . Theta, Alpha and Beta bands'};

% test that the subjects are sorted equaly when loading so loaded subjects agree
pownames = dir([pathdata '*itpcPowFre' fBand{1} '.mat']);
dwplinames = dir([pathdata '*dwpliFre' fBand{1} '.mat']);
if length(pownames) ~= length(dwplinames), error('Unbalancend number of files per subject'), end
for sj = 1:nSj
    idTest = strfind(dwplinames(sj).name ,'_');
    idpos = strfind(pownames(sj).name ,'_');
    if pownames(sj).name(1:idpos) ~= dwplinames(sj).name(1:idTest)
        error('Error: Subjects are not the same')
    end
end
clear idpos idTest

conAlTra = zeros(size(pairs,1),nBump,nSj,nBand);
conAlStg = zeros(size(pairs,1),nStg,nSj,nBand);
powAlTra = zeros(nCh,nBump,nSj,nBand);
powAlStg = zeros(nCh,nStg,nSj,nBand);
powInAlTra = zeros(nCh,nBump,nSj,nBand);
powInAlStg = zeros(nCh,nStg,nSj,nBand);
itpcAlTra = zeros(nCh,nBump,nSj,nBand);
itpcAlStg = zeros(nCh,nStg,nSj,nBand);

load onsetStgTimeInfo5Bu.mat
load([pathdata 'electrodes.mat'])

for fr = 1:nBand
    pownames = dir([pathdata '*itpcPowFre' fBand{fr} '.mat']);
    dwplinames = dir([pathdata '*dwpliFre' fBand{fr} '.mat']);
    for sj = 1:nSj
        % Connectivity analysis
        load([pathdata dwplinames(sj).name],'dwpliCIs', 'dwpliCIt', 'dwpliPreMed', 'dwpliStg','dwpliTra')
        % variavility on the connectivity window
        dwpliTraStd    = squeeze(std(dwpliTra));
        dwpliStgStd    = squeeze(std(dwpliStg));
        % median connectivity across network window
        dwpliStg    = squeeze(median(dwpliStg,1));
        dwpliTra    = squeeze(median(dwpliTra,1));
        % remove baseline & unceratanity on connectivity network
        dwpliStg    = dwpliStg - repmat(dwpliPreMed,nStg,1)' - dwpliStgStd;
        dwpliTra    = dwpliTra - repmat(dwpliPreMed,nBump,1)' - dwpliTraStd;    
        % Remove connectivity below confidence Interval
        noiseBL = dwpliStg <= dwpliCIs';
        dwpliStg(noiseBL) = 0;
        noiseBL = dwpliTra <= dwpliCIt';
        dwpliTra(noiseBL) = 0;
        % assign values to all subjects matrices
        conAlTra(:,:,sj,fr) = dwpliTra;
        conAlStg(:,:,sj,fr) = dwpliStg;
        clear dwpliStg dwpliTra noiseBL dwpliCIs dwpliCIt dwpliPreMed dwpliPreStd dwpliPre dwpliTraStd dwpliStgStd
        
        % Power and ITPC analysis
        load([pathdata pownames(sj).name])
%         % baseline 
%         itpcPre = median(itpcPre);
%         powPre = median(powPre);
%         % baseline correction at cognitive events (median)
         itpcAlStg(:,:,sj,fr) = itpcStg;
         itpcAlTra(:,:,sj,fr) = itpcTra;
         powAlStg(:,:,sj,fr) = powStg;
         powAlTra(:,:,sj,fr) = powTra;
         powInAlStg(:,:,sj,fr) = powInStg;
         powInAlTra(:,:,sj,fr) = powInTra;
         nTrm = nTrm + nTr;    % counts number of trial in all subjects
         clear itpcStg itpcTra powStg powTra powInTra powInStg nTr
    end
    % critical ITPCC
    nTrm = nTrm / nSj;  % mean number of trials across subjects
    crItpc = sqrt(-log(0.05/(nCh*11*2))/nTrm);
    % Global power and ITPC
    gPowStg = squeeze(mean(powAlStg,3));
    gPowTra = squeeze(mean(powAlTra,3));
    gPowInStg = squeeze(mean(powInAlStg,3));
    gPowInTra = squeeze(mean(powInAlTra,3));
    gItpcStg = squeeze(mean(itpcAlStg,3));
    gItpcTra = squeeze(mean(itpcAlTra,3));
    % Mean stage onset and mean point between onsets[ms]
    %timeT = [0, round(10 * cumsum([dur(21).stag1, dur(21).stag2, dur(21).stag3, dur(21).stag4, dur(21).stag5, dur(21).stag6]))];
    timeT = round([dur(1).stag1, dur(1).stag2, dur(1).stag3, dur(1).stag4, dur(1).stag5]);
    timeS = round(movmean([0, timeT, dur(1).avDur],2,'Endpoints','discard'));
    %timeT = timeT(2:end-1);
end

pltS = gPowInStg;
pltT = gPowInTra;

maxS = squeeze(max(max(max(pltS))));
minS = squeeze(min(min(min(pltS))));
maxT = squeeze(max(max(max(pltT))));
minT = squeeze(min(min(min(pltT))));
minZ = min([minS,minT]);
maxZ = max([maxS,maxT]);
[~, iMax ] = max([abs(minZ),maxZ]);
if iMax == 1
    maxZ = abs(minZ);
else
    minZ = - maxZ;
end

% global maximum & minimum connectivity strength for ploting
maxC = 0;
minC = 1;
for fr = 1:nBand
    for bu = 1:nStg
        gPairS = [];
        gStrenS = [];
        gPairT = [];
        gStrenT = [];
        thCo = 14;  % threshold connectivity consistency, best 14
        for pr = 1:size(pairs,1)
            if sum(logical(conAlStg(pr,bu,:,fr))) >= thCo
                gPairS(end+1,:) = pairs(pr,:);
                gStrenS(end+1) = median(conAlStg(pr,bu,logical(conAlStg(pr,bu,:,fr)),fr));
            end
            if bu ~= nStg
                if sum(logical(conAlTra(pr,bu,:,fr))) >= thCo
                    gPairT(end+1,:) = pairs(pr,:);
                    gStrenT(end+1) = median(conAlTra(pr,bu,logical(conAlTra(pr,bu,:,fr)),fr));
                end
            end
        end
        maxC = max([gStrenS,gStrenT,maxC]);
        minC = min([gStrenS,gStrenT,minC]);
    end
end

for fr = 1:nBand
      for bu = 1:nStg
        gPairS = [];
        gStrenS = [];
        gPairT = [];
        gStrenT = [];
        thCo = 14;  % threshold connectivity consistency, best 14
        for pr = 1:size(pairs,1)
            if sum(logical(conAlStg(pr,bu,:,fr))) >= thCo
                gPairS(end+1,:) = pairs(pr,:);
                gStrenS(end+1) = median(conAlStg(pr,bu,logical(conAlStg(pr,bu,:,fr)),fr));
            end
            if bu ~= nStg
                if sum(logical(conAlTra(pr,bu,:,fr))) >= thCo
                    gPairT(end+1,:) = pairs(pr,:);
                    gStrenT(end+1) = median(conAlTra(pr,bu,logical(conAlTra(pr,bu,:,fr)),fr));
                end
            end
        end
        % Global power and ITPC
        if plotNtw
            pltS = gPowInStg;
            pltT = gPowInTra;
            tit = [fBand{fr} ',Power'];
            
            figure(fr)
            
            colMap = colormap('jet');
%             colMap(1,:) = [1, 0, 1];
            
%             maxS = squeeze(max(max(pltS)));
%             minS = squeeze(min(min(pltS)));
%             maxT = squeeze(max(max(pltT)));
%             minT = squeeze(min(min(pltT)));
%             minZ = min([minS,minT]');
%             maxZ = max([maxS,maxT]');
            
%             minZ = [crItpc,crItpc,crItpc];
%             maxZ = max([maxZ;minZ]);
            
            % stages
            ds.chanPairs = gPairS;
            ds.connectStrength = gStrenS;
            ds.connectStrengthLimits = [0, maxC];
            subplot(2,6,bu)
            topoplot_connect(ds, locElect)
            topoplot(pltS(:,bu,fr),locElect,'maplimits',[minZ,maxZ], 'colormap', colMap,'hcolor','none','style','map','electrodes','off');
%             alpha(0.7)
            title([evenS{bu} ': ' num2str(timeS(bu)) ' ms'], 'FontSize', 8)
            if bu ~= nStg %& gPairT
                % transitions            
                ds.chanPairs = gPairT;
                ds.connectStrength = gStrenT;
                ds.connectStrengthLimits = [0, maxC];
                subplot(2,6,bu+nStg)
                topoplot_connect(ds, locElect)
                topoplot(pltT(:,bu,fr),locElect,'maplimits',[minZ,maxZ],'colormap', colMap,'hcolor','none','style','map','electrodes','off');
%                 alpha(0.7)
                title([evenT{bu} ': ' num2str(timeT(bu)) ' ms'], 'FontSize', 8)
            else
                subplot(2,6,bu+nStg)
                title(tit, 'FontSize', 8)
                cb = colorbar('Ticks',[0,0.5,1],'TickLabels',{num2str(minZ), '0', num2str(maxZ)});
                % Connectivity colorbar

%                 annotation('textbox',...
%                     cb.Position,...
%                     'FitBoxToText','on',...
%                     'FaceAlpha',0.5,...
%                     'EdgeColor',[1 1 1],...
%                     'BackgroundColor',[1 1 1]);
            end
        end
    end
end
        