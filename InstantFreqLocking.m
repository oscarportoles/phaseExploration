% It time-locks instantaneous phase to the expected bumps location. There
% is an interpolation and extraolation of trials to average all trials when
% time-locked to the location of bumps

clear all
close all
clc

nBump = 5;                          % number of bumps in HSMM
nCh = 32;                           % number of channels
freqBand = 'theta';                 % frequency band to be analyzed
nRep = 400;                         % number of repetition to build a random PDF
avStag = [11, 12,19, 46, 16, 12];   % average duration of each stage on samples
avDur = sum(avStag);                % average duration of all trials in samples: sum
fs = 100;                           % sampling frequency [Hz]              
% n_order = 10;                       % median filter parameters (Cohen X et.al. FReq sliding...)
% orders = round(linspace(0.02*fs,0.4*fs,n_order)); % recommended: 10 steps between 10 (here 20ms becoause of fs) and 400 ms
% orders = floor((orders-1)/2);         % pre/post halves

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

info = ['fs = 100, Band: theta. phase, angele, envelop and real part of the' ...
        ' analytic signal are time locked to bumps expected location. Instantaneous ferquency is not filtered'];

% load variables
load([pathdata 'electrodes.mat'])
load([pathdata namePDF], 'eventprobs')

% lockPha = {};
% lockUpha = {};
% lockIfre = {};
% lockRe = {};
% % lockEnv = {};
% lockPot = {};
meanPha = zeros(avDur, nCh, nSj);
meanUpha = zeros(avDur, nCh, nSj);
meanIfre = zeros(avDur, nCh, nSj);
meanRe = zeros(avDur, nCh, nSj);
% meanEnv = zeros(avDur, nCh, nSj);
meanPot = zeros(avDur, nCh, nSj);

meanGFiFr = zeros(avDur,nSj);
meanGFpha = zeros(avDur,nSj);
meanGFre = zeros(avDur,nSj);
meanGFpot = zeros(avDur,nSj);

for sj = 1:nSj
    % load datasets
    data = load([pathdata namesHilb(sj).name], freqBand, 'x', 'y');
    x = data.x;
    y = data.y;
    % do pahse angle from analytic signal
    dataPha = angle(data.(freqBand));
    dataRe = real(data.(freqBand));
%     dataEnv = abs(data.(freqBand));
    clear data
    
    instFrN = zeros(avDur, nCh, length(y));
    dataPhaN = zeros(avDur, nCh, length(y));
    dataReN = zeros(avDur, nCh, length(y));
%     dataEnvN = zeros(avDur, nCh, length(y));
    dataUphaN = zeros(avDur, nCh, length(y));
    dataPotN = zeros(avDur, nCh, length(y));
    
    dataGFiFr = zeros(avDur,length(y));
    dataGFpha = zeros(avDur,length(y));
    dataGFre = zeros(avDur,length(y));
    dataGFpot = zeros(avDur,length(y));
    % Concatenate trials vertically to remove zeros out of trial
    for tr = 1:length(y)
        yTr = y(tr) - x(tr) + 1;
        % Instantaneous frequency
        instFr = fs .* diff(unwrap(dataPha(x(tr):y(tr),:))) ./ (2*pi);
        instFr = vertcat(instFr, instFr(end,:)); % duplicate last sumple to compesate 'diff'
        % Add mirrowed tails at the edges of the trial to avoid edge artifacts
%         sMirr = orders(end) + 1;        % number of samples to mirrow
%         instFr = [flipud(instFr(1:sMirr,:)); instFr; flipud(instFr(end-sMirr:end,:))];
        % Compute median filtered instantaneous (Cohen X 2014 Freq sliding... 2014)
        % mean of the median of length varying sliding windows
%         phasedmed = zeros(length(orders),size(instFr,1),nCh);
%         for oi=1:n_order
%             for ti=1:size(instFr,1)
%                 for ch =1:nCh
%                     temp = sort(instFr(max(ti-orders(oi),1):min(ti+orders(oi),size(instFr,1)),ch));
%                     phasedmed(oi,ti,ch) = temp(floor(numel(temp)/2)+1);
%                 end
%             end
%         end
%         instFr = fs .* squeeze(mean(phasedmed(:,sMirr+1:end-sMirr-1,:),1)) ./ (2*pi);
        % end median filter
        trPha = dataPha(x(tr):y(tr),:);
        trRe = dataRe(x(tr):y(tr),:);
%         trEnv = dataEnv(x(tr):y(tr),:);
        trPot = dataRe(x(tr):y(tr),:).^2;
        
        % VAriability between all electrodes at each sample
        trGFpha = std(detrend(unwrap(dataPha(x(tr):y(tr),:))),[],2);    
        % detren unwrap to avoid variability arifacts in jumps from pi to -pi and detrend to correct long-term ofset
        trGFre = std(trRe,[],2);
        trGFpot = std(trPot,[],2);
        trGFiFr = std(instFr,[],2);
        % Expected location of bumps
        expLoc = round([1:yTr] * reshape(eventprobs(1:yTr,tr,:),yTr,nBump));
        % boundaries of a trial
        expLoc = [1 expLoc yTr-1];
        % shift oposite pahses to unique direction
        trUpha = trPha;
        trUpha(trUpha < 0) = trUpha(trUpha < 0) + pi;
        if ~isempty(find(trUpha < 0)), error('Not all phaes are above zero'), end
        
        trialInst = [];
        trialPha = [];
        trialUpha = [];
        trialRe = [];
%         trialEnv = [];
        trialPot = [];
        trialGFre = [];
        trialGFiFr = [];
        trialGFpot = [];
        trialGFpha = [];
        for st = 1:length(avStag),
            % Do new samples from originals samples
            ns = linspace(expLoc(st), expLoc(st+1) - 1, avStag(st));
            % Interpolate data to new samples
            nd = interp1(expLoc(st):expLoc(st+1)-1, instFr(expLoc(st):expLoc(st+1)-1,:) ,ns, 'spline', 'extrap');
            trialInst = vertcat(trialInst, nd);
            nd = interp1(expLoc(st):expLoc(st+1)-1, trGFiFr(expLoc(st):expLoc(st+1)-1,:) ,ns, 'spline', 'extrap');
            trialGFiFr = vertcat(trialGFiFr, nd');
       
            nd = interp1(expLoc(st):expLoc(st+1)-1, trPha(expLoc(st):expLoc(st+1)-1,:) ,ns, 'spline', 'extrap');
            trialPha = vertcat(trialPha, nd);            
            nd = interp1(expLoc(st):expLoc(st+1)-1, trGFpha(expLoc(st):expLoc(st+1)-1,:) ,ns, 'spline', 'extrap');
            trialGFpha = vertcat(trialGFpha, nd');  
            
            nd = interp1(expLoc(st):expLoc(st+1)-1, trRe(expLoc(st):expLoc(st+1)-1,:) ,ns, 'spline', 'extrap');
            trialRe = vertcat(trialRe, nd);  
            nd = interp1(expLoc(st):expLoc(st+1)-1, trGFre(expLoc(st):expLoc(st+1)-1,:) ,ns, 'spline', 'extrap');
            trialGFre = vertcat(trialGFre, nd'); 
            
%             ndEnv = interp1(expLoc(st):expLoc(st+1)-1, trEnv(expLoc(st):expLoc(st+1)-1,:) ,ns, 'spline', 'extrap');
%             trialEnv = vertcat(trialEnv, ndEnv);
            
            nd = interp1(expLoc(st):expLoc(st+1)-1, trUpha(expLoc(st):expLoc(st+1)-1,:) ,ns, 'spline', 'extrap');
            trialUpha = vertcat(trialUpha, nd);   
            
            nd = interp1(expLoc(st):expLoc(st+1)-1, trPot(expLoc(st):expLoc(st+1)-1,:) ,ns, 'spline', 'extrap');
            trialPot = vertcat(trialPot, nd);
            nd = interp1(expLoc(st):expLoc(st+1)-1, trGFpot(expLoc(st):expLoc(st+1)-1,:) ,ns, 'spline', 'extrap');
            trialGFpot = vertcat(trialGFpot, nd');            
        end
        instFrN(:,:,tr) = trialInst;
        dataPhaN(:,:,tr) = trialPha;
        dataUphaN(:,:,tr) = trialUpha;
        dataReN(:,:,tr) = trialRe;
%         dataEnvN(:,:,tr) = trialEnv;
        dataPotN(:,:,tr) = trialPot;
        
        dataGFiFr(:,tr) = trialGFiFr;
        dataGFpha(:,tr) = trialGFpha;
        dataGFre(:,tr) = trialGFre;
        dataGFpot(:,tr) = trialGFpot;
 
    end
    clear dataPha dataRe instFr trialInst trialPha trialUpha trialRe trialPot
    meanPha(:,:,sj) = mean(dataPhaN,3);
    meanUpha(:,:,sj) = mean(dataUphaN,3);
    meanIfre(:,:,sj) = mean(instFrN,3);
    meanRe(:,:,sj) = mean(dataReN,3);
%     meanEnv(:,:,sj) = mean(dataEnvN,3);
    meanPot(:,:,sj) = mean(dataPotN,3);

    meanGFiFr(:,sj) = mean(dataGFiFr,2);
    meanGFpha(:,sj) = mean(dataGFpha,2);
    meanGFre(:,sj) = mean(dataGFre,2);
%     gfEnv(:,sj) = mean(std(dataEnvN,[],2),3);
    meanGFpot(:,sj) = mean(dataGFpot,2);
    
%     lockPha{sj} = dataPhaN;
%     lockUpha{sj} = dataUphaN;
%     lockIfre{sj} = instFrN;
%     lockRe{sj} = dataReN;
% %     lockEnv{sj} = dataEnvN; 
%     lockPot{sj} = dataPotN;   
    
    % Remove trials from current subject
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --------------- ATTENTION!!!
    eventprobs(:,1:length(y),:) = [];
    %clear eventprobs % Only for the case of ONE Subject
    clear dataPhaN dataReN dataUphaN dataPotN instFrN
end  

pathasave = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/histPhase/';
% save([pathasave 'bumpLocking.mat'], 'info','lockPha','lockUpha','lockIfre','lockRe', ...
%     'meanPha','meanUpha','meanIfre','meanRe','gfInstF','gfPha', ...
%     'gfUpha','gfRe','lockPot','gfPot','meanPot')
save([pathasave 'bumpLocking.mat'], 'info', ...
    'meanPha','meanUpha','meanIfre','meanRe','meanGFiFr','meanGFpha', ...
    'meanGFre','meanGFpot','meanPot')