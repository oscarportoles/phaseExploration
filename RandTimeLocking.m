% It time-locks to random bump locations time series like phase, power,
% amplitde and scapl variavility of pahse, power and amplitude. The time
% locking points and the expected location of a bump are random

clear all
close all
clc

nBump = 5;                          % number of bumps in HSMM
nCh = 32;                           % number of channels
freqBand = 'theta';                 % frequency band to be analyzed
nRep = 5;                          % number of repetition to build a random PDF
avDur = 116;                        % average duration of all trials in samples. Same in empirical data for the sake of comparisson
tLock = zeros(nBump+1,nRep);          % RAndom time locking points for all subjects
fs = 100;                           % sampling frequency [Hz]
cr = 6;                             % constant. minimal distant between stages (max = 6)


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

info = ['fs = 100, Band: theta. angele, Power and real part of the' ...
        ' analytic signal are time locked to bumps with random expected location'];

% load variables
load([pathdata 'electrodes.mat'])

meanPha = zeros(avDur, nCh, nSj, nRep);
meanUpha = zeros(avDur, nCh, nSj, nRep);
meanRe = zeros(avDur, nCh, nSj, nRep);
meanPot = zeros(avDur, nCh, nSj, nRep);

meanGFre = zeros(avDur,nSj, nRep);
meanGFpot = zeros(avDur,nSj, nRep);
meanGFpha = zeros(avDur,nSj, nRep);

for re = 1:nRep,
    load([pathdata namePDF], 'eventprobs')
    % Random time locking points
    avStag = zeros(1,nBump+1);
    first = 1;
    while avStag(6) <= cr | sum(avStag) ~= avDur | first
        avStag = zeros(nBump+1,1);           % Random time locking points for one subjecte
        avStag(1) = randi([cr, round(avDur/4)]);
        avStag(2) = randi([cr, round(avDur/4)]);
        avStag(3) = randi([cr, round(avDur/4)]);
        avStag(4) = randi([cr, round(avDur/3)]);
        avStag(5) = randi([cr, round(avDur/4)]);
        avStag(6) = avDur - sum(avStag);
        first = 0;
    end
    tLock(:,re) = avStag;      % samples where data time-loced at.
    
    for sj = 1:nSj
        % load datasets
        data = load([pathdata namesHilb(sj).name], freqBand, 'x', 'y');
        x = data.x;
        y = data.y;
        % do pahse angle from analytic signal
        dataPha = angle(data.(freqBand));
        dataRe = real(data.(freqBand));  
        clear data

        dataPhaN = zeros(avDur, nCh, length(y));
        dataReN = zeros(avDur, nCh, length(y));
        dataUphaN = zeros(avDur, nCh, length(y));
        dataPotN = zeros(avDur, nCh, length(y));
        
        dataGFpha = zeros(avDur,length(y));
        dataGFre = zeros(avDur,length(y));
        dataGFpot = zeros(avDur,length(y));
        % Concatenate trials vertically to remove zeros out of trial
        for tr = 1:length(y)
            yTr = y(tr) - x(tr) + 1;
            trPha = dataPha(x(tr):y(tr),:);
            trRe = dataRe(x(tr):y(tr),:);
            trPot = dataRe(x(tr):y(tr),:).^2;
            % shift oposite phases to unique direction
            trUpha = trPha;
            trUpha(trUpha < 0) = trUpha(trUpha < 0) + pi;
            % global variavility
            % VAriability between all electrodes at each sample
            trGFpha = std(detrend(unwrap(dataPha(x(tr):y(tr),:))),[],2);    
            % detren unwrap to avoid variability arifacts in jumps from pi to -pi and detrend to correct long-term ofset
            trGFre = std(trRe,[],2);
            trGFpot = std(trPot,[],2);
            if ~isempty(find(trUpha < 0)), error('Not all phaes are above zero'), end
            % Random Expected location of bumps          
            expLoc = zeros(1, nBump+1);
            first = 1;
            while expLoc(6) <= cr | sum(expLoc) ~= yTr | first
                expLoc = zeros(1,nBump+1);           
                expLoc(1) = randi([cr, round(avDur/4)]);
                expLoc(2) = randi([cr, round(avDur/4)]);
                expLoc(3) = randi([cr, round(avDur/4)]);
                expLoc(4) = randi([cr, round(avDur/3)]);
                expLoc(5) = randi([cr, round(avDur/4)]);
                expLoc(6) = yTr - sum(expLoc);
                first = 0;
            end
            expLoc = [1, cumsum(expLoc)]; % Trial boundary
            trialPha = [];
            trialUpha = [];
            trialRe = [];
            trialPot = [];
            
            trialGFre = [];
            trialGFpot = [];
            trialGFpha = [];
            for st = 1:length(avStag),
                % Do new samples from originals samples
                ns = linspace(expLoc(st), expLoc(st+1) - 1, avStag(st));
                % Interpolate data to new samples
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
            dataPhaN(:,:,tr) = trialPha;
            dataUphaN(:,:,tr) = trialUpha;
            dataReN(:,:,tr) = trialRe;
    %         dataEnvN(:,:,tr) = trialEnv;
            dataPotN(:,:,tr) = trialPot;

            dataGFpha(:,tr) = trialGFpha;
            dataGFre(:,tr) = trialGFre;
            dataGFpot(:,tr) = trialGFpot;

        end
        clear dataPha dataRe trialPha trialUpha trialRe trialPot
        meanPha(:,:,sj,re) = mean(dataPhaN,3);
        meanUpha(:,:,sj,re) = mean(dataUphaN,3);
        meanRe(:,:,sj,re) = mean(dataReN,3);
    %     meanEnv(:,:,sj) = mean(dataEnvN,3);
        meanPot(:,:,sj,re) = mean(dataPotN,3);

        meanGFpha(:,sj,re) = mean(dataGFpha,2);
        meanGFre(:,sj,re) = mean(dataGFre,2);
    %     gfEnv(:,sj) = mean(std(dataEnvN,[],2),3);
        meanGFpot(:,sj,re) = mean(dataGFpot,2);

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
        clear dataPhaN dataReN dataUphaN dataPotN 
    end  
end


pathasave = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/histPhase/';
save([pathasave 'RandBumpLocking.mat'], 'info', ...
    'meanPha','meanUpha','meanRe','meanGFpha', ...
    'meanGFre','meanGFpot','meanPot')