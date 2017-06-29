% It time-locks ohase consistency across teh whoe scalp at the expected bumps location. There
% is an interpolation and extraolation of trials to average all trials when
% time-locked to the location of bumps

clear all
close all
clc

nBump = 5;                          % number of bumps in HSMM
nCh = 32;                           % number of channels
freqBand = 'theta';                 % frequency band to be analyzed
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

info = ['fs = 100, Band: theta. Intersite Phase consistency (ITPC across the sacal at each sample) of the' ...
        ' analytic signal are time locked to bumps expected location.'];

% load variables
load([pathdata 'electrodes.mat'])
load([pathdata namePDF], 'eventprobs')


% lockPot = {};
meanPhaCons = zeros(avDur, nSj);

for sj = 1:nSj
    % load datasets
    data = load([pathdata namesHilb(sj).name], freqBand, 'x', 'y');
    x = data.x;
    y = data.y;
    % do pahse angle from analytic signal
    dataPha = angle(data.(freqBand));
    dataRe = real(data.(freqBand));
    clear data
    
    dataCons = zeros(avDur,length(y));
    % Concatenate trials vertically to remove zeros out of trial
    for tr = 1:length(y)
        yTr = y(tr) - x(tr) + 1;
       
        trPha = dataPha(x(tr):y(tr),:);
        trRe = dataRe(x(tr):y(tr),:);
        % Phase consistency across the scalp
        dataConstr = squeeze(abs(mean(exp(1i*trPha),2)));
        % Expected location of bumps
        expLoc = round([1:yTr] * reshape(eventprobs(1:yTr,tr,:),yTr,nBump));
        % boundaries of a trial
        expLoc = [1 expLoc yTr-1];
        trialCons = [];
        for st = 1:length(avStag),
            % Do new samples from originals samples
            ns = linspace(expLoc(st), expLoc(st+1) - 1, avStag(st));
            % Interpolate data to new samples
            nd = interp1(expLoc(st):expLoc(st+1)-1, dataConstr(expLoc(st):expLoc(st+1)-1,:) ,ns, 'spline', 'extrap');
            trialCons = vertcat(trialCons, nd');
      
        end
        dataCons(:,tr) = trialCons;
    end
    %clear dataPha dataRe instFr trialInst trialPha trialUpha trialRe trialPot

    meanPhaCons(:,sj) = mean(dataCons,2);
    
    % Remove trials from current subject
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --------------- ATTENTION!!!
    eventprobs(:,1:length(y),:) = [];
    %clear eventprobs % Only for the case of ONE Subject
end  

pathasave = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/histPhase/';
% % save([pathasave 'bumpLocking.mat'], 'info','lockPha','lockUpha','lockIfre','lockRe', ...
% %     'meanPha','meanUpha','meanIfre','meanRe','gfInstF','gfPha', ...
% %     'gfUpha','gfRe','lockPot','gfPot','meanPot')
save([pathasave 'bumpLockingIallSPC.mat'], 'info', ...
    'meanPhaCons','avStag')