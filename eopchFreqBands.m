% It epochs band pass filtered data. The epochs are based on the 'bumps HSMM'
% segments. So that it takes exactly the same segments. NEw data is
% separeted by conditions and freuency bands.

clear all
close all
clc

% Full spectrum data set
pathdataO = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125bp05_30/events/';
snamesO = dir([pathdataO '*.set']);
NsjO = size(snamesO,1);
% Delta band data
pathdataD = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125_EEGbands/Delta/';
snamesD = dir([pathdataD '*.set']);
% Theta band data
pathdataT = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125_EEGbands/Theta/';
snamesT = dir([pathdataT '*.set']);
% Alpha band data
pathdataA = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125_EEGbands/Alpha/';
snamesA = dir([pathdataA '*.set']);
% Beta band data
pathdataB = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125_EEGbands/Beta/';
snamesB = dir([pathdataB '*.set']);
% Gamma band data
pathdataG = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125_EEGbands/Gamma/';
snamesG = dir([pathdataG '*.set']);


% Order subject files by subject number
subjctNO = zeros(NsjO,1);
for sj=1:NsjO,
    number = strsplit(snamesO(sj).name(2:3), '_');
    subjctNO(sj) = str2num(number{1});
    number = strsplit(snamesD(sj).name(2:3), '_');
    subjctND(sj) = str2num(number{1});
    number = strsplit(snamesT(sj).name(2:3), '_');
    subjctNT(sj) = str2num(number{1});
    number = strsplit(snamesA(sj).name(2:3), '_');
    subjctNA(sj) = str2num(number{1});
    number = strsplit(snamesB(sj).name(2:3), '_');
    subjctNB(sj) = str2num(number{1});
    number = strsplit(snamesG(sj).name(2:3), '_');
    subjctNG(sj) = str2num(number{1});
end
[~, Isort] = sort(subjctNO);
snamesO = snamesO(Isort);
[~, Isort] = sort(subjctND);
snamesD = snamesD(Isort);
[~, Isort] = sort(subjctNT);
snamesT = snamesT(Isort);
[~, Isort] = sort(subjctNA);
snamesA = snamesA(Isort);
[~, Isort] = sort(subjctNB);
snamesB = snamesB(Isort);
[~, Isort] = sort(subjctNG);
snamesG = snamesG(Isort);

% Check that all frequency bands have the same number of subjects and all
% subjects are the same
if (length(snamesO) == length(snamesD) == length(snamesT) == ...
    length(snamesA) == length(snamesB) == length(snamesG)) | ...
    (sum(subjctNO) == sum(subjctND) == sum(subjctNT) == ...
    sum(subjctNA) == sum(subjctNB) == sum(subjctNG))
    error('Not all frequency bands have the same subjects')
end


% SAve data to path
pathsave = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125bp05_58/events/';
namesave = 'SegbyBand125_0558.mat';
saveTo = [pathsave namesave];

info = ['segment bandpass filtered data like if it were for a bumps HSMM;' ...
    '125 Hz sampling freq.; filter 0.1 to 58 Hz;' ...
    'lenght trials[Onset, REsp], conditions: TargetFan1, TagetFan2, FoilFan1, FoilFan2, NewFan1;' ...
    'variables [npts,Ch]; Delta band 2 - 4 Hz, Theta band 4 - 9 Hz, ' ...
    'Alpha band 9 - 14 Hz, Beta band 14 - 30 Hz, Gamma band 30 - 58 Hz'];

% % Define epoching boundaries on time
% init = -0.2;    % time [s] beginin of trial from Onset
% fint = 0.16;    % time [s] end of trial from Response

% Define epoching boundaries on time
init = 0;    % time [s] beginin of trial from Onset
fint = 0;    % time [s] end of trial from Response

trCode =   [[2 3 1]  % first column: Jelmer's code for short words
            [4 5 2]  % second column: Jelmer's code for long words
            [6 7 3]  % third column: new code
            [8 9 4]
            [10 11 5]];
deltaEEG = [];
thetaEEG = [];
alphaEEG = [];
betaEEG = [];
gammaEEG = [];
        
for sj=1:NsjO,
    % Open subject dataset
    EEG = pop_loadset('filename',snamesO(sj).name ,'filepath',pathdataO);
    EEG = eeg_checkset( EEG );
    
    EEGD = pop_loadset('filename',snamesD(sj).name ,'filepath',pathdataD);
    EEGD = eeg_checkset( EEGD );
    EEGT = pop_loadset('filename',snamesT(sj).name ,'filepath',pathdataT);
    EEGT = eeg_checkset( EEGT );
    EEGA = pop_loadset('filename',snamesA(sj).name ,'filepath',pathdataA);
    EEGA = eeg_checkset( EEGA );
    EEGB = pop_loadset('filename',snamesB(sj).name ,'filepath',pathdataB);
    EEGB = eeg_checkset( EEGB );
    EEGG = pop_loadset('filename',snamesG(sj).name ,'filepath',pathdataG);
    EEGG = eeg_checkset( EEGG );

    % get indexes to don't iterate empty events
    idxEvnt = [];
    jump = 0;
    for i=1:length(EEG.event),
        if not(jump),
            % check that there is not a 'boundary' event within a trial
            if not(i == length(EEG.event)) & strcmp(EEG.event(i).type, 'boundary') & EEG.event(i-1).trialN == EEG.event(i+1).trialN,
                idxEvnt(end) = [];
                jump = 1;
                continue
            end 
            % include only non 'boundary' or 'NotFull' events
            if ~isempty(EEG.event(i).code),
                idxEvnt(end+1)=i;
            end
        else
            jump = 0;
        end
    end
    
    % Calculates the reaction time std per condition and per
    % subject to set the thresholds of accepted trials within
    % mean +- 3 * std  
    fa1ta = [];
    fa2ta = [];
    fa1fo = [];
    fa2fo = [];
    newpa = [];
    for i=1:length(EEG.event),
        if EEG.event(i).accur,
            if EEG.event(i).code == 2 | EEG.event(i).code == 3,
                fa1ta(end+1) = EEG.event(i).RT;
            elseif EEG.event(i).code == 4 | EEG.event(i).code == 5,
                fa2ta(end+1) = EEG.event(i).RT;
            elseif EEG.event(i).code == 6 | EEG.event(i).code == 7,
                fa1fo(end+1) = EEG.event(i).RT;
            elseif EEG.event(i).code == 8 | EEG.event(i).code == 9,
                fa2fo(end+1) = EEG.event(i).RT;
            elseif EEG.event(i).code == 10 | EEG.event(i).code == 11,
                newpa(end+1) = EEG.event(i).RT;
            end
        end
    end  
    
    % remove duplicated trials
    fa1ta(1:2:end) = [];
    fa2ta(1:2:end) = [];
    fa1fo(1:2:end) = [];
    fa2fo(1:2:end) = [];
    newpa(1:2:end) = [];

    % do std + -
    fa1taSD(1) = mean(fa1ta) - 3 * std(fa1ta);
    fa2taSD(1) = mean(fa2ta) - 3 * std(fa2ta);
    fa1foSD(1) = mean(fa1fo) - 3 * std(fa1fo);
    fa2foSD(1) = mean(fa2fo) - 3 * std(fa2fo);
    newpaSD(1) = mean(newpa) - 3 * std(newpa);

    fa1taSD(2) = mean(fa1ta) + 3 * std(fa1ta);
    fa2taSD(2) = mean(fa2ta) + 3 * std(fa2ta);
    fa1foSD(2) = mean(fa1fo) + 3 * std(fa1fo);
    fa2foSD(2) = mean(fa2fo) + 3 * std(fa2fo);
    newpaSD(2) = mean(newpa) + 3 * std(newpa);        
    
    % check if the number of event is odd -> if not a event is missing
    if mod(length(idxEvnt),2), error('The number of events is not even'); end
 
    for tr = 1:2:length(idxEvnt)-1, % Iterate all events to segment
        % Check that consecutive events belong to the same trial
        if EEG.event(idxEvnt(tr)).code == EEG.event(idxEvnt(tr+1)).code, 
            % get the condition code of this trial
            [row,col] = find((trCode(:,1:2) == EEG.event(idxEvnt(tr)).code));
            trialCode = trCode(row,3);
            
            % Define epoch boundaries on samples
            ini = floor(EEG.event(idxEvnt(tr)).latency + init * EEG.srate + 1); % attention! add +1 to get the same segments ans in Jelmer's data 
            fin = round(EEG.event(idxEvnt(tr+1)).latency + fint * EEG.srate);

            % Use the std values computed for Short land long words together  
            if [EEG.event(idxEvnt(tr):idxEvnt(tr+1)).accur] & ...
                    [EEG.event(idxEvnt(tr)).RT < 3000],
                
                if (trialCode == 1 & [EEG.event(idxEvnt(tr)).RT > fa1taSD(1)] & [EEG.event(idxEvnt(tr)).RT < fa1taSD(2)]) | ...
                        (trialCode == 2 & [EEG.event(idxEvnt(tr)).RT > fa2taSD(1)] & [EEG.event(idxEvnt(tr)).RT < fa2taSD(2)]) | ...
                        (trialCode == 3 & [EEG.event(idxEvnt(tr)).RT > fa1foSD(1)] & [EEG.event(idxEvnt(tr)).RT < fa1foSD(2)]) | ...
                        (trialCode == 4 & [EEG.event(idxEvnt(tr)).RT > fa2foSD(1)] & [EEG.event(idxEvnt(tr)).RT < fa2foSD(2)]) | ...
                        (trialCode == 5 & [EEG.event(idxEvnt(tr)).RT > newpaSD(1)] & [EEG.event(idxEvnt(tr)).RT < newpaSD(2)]),
                        % check std per condition
                    % Extract data epoch, detrend it and test amplitude artefacts
                    dataTest = detrend(EEG.data(:,ini:fin)', 'linear');
                    if not(any(any((dataTest > 75)))) && not(any(any((dataTest < -75)))),
                        % Segment & concatenate verticaly the data samples of a valid trial
                        deltaEEG = vertcat(deltaEEG, EEGD.data(:,ini:fin)');
                        thetaEEG = vertcat(thetaEEG, EEGT.data(:,ini:fin)');
                        alphaEEG = vertcat(alphaEEG, EEGA.data(:,ini:fin)');
                        betaEEG = vertcat(betaEEG, EEGB.data(:,ini:fin)');
                        gammaEEG = vertcat(gammaEEG, EEGG.data(:,ini:fin)');
                    end
                end
            end
        else
            error('Two consecutive trials dont have the same code')
        end
    end
    clear EEG EEGD EEGT EEGA EEGB EEGG
end


% save to a .mat file
save(saveTo, 'deltaEEG' ,'thetaEEG' ,'alphaEEG' ,'betaEEG' ,'gammaEEG' ,'info')
    
    
    