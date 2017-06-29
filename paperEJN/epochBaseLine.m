% It epochs the baseline of band pass filtered data. The epochs are choosen based on the 
% labels of each event on the .TXT document (preious step), the reaction time, the amplitude 
% threshods on raw data. Segments are the same for all frequency bands. New data is 
% separeted by conditions and freuency bands. Epochs of baseline are the same as epochs of 
% post stimulus data

clear all
close all
clc

% Full spectrum data set
pathdataO = '/Data_JB_AssoRec/DS100bp2_35/events/';
snamesO = dir([pathdataO '*.set']);
NsjO = size(snamesO,1);
% % Delta band data
% pathdataD = '/Data_JB_AssoRec/DS100_EEGbands/Delta/';
% snamesD = dir([pathdataD '*.set']);
% Theta band data
pathdataT = '/Data_JB_AssoRec/DS100_EEGbands/Theta/';
snamesT = dir([pathdataT '*.set']);
% Alpha band data
pathdataA = '/Data_JB_AssoRec/DS100_EEGbands/Alpha/';
snamesA = dir([pathdataA '*.set']);
% Beta band data
pathdataB = '/Data_JB_AssoRec/DS100_EEGbands/Beta/';
snamesB = dir([pathdataB '*.set']);
% % Gamma band data
% pathdataG = '/Data_JB_AssoRec/DS100_EEGbands/Gamma/';
% snamesG = dir([pathdataG '*.set']);

% SAve data to path
pathsave = '/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/';

% test that all subjects in all frequency bands are loadad correctly
% Order subject files by subject number
subjctNO = zeros(NsjO,1);
for sj=1:NsjO,
    number = strsplit(snamesO(sj).name(2:3), '_');
    subjctNO(sj) = str2num(number{1});
%     number = strsplit(snamesD(sj).name(2:3), '_');
%     subjctND(sj) = str2num(number{1});
    number = strsplit(snamesT(sj).name(2:3), '_');
    subjctNT(sj) = str2num(number{1});
    number = strsplit(snamesA(sj).name(2:3), '_');
    subjctNA(sj) = str2num(number{1});
    number = strsplit(snamesB(sj).name(2:3), '_');
    subjctNB(sj) = str2num(number{1});
%     number = strsplit(snamesG(sj).name(2:3), '_');
%     subjctNG(sj) = str2num(number{1});
end
[~, Isort] = sort(subjctNO);
snamesO = snamesO(Isort);
% [~, Isort] = sort(subjctND);
% snamesD = snamesD(Isort);
[~, Isort] = sort(subjctNT);
snamesT = snamesT(Isort);
[~, Isort] = sort(subjctNA);
snamesA = snamesA(Isort);
[~, Isort] = sort(subjctNB);
snamesB = snamesB(Isort);
% [~, Isort] = sort(subjctNG);
% snamesG = snamesG(Isort);

% Check that all frequency bands have the same number of subjects and all
% subjects are the same
% if (length(snamesO) == length(snamesD) == length(snamesT) == ...
%     length(snamesA) == length(snamesB) == length(snamesG)) | ...
%     (sum(subjctNO) == sum(subjctND) == sum(subjctNT) == ...
%     sum(subjctNA) == sum(subjctNB) == sum(subjctNG))
%     error('Not all frequency bands have the same subjects')
% end
if (length(snamesO) == length(snamesT) == ...
    length(snamesA) == length(snamesB)) | ...
    (sum(subjctNO) == sum(subjctNT) == ...
    sum(subjctNA) == sum(subjctNB)),
    error('Not all frequency bands have the same subjects')
end

info = ['segment bandpass filtered data like if it were for a bumps HSMM;' ...
    '100 Hz sampling freq.; filter 2 to 35 Hz;' ...
    'lenght trials[Onset, REsp], conditions: TargetFan1, TagetFan2, FoilFan1, FoilFan2;' ...
    'variables [npts,Ch]; Delta band 2 - 4 Hz, Theta band 4 - 9 Hz, ' ...
    'Alpha band 9 - 14 Hz, Beta band 14 - 30 Hz, Gamma band 30 - 58 Hz'];

% % Define epoching boundaries on time
% init = -0.2;    % time [s] beginin of trial from Onset
% fint = 0.16;    % time [s] end of trial from Response

% Define epoching boundaries on time
init = 0;    % time [s] beginin of trial from Onset
fint = 0;    % time [s] end of trial from Response
initBL = 0.3;   % time [s] begining of baseline period, the end is init - 1 sample

trCode =   [[2 3 1]  % first column: Jelmer's code for short words
            [4 5 2]  % second column: Jelmer's code for long words
            [6 7 3]  % third column: new code
            [8 9 4]
            [10 11 5]];

        
for sj=1:NsjO,
    % variable to keep the data
%     deltaEEG = [];
    thetaEEGbl = [];
    alphaEEGbl = [];
    betaEEGbl = [];
%     gammaEEG = [];
    x = [];
    y = [];
    x5 = cell(1,5);
    y5 = cell(1,5);
    conds = [];
    condsB = [];
    
        % Open subject dataset
    EEG = pop_loadset('filename',snamesO(sj).name ,'filepath',pathdataO);
    EEG = eeg_checkset( EEG );
    
%     EEGD = pop_loadset('filename',snamesD(sj).name ,'filepath',pathdataD);
%     EEGD = eeg_checkset( EEGD );
    EEGT = pop_loadset('filename',snamesT(sj).name ,'filepath',pathdataT);
    EEGT = eeg_checkset( EEGT );
    EEGA = pop_loadset('filename',snamesA(sj).name ,'filepath',pathdataA);
    EEGA = eeg_checkset( EEGA );
    EEGB = pop_loadset('filename',snamesB(sj).name ,'filepath',pathdataB);
    EEGB = eeg_checkset( EEGB );
%     EEGG = pop_loadset('filename',snamesG(sj).name ,'filepath',pathdataG);
%     EEGG = eeg_checkset( EEGG );

        % Save dataset as:
    idSj = regexp(snamesO(sj).name, '_DS1');
    namesave = [snamesO(sj).name(1:idSj) 'PreStimEpochBands100.mat'];
    saveTo = [pathsave namesave];

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
%     newpa = [];
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
%             elseif EEG.event(i).code == 10 | EEG.event(i).code == 11,
%                 newpa(end+1) = EEG.event(i).RT;
            end
        end
    end  
    
    % remove duplicated trials
    fa1ta(1:2:end) = [];
    fa2ta(1:2:end) = [];
    fa1fo(1:2:end) = [];
    fa2fo(1:2:end) = [];
%     newpa(1:2:end) = [];

    % do std + -
    fa1taSD(1) = mean(fa1ta) - 3 * std(fa1ta);
    fa2taSD(1) = mean(fa2ta) - 3 * std(fa2ta);
    fa1foSD(1) = mean(fa1fo) - 3 * std(fa1fo);
    fa2foSD(1) = mean(fa2fo) - 3 * std(fa2fo);
%     newpaSD(1) = mean(newpa) - 3 * std(newpa);

    fa1taSD(2) = mean(fa1ta) + 3 * std(fa1ta);
    fa2taSD(2) = mean(fa2ta) + 3 * std(fa2ta);
    fa1foSD(2) = mean(fa1fo) + 3 * std(fa1fo);
    fa2foSD(2) = mean(fa2fo) + 3 * std(fa2fo);
%     newpaSD(2) = mean(newpa) + 3 * std(newpa);        
    clear fa1ta fa2ta fa1fo fa2fo %newpa
    
    % check if the number of event is odd -> if not a event is missing
    if mod(length(idxEvnt),2), error('The number of events is not even'); end
 
    for tr = 1:2:length(idxEvnt)-1, % Iterate all events to segment
        % Check that consecutive events belong to the same trial
        if EEG.event(idxEvnt(tr)).code == EEG.event(idxEvnt(tr+1)).code, 
            % get the condition code of this trial
            [row,~] = find((trCode(:,1:2) == EEG.event(idxEvnt(tr)).code));
            trialCode = trCode(row,3);
            
            % Define epoch boundaries on samples
            ini = floor(EEG.event(idxEvnt(tr)).latency + init * EEG.srate + 1); % attention! add +1 to get the same segments ans in Jelmer's data 
            fin = round(EEG.event(idxEvnt(tr+1)).latency + fint * EEG.srate);
            iniBL = ini - round(initBL * EEG.srate);
            finBL = ini - 1;
            durBL = finBL - iniBL + 1;
            
            % Use the std values computed for Short land long words together  
            if [EEG.event(idxEvnt(tr):idxEvnt(tr+1)).accur] & ...
                    [EEG.event(idxEvnt(tr)).RT < 3000],
                
                if (trialCode == 1 & [EEG.event(idxEvnt(tr)).RT > fa1taSD(1)] & [EEG.event(idxEvnt(tr)).RT < fa1taSD(2)]) | ...
                        (trialCode == 2 & [EEG.event(idxEvnt(tr)).RT > fa2taSD(1)] & [EEG.event(idxEvnt(tr)).RT < fa2taSD(2)]) | ...
                        (trialCode == 3 & [EEG.event(idxEvnt(tr)).RT > fa1foSD(1)] & [EEG.event(idxEvnt(tr)).RT < fa1foSD(2)]) | ...
                        (trialCode == 4 & [EEG.event(idxEvnt(tr)).RT > fa2foSD(1)] & [EEG.event(idxEvnt(tr)).RT < fa2foSD(2)]),% | ...
%                         (trialCode == 5 & [EEG.event(idxEvnt(tr)).RT > newpaSD(1)] & [EEG.event(idxEvnt(tr)).RT < newpaSD(2)]),
                        % check std per condition
                    % Extract data epoch, detrend it and test amplitude artefacts
                    dataTest = detrend(EEG.data(:,ini:fin)', 'linear');
                    if not(any(any((dataTest > 80)))) && not(any(any((dataTest < -80)))),
                        x = vertcat(x, size(thetaEEGbl,1) + 1);   % first sample of this trial in x; initialy thetaEEGbl = []
                        x5{trialCode} = vertcat(x5{trialCode}, size(thetaEEGbl,1) + 1); % begining of trials by condition

                        % Segment & concatenate verticaly the data samples of a valid trial
                        beforSz = size(thetaEEGbl,1);
                        %deltaEEG = vertcat(deltaEEG, EEGD.data(:,ini:fin)');
                        thetaEEGbl = vertcat(thetaEEGbl, EEGT.data(:,iniBL:finBL)');
                        alphaEEGbl = vertcat(alphaEEGbl, EEGA.data(:,iniBL:finBL)');
                        betaEEGbl = vertcat(betaEEGbl, EEGB.data(:,iniBL:finBL)');
                        %gammaEEG = vertcat(gammaEEG, EEGG.data(:,ini:fin)');
                        % test that the baseline size is consistent
                        if size(thetaEEGbl,1) - beforSz ~= durBL,
                            error('The size of the baseline is not consisten')
                        end
                        
                        y = vertcat(y, size(thetaEEGbl,1));  % last sample of this trial
                        y5{trialCode} = vertcat(y5{trialCode}, size(thetaEEGbl,1)); % last sample of a trial by conditions
                        conds = vertcat(conds, trialCode);  %  condition of  each trial 
                        condsB =  vertcat(condsB, repmat(conds(end),size(dataTest,1),1)); % Condition of each sample
                        subject = sj;
                    end
                end
            end
        else
            error('Two consecutive trials dont have the same code')
        end
    end
    % save to a .mat file
    save(saveTo,'thetaEEGbl' ,'alphaEEGbl' ,'betaEEGbl','info', 'x', 'y', 'x5', 'y5', 'conds', 'condsB', 'subject')
    clear EEG EEGT EEGA EEGB
    %clear EEG EEGD EEGT EEGA EEGB EEGG
end



    
    
    
