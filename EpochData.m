% It epoches the data according to the events previously introduced.
% The new data will go to the folders - lockOnset and lockResponse. In
% addition the epochs are dretended, and rejected those with RT ouliers,
% response not accurate and an amilitude out of the range [-75, 75] micro V

clear all
close all
clc

% pathevents = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/Corrected/';
% fnames = dir([pathevents '*.txt']);
% Nsj = size(fnames,1);

pathdata = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/Downsamp125Lowpass58/events/';
snames = dir([pathdata '*.set']);
Nsj = size(snames,1);

pathsaveOnset = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/Downsamp125Lowpass58/events/lockOnset/';
pathsavResponse = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/Downsamp125Lowpass58/events/lockResponse/';
paths = {pathsaveOnset ;pathsavResponse };

epochsOn = [{'OnFoilFan1L'},{'OnFoilFan1S'},{'OnFoilFan2L'},{'OnFoilFan2S'}, ...
    {'OnTargetFan1L'},{'OnTargetFan1S'},{'OnTargetFan2L'},{'OnTargetFan2S'}, ...
    {'OnNewFan1S'}, {'OnNewFan1L'}];
epochsRes = [{'ResFoilFan1L'},{'ResFoilFan1S'},{'ResFoilFan2L'},{'ResFoilFan2S'}, ...
    {'ResTargetFan1L'},{'ResTargetFan1S'},{'ResTargetFan2L'},{'ResTargetFan2S'}, ...
    {'ResNewFan1S'}, {'ResNewFan1L'}];
epochnames = [epochsOn ;epochsRes ];

% % longOnset = [-0.7 2];  % duration [seconds] of a Onset lock trial from the locking event 
% % longResp = [-2 0.7];   % duration [seconds] of a Response lock trial from the locking event
% % longTrial = [longOnset;longResp];

baseLineOn = [-300 0]; % base line correction reference [millisecondes] Onset lock trials
baseLineRes = [0 300]; % base line correction reference [millisecondes] Response lock trials
baseLine = [baseLineOn;baseLineRes];

for sj=1:Nsj,
%sj = 16;
    % Open subject dataset
    EEG = pop_loadset('filename',snames(sj).name ,'filepath',pathdata);
    EEG = eeg_checkset( EEG );
    EEGo = EEG;   % Copy of the original EEG
    
    % Below it estimates epoch length according to the mean RT and the SD of the RT.
    % Each subject will have its own epoch length. Take into account for
    % Grand Average ERP. I do it because there is big variability in
    % between subjects RT
    meanRT = mean(single([EEG.event.RT]));
    stdRT = std(single([EEG.event.RT]));
    epochEnd = ((meanRT + 3 * stdRT) / 1000) + 0.5; % rounds-oof to upper sample 
    sjLengthOn = [-0.7 epochEnd];
    
    for lk = 1:1, % runs throught onset lock and Response lock epochs #### Now only locks to the Onset #####
        for ep = 1:length(epochsOn), % Iterate all epoch locked to the Onset of the stimuli
            nameEp = ['DS125_', epochnames{lk,ep}];
            % With first option commentd runs onlu through One locking
            %EEG = pop_epoch(EEGo, {epochnames{lk,ep}}, longTrial(lk,:), 'newname', nameEp, 'epochinfo', 'yes');
            EEG = pop_epoch(EEGo, {epochnames{lk,ep}}, sjLengthOn, 'newname', nameEp, 'epochinfo', 'yes');
            EEG = pop_rmbase(EEG, baseLine(lk,:));
            toRemove = [];
            for tr=1:EEG.trials, % iterate all trials in a newly created epoch
                if length(EEG.epoch(tr).event) == 1 || ...              % only two event per trial
                        not(any([EEG.epoch(tr).eventaccur{1,2}])) || ...   % The respose was right
                        any([EEG.epoch(tr).eventoutRT{1,2}]),              % RT whithin mean + 3*std
                    toRemove(end+1) = tr;
                else    % Pass first test
                    data = EEG.data(:,:,tr);
                    if any(any((data > 80))) || any(any((data < -80))),  % amplitude out of range (-75, 75)
                        toRemove(end+1) = tr;
                    else  % good trials
                        data = detrend(data', 'linear');    % remove linear trend over columns in a matrix
                        EEG.data(:,:,tr) = data';
                    end
                end
            end
            % Remove all elements that did't pass the test above
            EEG.epoch(toRemove) = [];
            EEG.data(:,:,toRemove) = [];
            %EEG.icaact(:,:,toRemove) = [];
            EEG.trials = EEG.trials - length(toRemove);
            EEG = eeg_checkset( EEG );
            EEG = pop_rmbase(EEG, baseLine(lk,:));
            
            name = strsplit(snames(sj).name,'_');
            newname = [name{1,1} '_DS125_' epochnames{lk,ep} '.set'];
            EEG.filename = newname;
            EEG = pop_saveset( EEG, 'filename',newname,'filepath',paths{lk});
            EEG = eeg_checkset( EEG );
        end
    end

    clear EEG
    
end         
    
    
    