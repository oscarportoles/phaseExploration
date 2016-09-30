% It epoches the data respect to conditions: Target, Foil, New, Fan 1 and Fan 2.
% The trials are epoched one by one and comberted into Qiong's HSMM format for 
% bups analysis (The discovery of... Qiong, Borst, Anderson 2016)
% The new data will go to the folder - ForBupsA-200-160. In
% addition the epochs are dretended, and rejected those with RT ouliers,
% response not accurate and an amilitude out of the range [-80, 80] micro V

clear all
close all
clc

pathdata = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125bp05_58/events/';
snames = dir([pathdata '*.set']);
Nsj = size(snames,1);

% Order subject files by subject number
subjctN = zeros(Nsj,1);
for sj=1:Nsj,
    number = strsplit(snames(sj).name(2:3), '_');
    subjctN(sj) = str2num(number{1});
end
[~, Isort] = sort(subjctN);
snames = snames(Isort);

pathsave = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125bp05_58/events/ForBumps_OnRes_no20/';
namesave = 'varForBumpsOn_Res12575.mat';
saveTo = [pathsave namesave];

info = ['varibles for HSMM code from Qiong; 125 Hz sampling freq.; filter 0.5 to 58 Hz;' ...
    'lenght trials[Onset, sep], conditions: TargetFan1, TagetFan2, FoilFan1, FoilFan2' ...
    ' NewFan1, NewFan2'];

epochsOn = [[{'OnFoilFan1L'},{'OnFoilFan1S'}];[{'OnFoilFan2L'},{'OnFoilFan2S'}]; ... 
    [{'OnTargetFan1L'},{'OnTargetFan1S'}];[{'OnTargetFan2L'},{'OnTargetFan2S'}]; ...
    [{'OnNewFan1S'}, {'OnNewFan1L'}]];


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

% variables that HSMM needs
data = [];
x = [];
y = [];
x5 = cell(1,5);
y5 = cell(1,5);
x20 = cell(1,Nsj);
y20 = cell(1,Nsj);
conds = [];
condsB = [];
subjects = [];
subjectsF = [];

for sj=1:Nsj,
%sj = 16;
    % Open subject dataset
    EEG = pop_loadset('filename',snames(sj).name ,'filepath',pathdata);
    EEG = eeg_checkset( EEG );

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
            

            % check that both events are not RT ouliers, the response was
            % inaccurate in any of them and the reaction time is less than
            % 3 seconds
% Selection based on .outRT, it was computed separating short and long
% words
%             if [EEG.event(idxEvnt(tr):idxEvnt(tr+1)).accur] & ...    
%                     not([EEG.event(idxEvnt(tr):idxEvnt(tr+1)).outRT]) &...
%                     [EEG.event(idxEvnt(tr)).RT < 3000],

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
                        % Do all variable for the HSMM
                        x = vertcat(x, size(data,1) + 1);   % first samle of this trial in x
                        x20{sj} = vertcat(x20{sj}, size(data,1)+1);     %  begining of trials by subject
                        x5{trialCode} = vertcat(x5{trialCode}, size(data,1) + 1); % begining of trials by condition

                        data = vertcat(data, dataTest);   % concatenate verticaly the data samples of a valid trial

                        y = vertcat(y, size(data,1));  % last sample of this trial
                        y20{sj} = vertcat(y20{sj}, size(data,1));   % last sample of a trial by subjects
                        y5{trialCode} = vertcat(y5{trialCode}, size(data,1)); % last sample of a trial by conditions
                        conds = vertcat(conds, trialCode);  %  condition of  each trial 
                        subjects = vertcat(subjects, sj);   % subject of each trial

                        condsB =  vertcat(condsB, repmat(conds(end),size(dataTest,1),1)); % Condition of each sample
                        subjectsF = vertcat(subjectsF, repmat(sj,size(dataTest,1),1));  % subject of each sample 
                    end
                end
            end
        else
            error('Two consecutive trials dont have the same code')
        end
    end
    clear EEG
end

% Do PCA and z-scoring

[coeff10 score10 latent10]= pca(data);

normedscore10=zeros(size(score10));
for i = 1:length(x)
    normedscore10(x(i):y(i),:)=zscore(score10(x(i):y(i),:)); % normalize by trial
end
normedscore10=normedscore10(:,1:10); % take first 10 components

% save to a .mat file
save(saveTo, 'info','data','x','x5','x20','y','y5','y20','conds','condsB', ...
    'subjects','subjectsF','coeff10','score10','latent10','normedscore10');
    
    
    