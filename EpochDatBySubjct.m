% It epoches the data respect to conditions: Target, Foil, New, Fan 1 and Fan 2.
% The trials are epoched one by one and comberted into Qiong's HSMM format for 
% bups analysis (The discovery of... Qiong, Borst, Anderson 2016)
% The new data will go to the folder - ForBupsA-200-160. In
% addition the epochs are dretended, and rejected those with RT ouliers,
% response not accurate and an amilitude out of the range [-80, 80] micro V

clear all
close all
clc

pathdata = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/';
snames = dir([pathdata '*.set']);
Nsj = size(snames,1);
pathsave = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/';
info = ['varibles for HSMM code from Qiong; 100 Hz sampling freq.; filter 1 to 35 Hz;' ...
    'lenght trials[Onset, sep], conditions: TargetFan1, TagetFan2, FoilFan1, FoilFan2;' ...
    'It epochs the data per subjects and does PCA per subject with a covariance matrix per trial'];
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
  
% Order subject files by subject number
subjctN = zeros(Nsj,1);
for sj=1:Nsj,
    number = strsplit(snames(sj).name(2:3), '_');
    subjctN(sj) = str2num(number{1});
end
[~, Isort] = sort(subjctN);
snames = snames(Isort);


for sj=1:Nsj,
%sj = 16;
    % variables that HSMM needs
    data = [];
    x = [];
    y = [];
    x5 = cell(1,5);
    y5 = cell(1,5);
    conds = [];
    condsB = [];
    covar = 0;

    % Open subject dataset
    EEG = pop_loadset('filename',snames(sj).name ,'filepath',pathdata);
    EEG = eeg_checkset( EEG );
        % Save dataset as:
    namesave = strrep(snames(sj).name,'Events235.set','epochs235.mat');
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
    %newpa = [];
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
    %newpa(1:2:end) = [];

    % do std + -
    fa1taSD(1) = mean(fa1ta) - 3 * std(fa1ta);
    fa2taSD(1) = mean(fa2ta) - 3 * std(fa2ta);
    fa1foSD(1) = mean(fa1fo) - 3 * std(fa1fo);
    fa2foSD(1) = mean(fa2fo) - 3 * std(fa2fo);
    %newpaSD(1) = mean(newpa) - 3 * std(newpa);

    fa1taSD(2) = mean(fa1ta) + 3 * std(fa1ta);
    fa2taSD(2) = mean(fa2ta) + 3 * std(fa2ta);
    fa1foSD(2) = mean(fa1fo) + 3 * std(fa1fo);
    fa2foSD(2) = mean(fa2fo) + 3 * std(fa2fo);
    %newpaSD(2) = mean(newpa) + 3 * std(newpa);        
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

            % check that both events are not RT ouliers, the response was
            % inaccurate in any of them and the reaction time is less than
            % 3 seconds
            % Use the std values computed for Short land long words together
            
            if [EEG.event(idxEvnt(tr):idxEvnt(tr+1)).accur] & ...
                    [EEG.event(idxEvnt(tr)).RT < 3000],
                
                if (trialCode == 1 & [EEG.event(idxEvnt(tr)).RT > fa1taSD(1)] & [EEG.event(idxEvnt(tr)).RT < fa1taSD(2)]) | ...
                        (trialCode == 2 & [EEG.event(idxEvnt(tr)).RT > fa2taSD(1)] & [EEG.event(idxEvnt(tr)).RT < fa2taSD(2)]) | ...
                        (trialCode == 3 & [EEG.event(idxEvnt(tr)).RT > fa1foSD(1)] & [EEG.event(idxEvnt(tr)).RT < fa1foSD(2)]) | ...
                        (trialCode == 4 & [EEG.event(idxEvnt(tr)).RT > fa2foSD(1)] & [EEG.event(idxEvnt(tr)).RT < fa2foSD(2)]),% | ...
                        %(trialCode == 5 & [EEG.event(idxEvnt(tr)).RT > newpaSD(1)] & [EEG.event(idxEvnt(tr)).RT < newpaSD(2)]),
                        % check std per condition
                
                    % Extract data epoch, detrend it and test amplitude artefacts
                    dataTest = detrend(EEG.data(:,ini:fin)', 'linear');
                    if not(any(any((dataTest > 80)))) && not(any(any((dataTest < -80)))),
                        % Do all variable for the HSMM
                        x = vertcat(x, size(data,1) + 1);   % first samle of this trial in x; initialy data = []
                        x5{trialCode} = vertcat(x5{trialCode}, size(data,1) + 1); % begining of trials by condition

                        data = vertcat(data, dataTest);   % concatenate verticaly the data samples of a valid trial

                        y = vertcat(y, size(data,1));  % last sample of this trial
                        y5{trialCode} = vertcat(y5{trialCode}, size(data,1)); % last sample of a trial by conditions
                        conds = vertcat(conds, trialCode);  %  condition of  each trial 
                        condsB =  vertcat(condsB, repmat(conds(end),size(dataTest,1),1)); % Condition of each sample
                        subject = sj;
                        % Trial covariance matrix
                        dataTest = dataTest'; % transposes data test to do PCA like in Cohens's book
                        covar = covar + (dataTest*dataTest')./(size(dataTest,2)-1);
                    end
                end
            end
        else
            error('Two consecutive trials dont have the same code')
        end
    end
    % PCA
    covar = covar./length(x);   % mean covariance of all trials in a subject
    [coeff10,latent10] = pcacov(covar);
    %coeff10 = coeff(:,1:10);
    %latent10 = latent(1:10);
    if size(coeff10,1) ~= size(data,2), error('Error: PC matrix does nto macht data'), end
    % PCA time seriers
    score10 = data * coeff10;
    normedscore10=zeros(size(score10));
    means10 = mean(data,1);
     % normalize (zscore) by trial
    for i = 1:length(x)
        normedscore10(x(i):y(i),:)=zscore(score10(x(i):y(i),:));
    end
    normedscore10 = normedscore10(:,1:10);
    % save to a .mat file
    save(saveTo, 'info','data','x','x5','y','y5','conds','condsB', ...
        'subject','coeff10','score10','latent10','normedscore10', 'means10');
    clear EEG normedscore10 dataTest coeff coeff10 latent latent10 score10
end

    
    
    