% Estimates bumps location and creates bump index matrices with infor about
% bump location time left after bump, subject and condition

clear all
%clearvars
close all
clc

load varForBumpsOn_Res12575.mat x y x20 y20 subjectsF conds subjects
%load SegByBand125_0558.mat -regexp ^(?!info$).

BuModls = 8;                    % Number of models generated and considered
nCon = 5;                       % Number of conditions
nSub = 20;                      % number of subjects
nTrial = length(y);             % Total number of trials
%nTrialS = zeros(nSub,1);        % number fo trial per subject
lenS = cell(nSub,1);            % lenght of trials per subjects
maxId = cell(nSub,1);           % [Ntrials,Nbump] Index of samples with maximum likelihood of having a bump
bumpLoc = cell(BuModls,1);      % cell{Model}[TotalNtrials, Nbumps, (idxBump, t_pos, subject, condition)]
nTrConS = zeros(nCon, nSub);    % number of trials per condition per subject
sampl = 375;                    % number of samples per trial in HSMM, if the trial was shorter it si filled with zeros
nCh = 32;                       % number of channel
noNaNtr = zeros(BuModls, nSub); % number of non NaN elements per subject

 % search for Bump locations
for buModl = 1:BuModls % iterates the files with different number of bumps models
    filename = ['stages125_0558_' num2str(buModl) '_Bu.mat'];
    load(filename, 'eventprob20i'); %,'eventprob20iGF') % eventprob20i []
                                        % data values are: {sj}[samples, Ntrials, Nbumps]
    conTr = 0;                      % trail counter to get the conditions per trial
    %%%% ------ Find bumps locations ----
    for sj = 1:nSub,
        nTrialS = size(x20{sj},1);          % number of trials per subject
        noNaNtr(1,sj) = nTrialS;            % number of trials per subject
        maxId = NaN(nTrialS, buModl, 4);
        lenS = y20{sj} - x20{sj} + 1;       % trials lenght
        for tr = 1:nTrialS,
            conTr = conTr + 1;              % counts trial
            buOut = 0;
            for bu = 1:buModl,
                % Don't take PDF with ouliers or many peaks
                if not(buOut) && kurtosis(eventprob20i{sj}(:,tr,bu)) > 35 % before 30
                    % take maximaly likely sample
                    [~, maxId(tr,bu,1)] = max(eventprob20i{sj}(:,tr,bu));
                else
                    buOut = 1;
                    buNaN = 1:buModl;
                    maxId(tr,buNaN,1) = NaN;
                end
                maxId(tr,bu,2) = lenS(tr) - maxId(tr,bu,1);   % samples from bump location to end of trial
                maxId(tr,bu,3) = sj;                           % Subject code
                maxId(tr,bu,4) = conds(tr);                 % condition co
            end
        end
        bumpLoc{buModl} = cat(1,bumpLoc{buModl}, maxId);
        noNaNtr(buModl+1,sj) = sum(~isnan(maxId(:,1,1)), 1);
   end
   if conTr ~= length(x),
       error('Number of trials counted does not match')
   end
    clear eventprob20i
end
save bumpLocationsMax35.mat bumpLoc noNaNtr
    
