% Estimates bumps location and creates bump index matrices with infor about
% bump location time left after bump, subject and condition

%clear all
clearvars
close all
clc

load varForBumpsOn_Res12575.mat x y x20 y20 subjectsF conds subjects
%load SegByBand125_0558.mat -regexp ^(?!info$).

BuModls = 8;                    % Number of models generated and considered
nCon = 5;                       % Number of conditions
nSub = 20;                      % number of subjects
nTrial = length(y);             % Total number of trials
nTrialS = zeros(nSub,1);        % number fo trial per subject
lenS = cell(nSub,1);            % lenght of trials per subjects
%maxId = cell(nSub,1);           % [Ntrials,Nbump] Index of samples with maximum likelihood of having a bump
meanId = [];                    % index of samples with mean probability of having abump [TotalNtrials, Nbumps, (idxBump, t_pos, subject, condition)]
bumpLoc = cell(BuModls,1);      % cell{Model}[TotalNtrials, Nbumps, (idxBump, t_pos, subject, condition)]
offset = 3;                     % mimimal distance between two consecutive bumps
nTrConS = zeros(nCon, nSub);    % number of trials per condition per subject
dataM = cell(BuModls,1);        % [BumModel, (mean, max Prob)] contains for al freq. bands and all conditions loked to bumps
sampl = 375;                    % number of samples per trial in HSMM, if the trial was shorter it si filled with zeros
nCh = 32;                       % number of channel

for c = 1:nCon
    for sj = 1:nSub
        % Counts how many trial has each condition [nCond, nSub]
        nTrConS(c,sj) = length(find(conds == c & subjects == sj));
    end
end


 % search for Bump locations
for buModl = 1:BuModls % iterates the files with different number of bumps models
    filename = ['stages125_0558_' num2str(buModl) '_Bu.mat'];
    load(filename, 'eventprob20i'); %,'eventprob20iGF') % eventprob20i []
                                        % data values are: {sj}[samples, Ntrials, Nbumps]
    conTr = 0;      % trail counter to get the conditions per trial
    %%%% ------ Find bumps locations ----
    for sj = 1:nSub,
        nTrialS = size(x20{sj},1);          % number of trials per subject
        meanId = zeros(nTrialS, buModl, 4);
        lenS = y20{sj} - x20{sj} + 1;       % trials lenght
        for tr = 1:nTrialS
            conTr = conTr + 1;              % counts trial
            % take mean prboability (Expected value) that a bump is in a sample
            meanId(tr,:,1) = round([1:lenS(tr)] * ...
                    reshape(eventprob20i{sj}(1:lenS(tr),tr,:),lenS(tr),buModl));
            for bu=1:buModl
                if  bu > 1 && meanId(tr,bu,1) <= meanId(tr,bu-1,1) + offset,
                    % Pops window if two bumps overlap
                    [answer] = pop_locMismatch(tr,lenS,eventprob20i{sj},squeeze(meanId(:,:,1)),bu)
                    meanId(tr,bu-1,1) = str2num(answer{1});
                    meanId(tr,bu,1) = str2num(answer{2});
                end
%                 if  meanId{sj}(tr,bu) <= meanId{sj}(tr,bu-1) + offset,
%                     % looks for the next sample with higest probability of
%                     % having a bump. 
%                     [~, newBu] = max(eventprob20i{sj}(meanId{sj}(tr,bu-1)+offset:lenS{sj}(tr),tr,bu));
%                     maxId{sj}(tr,bu) = newBu + meanId{sj}(tr,bu-1)+offset;
%                     dataM{buModl}.overlap = vertcat(dataM{buModl}.overlap, ...    % keeps a record of the bumps that overlap
%                         [sj, tr, bu, bu-1, meanId{sj}(tr,bu), meanId{sj}(tr,bu-1)]);  % 
%                     clear newBu
%                end
                meanId(tr,bu,2) = lenS(tr) - meanId(tr,bu,1);   % samples from bump location to end of trial
                meanId(tr,bu,3) = sj;                           % Subject code
                meanId(tr,bu,4) = conds(conTr);                 % condition code
            end
        end
        bumpLoc{buModl} = cat(1,bumpLoc{buModl}, meanId);
    end
   if conTr ~= length(x),
       error('Number of trials counted does not match')
   end
    clear eventprob20i
end
save bumpLocations.mat bumpLoc
    
