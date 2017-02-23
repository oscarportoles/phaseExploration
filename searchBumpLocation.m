% Estimates bumps location and creates bump index matrices.
% It genrates an enormous amount of data so I rejected this approach

%clear all
clearvars
close all
clc

load varForBumpsOn_Res12575.mat x y x20 y20 subjectsF conds subjects
load SegByBand125_0558.mat -regexp ^(?!info$).

BuModls = 8;                    % Number of models generated and considered
nCon = 5;                       % Number of conditions
nSub = 20;                      % number of subjects
nTrial = length(y);             % Total number of trials
nTrialS = zeros(nSub,1);        % number fo trial per subject
lenS = cell(nSub,1);            % lenght of trials per subjects
%maxId = cell(nSub,1);           % [Ntrials,Nbump] Index of samples with maximum likelihood of having a bump
meanId = cell(nSub,1);          % [Ntrials,Nbumps] index of samples with mean probability of having abump
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

% Creata data structures with zeros dataM{i}.bump.band.condition
for c = 1:nCon
    fieldC = ['co_' num2str(c)];
    value = zeros(0, 2*sampl, nCh, nSub);
    alpha.(fieldC) =  value;
    beta.(fieldC) = value;
    delta.(fieldC) = value;
    theta.(fieldC) = value;
    gamma.(fieldC) = value;
end
temst.('alpha') = alpha;
temst.('beta') = beta;
temst.('delta') = delta;
temst.('theta') = theta;
temst.('gamma') = gamma;
for buModl = 1:BuModls
    for bu = 1:buModl
        fieldB = ['bu_' num2str(bu)];
        dataM{buModl,1}.(fieldB) =  temst;      % locked to mean location
        %dataM{buModl,2}.(fieldB) =  temst;      % lock to max likely location
    end 
end
clear fieldB fieldC temst value alpha beta delta theta gamma bu buModl c nCh nCon


 % search for Bump locations
for buModl = 1:BuModls % iterates the files with different number of bumps models
    filename = ['stages125_0558_' num2str(buModl) '_Bu.mat'];
    load(filename, 'eventprob20i'); %,'eventprob20iGF') % eventprob20i []
    dataM{buModl}.overlap = [];         % keeps track of the bumps that overlap and are corrected 
                                        % data values are: {sj}[samples, Ntrials, Nbumps]
    
    %%%% ------ Find bumps locations ----
    for sj = 1:nSub,
        nTrialS(sj) = size(x20{sj},1);          % number of trials per subject
        lenS{sj} = y20{sj} - x20{sj} + 1;       % trials lenght
% We will only use the mean expected locations        
%         for bu = 1:buModl,
%             % take maximaly likely sample
%             [~, maxId{sj}(:,bu)] = max(eventprob20i(:,:,bu,sj));
%             for tr = 1:nTrialS(sj),
%                 if bu > 1 && maxId{sj}(tr,bu) <= maxId{sj}(tr,bu-1) + offset,
%                     % Pops window if two bumps overlap 
%                     [answer] = pop_locMismatch(sj,tr,lenS,eventprob20i,maxId)
%                     maxId{sj}(tr,bu-1) = str2num(answer{1});
%                     maxId{sj}(tr,bu) = str2num(answer{2});
%                 end
%             end
%         end
        for tr = 1:nTrialS(sj)
            % take mean prboability sample
            meanId{sj}(tr,:) = round([1:lenS{sj}(tr)] * ...
                    reshape(eventprob20i{sj}(1:lenS{sj}(tr),tr,:),lenS{sj}(tr),buModl));
            for bu=2:buModl
%                 if  meanId{sj}(tr,bu) <= meanId{sj}(tr,bu-1) + offset,
%                     % Pops window if two bumps overlap
%                     [answer] = pop_locMismatch(tr,lenS{sj},eventprob20i{sj},meanId{sj})
%                     meanId{sj}(tr,bu-1) = str2num(answer{1});
%                     meanId{sj}(tr,bu) = str2num(answer{2});
%                 end
                if  meanId{sj}(tr,bu) <= meanId{sj}(tr,bu-1) + offset,
                    % looks for the next sample with higest probability of
                    % having a bump. 
                    [~, newBu] = max(eventprob20i{sj}(meanId{sj}(tr,bu-1)+offset:lenS{sj}(tr),tr,bu));
                    maxId{sj}(tr,bu) = newBu + meanId{sj}(tr,bu-1)+offset;
                    dataM{buModl}.overlap = vertcat(dataM{buModl}.overlap, ...    % keeps a record of the bumps that overlap
                        [sj, tr, bu, bu-1, meanId{sj}(tr,bu), meanId{sj}(tr,bu-1)]);  % 
                    clear newBu
                end
            end
        end 
    end
    %%%%%%%%%%%% Hilbert transform and bump-lock the data %%%%%%%%%%%%%%%%%
    meanSj = [];
    maxSj = [];
    for sj = 1:nSub
        % searches the sample where the trials will be locked to
        meanSj = vertcat(meanSj,mean(meanId{sj}));
        maxSj = vertcat(maxSj, max(meanId{sj}));
    end
    lockIx = floor(mean(meanSj) + max(maxSj));  % loIckxMean[buModl]    
    dataM{buModl}.lockIdx = lockIx;             % put locking indexes in the data struct
    dataM{buModl}.bumpLoc = meanId;             % {sj}[tr,bu]Save bumps location per subject, trial and bump             
    for sj = 1:nSub   
        for bu = 1:buModl
            for tr = 1:length(y20{sj}), % Trials per subject
                % hilbert transform each trial. It trasnforms columns wise
                dataA = hilbert(alphaEEG(x20{sj}(tr):y20{sj}(tr),:));
                dataB = hilbert(betaEEG(x20{sj}(tr):y20{sj}(tr),:));
                dataD = hilbert(deltaEEG(x20{sj}(tr):y20{sj}(tr),:));
                dataT = hilbert(thetaEEG(x20{sj}(tr):y20{sj}(tr),:));
                dataG = hilbert(gammaEEG(x20{sj}(tr):y20{sj}(tr),:));
                %%% Put data bump locked to the new data variables
                % Alpha band
                pre = dataA(1:meanId{sj}(tr,bu), :);
                pos = dataA(meanId{sj}(tr,bu)+1:end,:);
                %dataM{buModl,1}.bu_1.alpha.c_1(end+1,lockIx(bu)-size(preA,1)+1:lockIx(bu)) = preA;
                %dataM{buModl,1}.bu_1.alpha.c_1(end, lockIx(bu)+1:lockIx(bu)+size(posA,1)) = posA;
                eval(['dataM{buModl,1}.bu_' num2str(bu) '.alpha.co_' num2str(conds(tr)) ...
                    '(end+1,lockIx(bu)-size(pre,1)+1:lockIx(bu),:,sj) = pre;']);
                eval(['dataM{buModl,1}.bu_' num2str(bu) '.alpha.co_' num2str(conds(tr)) ...
                    '(end, lockIx(bu)+1:lockIx(bu)+size(pos,1),:,sj) = pos;']);
                % Beta band
                pre = dataB(1:meanId{sj}(tr,bu), :);
                pos = dataB(meanId{sj}(tr,bu)+1:end,:);
                eval(['dataM{buModl,1}.bu_' num2str(bu) '.beta.co_' num2str(conds(tr)) ...
                    '(end+1,lockIx(bu)-size(pre,1)+1:lockIx(bu),:,sj) = pre;']);
                eval(['dataM{buModl,1}.bu_' num2str(bu) '.beta.co_' num2str(conds(tr)) ...
                    '(end, lockIx(bu)+1:lockIx(bu)+size(pos,1),:,sj) = pos;']);
                % Delta band
                pre = dataD(1:meanId{sj}(tr,bu), :);
                pos = dataD(meanId{sj}(tr,bu)+1:end,:);
                eval(['dataM{buModl,1}.bu_' num2str(bu) '.delta.co_' num2str(conds(tr)) ...
                    '(end+1,lockIx(bu)-size(pre,1)+1:lockIx(bu),:,sj) = pre;']);
                eval(['dataM{buModl,1}.bu_' num2str(bu) '.delta.co_' num2str(conds(tr)) ...
                    '(end, lockIx(bu)+1:lockIx(bu)+size(pos,1),:,sj) = pos;']);            
                % Theta band
                pre = dataT(1:meanId{sj}(tr,bu), :);
                pos = dataT(meanId{sj}(tr,bu)+1:end,:);
                eval(['dataM{buModl,1}.bu_' num2str(bu) '.theta.co_' num2str(conds(tr)) ...
                    '(end+1,lockIx(bu)-size(pre,1)+1:lockIx(bu),:,sj) = pre;']);
                eval(['dataM{buModl,1}.bu_' num2str(bu) '.theta.co_' num2str(conds(tr)) ...
                    '(end, lockIx(bu)+1:lockIx(bu)+size(pos,1),:,sj) = pos;']);
                % Gamma band
                pre = dataG(1:meanId{sj}(tr,bu), :);
                pos = dataG(meanId{sj}(tr,bu)+1:end,:);
                eval(['dataM{buModl,1}.bu_' num2str(bu) '.gamma.co_' num2str(conds(tr)) ...
                    '(end+1,lockIx(bu)-size(pre,1)+1:lockIx(bu),:,sj) = pre;']);
                eval(['dataM{buModl,1}.bu_' num2str(bu) '.gamma.co_' num2str(conds(tr)) ...
                    '(end, lockIx(bu)+1:lockIx(bu)+size(pos,1),:,sj) = pos;']);
            end
        end
    end
end
save bumpLockData.mat dataM meanId lockIx
