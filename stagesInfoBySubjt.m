% Estimate the temporal duration of each stage per condition (Target and Foils)
% and for all. To include 'New' condition the code needs changing.
% Estimation is done subject by subject.

clear all
close all
clc

nBump = 5;
fs = 100; % sampling frequency [Hz]
Nsj = 20;

% pathdata = '/home/oscar/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/';
% filecode = ['*_DS100HSMM' num2str(nBump) 'Bout235.mat'];
% snames = dir([pathdata filecode]);
% Nsj = length(snames);
% varnames = dir([pathdata '*ds100.mat']);
if nBump == 5,
    dur = struct('name',{},'nTrial',{},'nTrialTarg',{},'nTrialFoil',{},'avDur',{},'avDurTarg',{},'avDurFoil',{}, ...
        'stag1',{},'stag2',{},'stag3',{},'stag4',{},'stag5',{},'stag6',{}, ...
        'stag1Targ',{},'stag2Targ',{},'stag3Targ',{},'stag4Targ',{},'stag5Targ',{},'stag6Targ',{}, ...
        'stag1Foil',{},'stag2Foil',{},'stag3Foil',{},'stag4Foil',{},'stag5Foil',{},'stag6Foil',{}, ...
        'modLikeli',{}, 'fs', {});
elseif nBump == 6
    dur = struct('name',{},'nTrial',{},'nTrialTarg',{},'nTrialFoil',{},'avDur',{},'avDurTarg',{},'avDurFoil',{}, ...
        'stag1',{},'stag2',{},'stag3',{},'stag4',{},'stag5',{},'stag6',{},'stag7',{}, ...
        'stag1Targ',{},'stag2Targ',{},'stag3Targ',{},'stag4Targ',{},'stag5Targ',{},'stag6Targ',{},'stag7Targ',{}, ...
        'stag1Foil',{},'stag2Foil',{},'stag3Foil',{},'stag4Foil',{},'stag5Foil',{},'stag6Foil',{},'stag7Foil',{}, ...
        'modLikeli',{}, 'fs', {});
end

% for sj = 1:Nsj,
%     % Load data sets
%     fileload = [pathdata snames(sj).name];
%     load(fileload);
%     clear bumpMag gammPara
%     fileload = [pathdata varnames(sj).name];
%     load(fileload,'x','y','x5','y5','conds','condsB');
%     % Keep infromation
%     ids = strfind(snames(sj).name, '_DS');
%     dur(sj).name = snames(sj).name(1:ids-1);
%     lens = y - x;
%     dur(sj).nTrial = length(x);
%     dur(sj).nTrialTarg = length(x5{1}) + length(x5{2});
%     dur(sj).nTrialFoil = length(x5{3}) + length(x5{4});
%     dur(sj).avDur = mean(lens);
%     dur(sj).avDurTarg = mean([y5{1}-x5{1};y5{2}-x5{2}]);
%     dur(sj).avDurFoil = mean([y5{3}-x5{3};y5{4}-x5{4}]);
%     dur(sj).modLikeli = likehood;
%     dur(sj).fs = fs;
%     
%     expLoc = zeros(dur(sj).nTrial,nBump);
%     expLocTarg = [];
%     expLocFoil = [];
%     lensTarg = [];
%     lensFoil = [];
%     for tr = 1:dur(sj).nTrial,
%         % Expected bump loacations
%         expLoc(tr,:) = [1:lens(tr)] * reshape(eventprobs(1:lens(tr),tr,:),lens(tr),nBump);
%         if conds(tr) == 1 | conds(tr) == 2,
%             expLocTarg(end+1,:) = expLoc(tr,:);
%             lensTarg(end+1) = lens(tr);
%         elseif conds(tr) == 3 | conds(tr) == 4,
%             expLocFoil(end+1,:) = expLoc(tr,:);
%             lensFoil(end+1) = lens(tr);
%         else
%             error('Error: There is an unexpected sort of condition')
%         end
%     end
%     clear eventprobs conds condsB x5 y5 x y
%     % Stages durations
%         % all Trial
%     duration = zeros(dur(sj).nTrial,nBump+1);
%     duration = expLoc(:,1);
%     duration(:,2:nBump) = expLoc(:,2:end) - expLoc(:,1:nBump-1);
%     duration(:,nBump+1) = lens - expLoc(:,nBump);
%     dur(sj).stag1 = mean(duration(:,1));
%     dur(sj).stag2 = mean(duration(:,2));
%     dur(sj).stag3 = mean(duration(:,3));
%     dur(sj).stag4 = mean(duration(:,4));
%     dur(sj).stag5 = mean(duration(:,5));
%     dur(sj).stag6 = mean(duration(:,6));
%     if nBump == 6, dur(sj).stag7 = mean(duration(:,7)); end
%         % Per conditions: Target
%     duration = zeros(dur(sj).nTrialTarg,nBump+1);
%     duration = expLocTarg(:,1);
%     duration(:,2:nBump) = expLocTarg(:,2:end) - expLocTarg(:,1:nBump-1);
%     duration(:,nBump+1) = lensTarg' - expLocTarg(:,nBump);
%     dur(sj).stag1Targ = mean(duration(:,1));
%     dur(sj).stag2Targ = mean(duration(:,2));
%     dur(sj).stag3Targ = mean(duration(:,3));
%     dur(sj).stag4Targ = mean(duration(:,4));
%     dur(sj).stag5Targ = mean(duration(:,5));
%     dur(sj).stag6Targ = mean(duration(:,6));
%     if nBump == 6, dur(sj).stag7Targ  = mean(duration(:,7)); end
%         % Per conditions: Foil
%     duration = zeros(dur(sj).nTrialFoil,nBump+1);
%     duration = expLocFoil(:,1);
%     duration(:,2:nBump) = expLocFoil(:,2:end) - expLocFoil(:,1:nBump-1);
%     duration(:,nBump+1) = lensFoil' - expLocFoil(:,nBump);
%     dur(sj).stag1Foil = mean(duration(:,1));
%     dur(sj).stag2Foil = mean(duration(:,2));
%     dur(sj).stag3Foil = mean(duration(:,3));
%     dur(sj).stag4Foil = mean(duration(:,4));
%     dur(sj).stag5Foil = mean(duration(:,5));
%     dur(sj).stag6Foil = mean(duration(:,6));
%     if nBump == 6, dur(sj).stag7Foil  = mean(duration(:,7)); end
%     clear expLoc expLocTarg expLocFoil
% end
% % Global average results
% dur(Nsj+1).name = 'GlobalMean';
% dur(Nsj+1).nTrial = mean([dur.nTrial]);
% dur(Nsj+1).nTrialTarg = mean([dur.nTrialTarg]);
% dur(Nsj+1).nTrialFoil = mean([dur.nTrialFoil]);
% dur(Nsj+1).avDur = mean([dur.avDur]);
% dur(Nsj+1).avDurTarg = mean([dur.avDurTarg]);
% dur(Nsj+1).avDurFoil = mean([dur.avDurFoil]);
% dur(Nsj+1).stag1 = mean([dur.stag1]);
% dur(Nsj+1).stag2 = mean([dur.stag2]);
% dur(Nsj+1).stag3 = mean([dur.stag3]);
% dur(Nsj+1).stag4 = mean([dur.stag4]);
% dur(Nsj+1).stag5 = mean([dur.stag5]);
% dur(Nsj+1).stag6 = mean([dur.stag6]);
% if nBump == 6, dur(Nsj+1).stag7 = mean([dur.stag7]); end
% dur(Nsj+1).stag1Targ = mean([dur.stag1Targ]);
% dur(Nsj+1).stag2Targ = mean([dur.stag2Targ]);
% dur(Nsj+1).stag3Targ = mean([dur.stag3Targ]);
% dur(Nsj+1).stag4Targ = mean([dur.stag4Targ]);
% dur(Nsj+1).stag5Targ = mean([dur.stag5Targ]);
% dur(Nsj+1).stag6Targ = mean([dur.stag6Targ]);
% if nBump == 6, dur(Nsj+1).stag7Targ = mean([dur.stag7Targ]); end
% dur(Nsj+1).stag1Foil = mean([dur.stag1Foil]);
% dur(Nsj+1).stag2Foil = mean([dur.stag2Foil]);
% dur(Nsj+1).stag3Foil = mean([dur.stag3Foil]);
% dur(Nsj+1).stag4Foil = mean([dur.stag4Foil]);
% dur(Nsj+1).stag5Foil = mean([dur.stag5Foil]);
% dur(Nsj+1).stag6Foil = mean([dur.stag6Foil]);
% if nBump == 6, dur(Nsj+1).stag7Foil = mean([dur.stag7Foil]); end

%%%% -------------------------------------------
% Repeat stage duration operations for HSMM computed from all subjects concatenated
%%% ---------------------------------
% pathdata = '/home/oscar/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/';
% filecode = ['AllSj100HSMM' num2str(nBump) 'Bout235.mat'];
% 
% pathvars = '/home/oscar/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/DS100epoch235allSj.mat';
% 
% % Load data sets
% fileload = [pathdata filecode];
% load(fileload);
% clear bumpMag gammPara
% load(pathvars,'x','y','x5','y5','conds');
% % Keep infromation
% name = ['allConct' num2str(nBump) 'B'];
% dur(Nsj+2).name = name;
% lens = y - x;
% dur(Nsj+2).nTrial = length(x);
% dur(Nsj+2).nTrialTarg = length(x5{1}) + length(x5{2});
% dur(Nsj+2).nTrialFoil = length(x5{3}) + length(x5{4});
% dur(Nsj+2).avDur = mean(lens);
% dur(Nsj+2).avDurTarg = mean([y5{1}-x5{1};y5{2}-x5{2}]);
% dur(Nsj+2).avDurFoil = mean([y5{3}-x5{3};y5{4}-x5{4}]);
% dur(Nsj+2).modLikeli = likehood;
% dur(Nsj+2).fs = fs;
% 
% expLoc = zeros(dur(Nsj+2).nTrial,nBump);
% expLocTarg = [];
% expLocFoil = [];
% lensTarg = [];
% lensFoil = [];
% for tr = 1:dur(Nsj+2).nTrial,
%     % Expected bump loacations
%     expLoc(tr,:) = [1:lens(tr)] * reshape(eventprobs(1:lens(tr),tr,:),lens(tr),nBump);
%     if conds(tr) == 1 | conds(tr) == 2,
%         expLocTarg(end+1,:) = expLoc(tr,:);
%         lensTarg(end+1) = lens(tr);
%     elseif conds(tr) == 3 | conds(tr) == 4,
%         expLocFoil(end+1,:) = expLoc(tr,:);
%         lensFoil(end+1) = lens(tr);
%     else
%         error('Error: There is an unexpected sort of condition')
%     end
% end
% clear eventprobs conds condsB x5 y5 x y
% % Stages durations
%     % all Trial
% duration = zeros(dur(Nsj+2).nTrial,nBump+1);
% duration = expLoc(:,1);
% duration(:,2:nBump) = expLoc(:,2:end) - expLoc(:,1:nBump-1);
% duration(:,nBump+1) = lens - expLoc(:,nBump);
% dur(Nsj+2).stag1 = mean(duration(:,1));
% dur(Nsj+2).stag2 = mean(duration(:,2));
% dur(Nsj+2).stag3 = mean(duration(:,3));
% dur(Nsj+2).stag4 = mean(duration(:,4));
% dur(Nsj+2).stag5 = mean(duration(:,5));
% dur(Nsj+2).stag6 = mean(duration(:,6));
% if nBump == 6, dur(Nsj+2).stag7 = mean(duration(:,7)); end
%     % Per conditions: Target
% duration = zeros(dur(Nsj+2).nTrialTarg,nBump+1);
% duration = expLocTarg(:,1);
% duration(:,2:nBump) = expLocTarg(:,2:end) - expLocTarg(:,1:nBump-1);
% duration(:,nBump+1) = lensTarg' - expLocTarg(:,nBump);
% dur(Nsj+2).stag1Targ = mean(duration(:,1));
% dur(Nsj+2).stag2Targ = mean(duration(:,2));
% dur(Nsj+2).stag3Targ = mean(duration(:,3));
% dur(Nsj+2).stag4Targ = mean(duration(:,4));
% dur(Nsj+2).stag5Targ = mean(duration(:,5));
% dur(Nsj+2).stag6Targ = mean(duration(:,6));
% if nBump == 6, dur(Nsj+2).stag7Targ  = mean(duration(:,7)); end
%     % Per conditions: Foil
% duration = zeros(dur(Nsj+2).nTrialFoil,nBump+1);
% duration = expLocFoil(:,1);
% duration(:,2:nBump) = expLocFoil(:,2:end) - expLocFoil(:,1:nBump-1);
% duration(:,nBump+1) = lensFoil' - expLocFoil(:,nBump);
% dur(Nsj+2).stag1Foil = mean(duration(:,1));
% dur(Nsj+2).stag2Foil = mean(duration(:,2));
% dur(Nsj+2).stag3Foil = mean(duration(:,3));
% dur(Nsj+2).stag4Foil = mean(duration(:,4));
% dur(Nsj+2).stag5Foil = mean(duration(:,5));
% dur(Nsj+2).stag6Foil = mean(duration(:,6));
% if nBump == 6, dur(Nsj+2).stag7Foil  = mean(duration(:,7)); end
% clear expLoc expLocTarg expLocFoil
% 
% dur(Nsj+3).name = [num2str(nBump) '_bumpModel'];
% 
% savename = ['satageTimeInfo' num2str(nBump) 'Bu.mat'];
% savepath = [pathdata savename];
% save(savepath, 'dur')

%%%% -------------------------------------------
% Do stage onset operations for HSMM computed from all subjects concatenated
%%% ---------------------------------
pathdata = '/home/oscar/Documents/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/';
filecode = ['AllSj100HSMM' num2str(nBump) 'Bout235.mat'];

pathvars = '/home/oscar/Documents/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/DS100epoch235allSj.mat';

% Load data sets
fileload = [pathdata filecode];
load(fileload);
clear bumpMag gammPara
load(pathvars,'x','y','x5','y5','conds');
lens = y - x;

expLoc = zeros(length(y),nBump);
for tr = 1:length(y)
    % Expected bump loacats
    expLoc(tr,:) = [1:lens(tr)] * reshape(eventprobs(1:lens(tr),tr,:),lens(tr),nBump);
end
clear eventprobs conds condsB x5 y5 x
% pass to milliseconds
expLoc = 10 * expLoc;
% Stages durations
% cumulative sum of time
dur(1).stag1 = mean(expLoc(:,1));
dur(1).stag2 = mean(expLoc(:,2));
dur(1).stag3 = mean(expLoc(:,3));
dur(1).stag4 = mean(expLoc(:,4));
dur(1).stag5 = mean(expLoc(:,5));
dur(1).avDur = 10 * mean(lens);

dur(4).stag1 = std([expLoc(:,1)]);
dur(4).stag2 = std([expLoc(:,2)]);
dur(4).stag3 = std([expLoc(:,3)]);
dur(4).stag4 = std([expLoc(:,4)]);
dur(4).stag5 = std([expLoc(:,5)]);

dur(5).stag1 = dur(4).stag1;
dur(5).stag2 = dur(4).stag2 - dur(4).stag1;
dur(5).stag3 = dur(4).stag3 - dur(4).stag2;
dur(5).stag4 = dur(4).stag4 - dur(4).stag3;
dur(5).stag5 = dur(4).stag5 - dur(4).stag4;

onsets = diff([repmat(0,length(lens),1), expLoc, 10*lens],1,2);

dur(2).stag1 = mean(onsets(:,1));
dur(2).stag2 = mean(onsets(:,2));
dur(2).stag3 = mean(onsets(:,3));
dur(2).stag4 = mean(onsets(:,4));
dur(2).stag5 = mean(onsets(:,5));
dur(2).stag6 = mean(onsets(:,6));

% standard deviation
dur(3).stag1 = std(onsets(:,1));
dur(3).stag2 = std(onsets(:,2));
dur(3).stag3 = std(onsets(:,3));
dur(3).stag4 = std(onsets(:,4));
dur(3).stag5 = std(onsets(:,5));
dur(3).stag6 = std(onsets(:,6));

dur(1).name = [num2str(nBump) '_meanCumONSETs'];
dur(2).name = [num2str(nBump) '_meanDurStags'];
dur(3).name = [num2str(nBump) '_stdDurStags'];
dur(4).name = [num2str(nBump) '_stdCumOnset'];
dur(5).name = [num2str(nBump) '_stdCumOnCorrect'];

savename = ['onsetStgTimeInfo' num2str(nBump) 'Bu.mat'];
savepath = [pathdata savename];
save(savepath, 'dur')
