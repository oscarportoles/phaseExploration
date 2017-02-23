% Estimate the temporal duration of each stage (without sidtinguishing conditions)
% with HSMM computed from all subjects concatenated. Condition 'New' is not
% in the data, It was not included in the HSMM

clear all
close all
clc

nBump = 6;
fs = 100; % sampling frequency

pathdata = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/';
filecode = ['AllSj100HSMM' num2str(nBump) 'Bout235.mat'];

pathvars = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/DS100epoch235allSj.mat';

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


% Load data sets
fileload = [pathdata filecode];
load(fileload);
clear bumpMag gammPara
load(pathvars,'x','y','x5','y5','conds');
% Keep infromation
name = ['allConct' num2str(nBump) 'B'];
dur(1).name = name;
lens = y - x;
dur(1).nTrial = length(x);
dur(1).nTrialTarg = length(x5{1}) + length(x5{2});
dur(1).nTrialFoil = length(x5{3}) + length(x5{4});
dur(1).avDur = mean(lens);
dur(1).avDurTarg = mean([y5{1}-x5{1};y5{2}-x5{2}]);
dur(1).avDurFoil = mean([y5{3}-x5{3};y5{4}-x5{4}]);
dur(1).modLikeli = likehood;
dur(1).fs = fs;

expLoc = zeros(dur(1).nTrial,nBump);
expLocTarg = [];
expLocFoil = [];
lensTarg = [];
lensFoil = [];
for tr = 1:dur(1).nTrial,
    % Expected bump loacations
    expLoc(tr,:) = [1:lens(tr)] * reshape(eventprobs(1:lens(tr),tr,:),lens(tr),nBump);
    if conds(tr) == 1 | conds(tr) == 2,
        expLocTarg(end+1,:) = expLoc(tr,:);
        lensTarg(end+1) = lens(tr);
    elseif conds(tr) == 3 | conds(tr) == 4,
        expLocFoil(end+1,:) = expLoc(tr,:);
        lensFoil(end+1) = lens(tr);
    else
        error('Error: There is an unexpected sort of condition')
    end
end
clear eventprobs conds condsB x5 y5 x y
% Stages durations
    % all Trial
duration = zeros(dur(1).nTrial,nBump+1);
duration = expLoc(:,1);
duration(:,2:nBump) = expLoc(:,2:end) - expLoc(:,1:nBump-1);
duration(:,nBump+1) = lens - expLoc(:,nBump);
dur(1).stag1 = mean(duration(:,1));
dur(1).stag2 = mean(duration(:,2));
dur(1).stag3 = mean(duration(:,3));
dur(1).stag4 = mean(duration(:,4));
dur(1).stag5 = mean(duration(:,5));
dur(1).stag6 = mean(duration(:,6));
if nBump == 6, dur(1).stag7 = mean(duration(:,7)); end
    % Per conditions: Target
duration = zeros(dur(1).nTrialTarg,nBump+1);
duration = expLocTarg(:,1);
duration(:,2:nBump) = expLocTarg(:,2:end) - expLocTarg(:,1:nBump-1);
duration(:,nBump+1) = lensTarg' - expLocTarg(:,nBump);
dur(1).stag1Targ = mean(duration(:,1));
dur(1).stag2Targ = mean(duration(:,2));
dur(1).stag3Targ = mean(duration(:,3));
dur(1).stag4Targ = mean(duration(:,4));
dur(1).stag5Targ = mean(duration(:,5));
dur(1).stag6Targ = mean(duration(:,6));
if nBump == 6, dur(1).stag7Targ  = mean(duration(:,7)); end
    % Per conditions: Foil
duration = zeros(dur(1).nTrialFoil,nBump+1);
duration = expLocFoil(:,1);
duration(:,2:nBump) = expLocFoil(:,2:end) - expLocFoil(:,1:nBump-1);
duration(:,nBump+1) = lensFoil' - expLocFoil(:,nBump);
dur(1).stag1Foil = mean(duration(:,1));
dur(1).stag2Foil = mean(duration(:,2));
dur(1).stag3Foil = mean(duration(:,3));
dur(1).stag4Foil = mean(duration(:,4));
dur(1).stag5Foil = mean(duration(:,5));
dur(1).stag6Foil = mean(duration(:,6));
if nBump == 6, dur(1).stag7Foil  = mean(duration(:,7)); end
clear expLoc expLocTarg expLocFoil


savename = ['satageTimeInfo' num2str(nBump) 'BuAllConct.mat'];
savepath = [pathdata savename];
save(savepath, 'dur')
    
            
