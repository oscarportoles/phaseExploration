% Plot the scalp topografies of the newly calculated bumps, same as Figure
% 5 in the paper from Qiong 2016
clear all
close all
clc

% variables to load from path:
% '/MATLAB/Data_JB_AssoRec/DS125bp05_30/events/ForBumps_OnRes_no20'

% load electrodes.mat
% load varForBumpsOn_Res75dtr.mat subjects coeff10 latent10 data x y
% %load forLOOCV3075N.mat mags8
% means10 =  mean(data,1); % Grand mean accross samples, trials and subjects of the 32 electrode [1 x Ch]
% clear data
% flagTime = 1; % plots temporal location of bumps


% variables to load from path: 
% '/MATLAB/newbumps'
load starter.mat cl means10 coeff10 latent10 x y
% load fit8.mat mags8
locElect = cl;
flagTime = 1; % Not plot bumps temporal lication

%%% Core

d = y - x;
meanTime = floor(mean(d));
electrodes8 = cell(1,8);
timeBump = cell(1, 8); % mean time location of the bumpsclear all
samplBump = cell(1, 8); % samples with highest probability of having a bump
mags8 = cell(1,8);
bT = ones(8, meanTime);
t = 1:meanTime;

% create bumps topografy and other variable, the values are calculated as
% mean for all subjects and all trials
for i = 1:8
    if flagTime,
        filename = ['iniCondOutQmanyGF_' num2str(i) '_Bu.mat'];
        probs = load(filename,'eventprobs2', 'bumpMag2'); %
        mags8{i} = probs.bumpMag2;
        timeBump{i} = squeeze(mean(probs.eventprobs2, 2)); % mean probability a bump in a sample
        [~, mIdx] = max(timeBump{i}); % median?
        samplBump{i} = mIdx;
        
        bT(i,:) = bT(i,:) .* (9 - i); % for plotting the bumps time points
        electrodes8{i} = reconstruct(mags8{i},coeff10,latent10,means10); % rconstruction of mean bump topology
    end
end


% Plot of bumps' scalp maps
for i = 1:8,
    [rowB,colB] = size(mags8{i});
    for bump = 1:colB,
        figure(i)
        subplot(1, colB, bump)
        topoplot(electrodes8{i}(bump,:), locElect, 'maplimits', [-12 12],'style', 'map', 'electrodes','off', 'shading', 'interp');
    
    end
    if flagTime
        figure(100)
        plot(t, bT(i,:),'-b', t(samplBump{i}),bT(i,samplBump{i}),'ro')
        hold on
    end
end

lkmanyGF = zeros(1,8);
lkmany = zeros(1,8);
lk10gF = zeros(1,8);
lkN10 = zeros(1,8);
lkN10gF = zeros(1,8);
tl = 1:8;
for i = 1:8
    filename = ['iniCondOutQmanyGF_' num2str(i) '_Bu.mat'];
    load(filename,'likehood2');
    lkmanyGF(i) = likehood2;
    clear likehood2
    filename = ['iniCondOutQmany_' num2str(i) '_Bu.mat'];
    load(filename,'likehood2');
    lkmany(i) = likehood2;
    clear likehood2
    filename = ['iniCondOutQ10gF_' num2str(i) '_Bu.mat'];
    load(filename,'likehood2');
    lk10gF(i) = likehood2;
    clear likehood2
    filename = ['iniCondOutQn10_' num2str(i) '_Bu.mat'];
    load(filename,'likehood2');
    lkN10(i) = likehood2;
    clear likehood2
    filename = ['iniCondOutQn10gF_' num2str(i) '_Bu.mat'];
    load(filename,'likehood2');
    lkN10gF(i) = likehood2;
    clear likehood2
end
plot(tl, lkmanyGF, 'o-',tl, lkmany, 'o-',tl, lk10gF, 'o-',tl, lkN10, 'o-',tl, lkN10gF, 'o-')
