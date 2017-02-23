% plot bump's topologies per subjects an all subjects concatenated. 
% It allows to compare the estimated
% bumps tolologies for multiple subjectsa as well as all subejcts
% concatenated

clear all
close all
clc

compare = 0;            % whether to plot to Plots to copare topologies btwn same subject
allConc = 0;            % plot topologies for all PCs subjects concatenated
eachSj = 0;             % plot topologies for PCs of each subject
topoStds = 1;           % plot standard deviation on bumps topologies
nBump = 5;
sjPlot = [1:20];        % Choose which subjects' topologies will be ploted
%sjPlot = [3,4,5,6,9,10,12,13,14,19,20];
%sjPlot = [1];
maxTime = 300;          % maximum number of samples: time = samples * (1/fs)
nCh = 32;

if compare, allConc = 0; eachSj = 0; end

% prepare data paths to load
pathPCA = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/';
namesPCA = dir([pathPCA '*epochs235.mat']);
pathdata = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/';
varElec = [pathdata 'electrodes.mat'];
code= ['*DS100HSMM' num2str(nBump) 'Bout235.mat'];
namesBu = dir([pathdata code]);
varAll = [pathdata 'AllSj100HSMM' num2str(nBump) 'Bout235.mat'];
code = ['satageTimeInfo' num2str(nBump) 'Bu.mat'];
nameT = [pathdata code];

% load variables
load(varElec);
allSj = load(varAll, 'bumpMag');
load(nameT);
    % Get stages durations
if nBump == 5,
    durStagAll = [dur(22).stag1, dur(22).stag2, dur(22).stag3, dur(22).stag4, dur(22).stag5, dur(22).stag6];
elseif nBump == 6,
    durStagAll = [dur(22).stag1, dur(22).stag2, dur(22).stag3, dur(22).stag4, dur(22).stag5, dur(22).stag6, dur(22).stag6]; 
else
    error('Error: nBump should be 5 ot 6')
end
timeStagAll = cumsum(durStagAll);
timeStagAll = round(timeStagAll .* (maxTime / timeStagAll(end)));
% plotting variable
locStagAll = NaN(1,maxTime);
locStagAll(timeStagAll) = 1;
locStagAll(end) = NaN;

topoSj  = zeros(nBump, nCh, length(sjPlot));
topoAll  = zeros(nBump, nCh, length(sjPlot));

for sj = 1:length(sjPlot),
    % Get stages durations
    if nBump == 5,
        durStag = [dur(sj).stag1, dur(sj).stag2, dur(sj).stag3, dur(sj).stag4, dur(sj).stag5, dur(sj).stag6];
    elseif nBump == 6,
        durStag = [dur(sj).stag1, dur(sj).stag2, dur(sj).stag3, dur(sj).stag4, dur(sj).stag5, dur(sj).stag6, dur(sj).stag6]; 
    end
    timeStag = cumsum(durStag);
    timeStag = round(timeStag .* (maxTime / timeStag(end)));
    % plotting variable
    locStag = NaN(1,maxTime);
    locStag(timeStag) = 1;
    locStag(end) = NaN;
    % load and do variable
    pcs = load([pathPCA namesPCA(sjPlot(sj)).name], 'means10', 'coeff10', 'latent10');
    sjBump = load([pathdata namesBu(sjPlot(sj)).name], 'bumpMag');
    % reconstruction of mean bump topology
    bumpTopoSj = reconstruct(sjBump.bumpMag,pcs.coeff10,pcs.latent10,pcs.means10);
    bumpTopoAll = reconstruct(allSj.bumpMag,pcs.coeff10,pcs.latent10,pcs.means10);
    clear pcs sjBump
    %figure( 'Name', ['Subject: ' num2str(sjPlot(sj))]);
    %title(['Subject: ' num2str(sjPlot(sj))])
    for bu = 1:nBump,
        if compare
            subplot(3, nBump, bu)
            topoplot(bumpTopoSj(bu,:), locElect, 'maplimits', [-12 12],'style', 'map', 'electrodes','off', 'shading', 'interp');
            if bu == 1, title('HSMM on Single Subject'), end
            subplot(3, nBump, nBump+bu)
            topoplot(bumpTopoAll(bu,:), locElect, 'maplimits', [-12 12],'style', 'map', 'electrodes','off', 'shading', 'interp');
            if bu == 1, title('HSMM on All Subject Concatenated'), end
        elseif allConc
            subplot(1, nBump, bu)
            topoplot(bumpTopoAll(bu,:), locElect, 'maplimits', [-12 12],'style', 'map', 'electrodes','off', 'shading', 'interp');
        elseif eachSj
            subplot(1, nBump, bu)
            topoplot(bumpTopoSj(bu,:), locElect, 'maplimits', [-12 12],'style', 'map', 'electrodes','off', 'shading', 'interp');
        end 
    end
    if compare,
        subplot(3, nBump, [nBump+bu+1:nBump+2*bu])
        plot(locStag, '^','LineWidth',2)
        axis([0 maxTime 0.5 1.5])
        %legend('SingleSj')
        hold on
        plot(locStagAll, 'v','LineWidth',2)
        legend('SingleSj','AllConcat')
        title('Temporal location of bumps (normalized)')
    end
    topoSj(:,:,sj) = bumpTopoSj;
    topoAll(:,:,sj) = bumpTopoAll;
end

% standard deiation of bups topologies of individual subjects
if topoStds,
    stdTopoSj(:,:,1) = mean(topoSj,3) + std(topoSj,[],3);
    stdTopoSj(:,:,2) = mean(topoSj,3) - std(topoSj,[],3);
    stdTopoAll(:,:,1) = mean(topoAll,3) + std(topoAll,[],3);
    stdTopoAll(:,:,2) = mean(topoAll,3) - std(topoAll,[],3);

    figure('Name', ['Deviation on bumps ' num2str(nBump) ' topology with a HSMM per subject']);
    for bu = 1:nBump,
        subplot(3, nBump, bu)
        title('mean + std')
        topoplot(squeeze(stdTopoSj(bu,:,1)), locElect, 'maplimits', [-12 12],'style', 'map', 'electrodes','off', 'shading', 'interp');
        subplot(3, nBump, bu+nBump)
        title('mean')
        topoplot(squeeze(mean(topoSj(bu,:,:),3)), locElect, 'maplimits', [-12 12],'style', 'map', 'electrodes','off', 'shading', 'interp');
        subplot(3, nBump, bu+2*nBump)
        title('mean - std')
        topoplot(squeeze(stdTopoSj(bu,:,2)), locElect, 'maplimits', [-12 12],'style', 'map', 'electrodes','off', 'shading', 'interp');
    end
    cbar('horiz',0,[-12 12]);
    figure('Name', ['Deviation on bumps ' num2str(nBump) ' topology with a HSMM for all subject']);
    for bu = 1:nBump,
        subplot(3, nBump, bu)
        title('mean + std')
        topoplot(squeeze(stdTopoAll(bu,:,1)), locElect, 'maplimits', [-12 12],'style', 'map', 'electrodes','off', 'shading', 'interp');
        subplot(3, nBump, bu+nBump)
        title('mean')
        topoplot(squeeze(mean(topoAll(bu,:,:),3)) , locElect, 'maplimits', [-12 12],'style', 'map', 'electrodes','off', 'shading', 'interp');
        subplot(3, nBump, bu+2*nBump)
        title('mean - std')
        topoplot(squeeze(stdTopoAll(bu,:,2)), locElect, 'maplimits', [-12 12],'style', 'map', 'electrodes','off', 'shading', 'interp');
    end
    cbar('horiz',[0],[-12 12]);
end


