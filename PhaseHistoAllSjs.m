% makes a histogram of phases with the bump PDF and the pahse of the
% analytic signal at a frequency band. It works fro all subjects and one
% bump

clear all
close all
clc

bumpModl = 5;
nBins = 20;                        % number of bins in the phase histogram
edges = linspace(-pi, pi, nBins+1);  % Does nBins from -pi to pi Bin edges
width = abs(edges(1) - edges(2));        % width of bins
nCh = 32;                           % number of channels
%sj = 1;                             % subject to be analyzed
nSj = 20;                           % number of subjects
freqBand = 'thetaH.mat';            % frequency band to be analyzed

pathPDF = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125bp05_58/events/ForBumps_OnRes_no20/';
pathOther = [pathPDF 'varForBumpsOn_Res12575.mat'];
pathElec = [pathPDF 'electrodes.mat'];
pathPDF = [pathPDF 'stages125_0558_' num2str(bumpModl) '_Bu.mat'];


pathFreq = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125bp05_58/events/ForBumps_OnRes_no20/LocBandsHilbert/';
pathFreq = [pathFreq freqBand];

load(pathFreq)
load(pathPDF, 'eventprob20i')
load(pathOther, 'x20', 'y20', 'subjects')
load(pathElec)

% do pahse angle from analytic signal
dataAll = angle(dataAll);
% check that there is not values out of range [-pi, pi] due to rounding
dataAll(dataAll>pi) = pi - width^2;
dataAll(dataAll<-pi) = -pi + width^2;

for sj = 1:nSj
    probump = [];
        % Concatenate trials vertically to remove zeros out of trial
        
    %%% all subjects HSMM computed at the same time
%     eventprobSj = eventprobs2(:,find(subjects == sj),:);
%     for tr = 1:length(y20{sj})
%         y = y20{sj}(tr) - x20{sj}(tr) + 1;
%         probump = vertcat(probump, squeeze(eventprobSj(1:y,tr,:)));
%     end
    %clear eventprobSj
    %%% Independently HSMM computed subects
    for tr = 1:length(y20{sj})
        y = y20{sj}(tr) - x20{sj}(tr) + 1;
        probump = vertcat(probump, squeeze(eventprob20i{sj}(1:y,tr,:)));
    end

    binIdx = zeros(size(dataAll(x20{sj}(1):y20{sj}(end),:)));       % [All samples x channels]
        % assign to each sample a phase bin
    for ch = 1:nCh
        binIdx(:,ch) = discretize(dataAll(x20{sj}(1):y20{sj}(end),ch), edges);
    end

    phasePDF = zeros(nCh,nBins,bumpModl);
    for bu = 1:bumpModl
        for ch = 1:nCh
            % sums the prob of a bump at each phase bin for all saples
            phasePDF(ch,:,bu) = accumarray(binIdx(:,ch),probump(:,bu));
            %phasePDF(ch,:,bu) = output ./ sum(output);
        end
    end
    % normalization by number of trials
    phasePDF = phasePDF ./ length(y20{sj});
    
    %clear binIdx probump

    % visual representation
    fig = figure(1);
    titlen = ['Phase distribution of: ' freqBand(1:end-5) ', subject: ' num2str(sj)];
    colorB = {'b', 'r','c','m','g','k','y','b'};
    colorB(bumpModl+1:end) = [];
    legendB = {'Bump 1','Bump 2','Bump 3','Bump 4','Bump 5','Bump 6','Bump 7','Bump 8'};
    legendB(bumpModl+1:end) = [];
    ymax = max(max(max(phasePDF))) + 0.1 * max(max(max(phasePDF)));
    plottopo(phasePDF, 'chanlocs', locElect, 'ylim', [0,ymax], 'title',titlen, 'colors', colorB, ...
                        'ydir', 1,'legend', legendB, 'showleg','on')

    % save figure to a .fig
    namesave = ['sj' num2str(sj) freqBand(1:end-5) 'Modl' num2str(bumpModl) 'HSMMind.fig' ];
    savepath = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125bp05_58/events/ForBumps_OnRes_no20/histPhase/';
    tosave = [savepath namesave ];
    savefig(fig,tosave)
    clear fig
    close all
end