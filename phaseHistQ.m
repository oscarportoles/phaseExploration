% % % TEst with Qiong's data% 
% % % makes a histogram of phases with the bump PDF and the pahse of the
% % % analytic signal at a frequency band. It works fro all subjects and one
% % % bump
% 
% clear all
% close all
% clc
% 
% bumpModl = 5;
% nBins = 20;                        % number of bins in the phase histogram
% edges = linspace(-pi, pi, nBins+1);  % Does nBins from -pi to pi Bin edges
% width = abs(edges(1) - edges(2));        % width of bins
% nCh = 32;                           % number of channels
% %sj = 1;                             % subject to be analyzed
% nSj = 20;                           % number of subjects
% fs = 100;                           % samping ration
% locutoff = 4;
% hicutoff = 9;
% 
% path = '/Users/lab/Documents/MATLAB/newbumps/';
% pathStarter = [path 'starter.mat'];
% pathProbs = [path 'fit8.mat'];
% 
% pathElec = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125bp05_58/events/ForBumps_OnRes_no20/electrodes.mat';
% 
% load(pathElec)
% 
% load(pathStarter, 'x20', 'y20', 'subjects', 'subjectsF', 'x', 'y', 'compressed10')
% 
% % filter the data
% options = { fs, locutoff, hicutoff};
% data = iirfilt( compressed10', options{:}); % data(ch, pnts)
% clear compressed10
% 
% load(pathProbs, 'eventprobs5')
% 
% % Hilbert transform
% dataH = [];
% probump = [];
% for tr = 1:length(y20{1})
%     dataout = angle(hilbert(data(:,x20{1}(tr):y20{1}(tr))'));  % hilbert() computes H.T. columnswise
%     dataH = vertcat(dataH, dataout);           % dataout [samples, ch]
%     yy = y20{1}(tr) - x20{1}(tr) + 1;
%     probump = vertcat(probump, squeeze(eventprobs5(1:yy,tr,:)));
% end
% clear data eventprobs5
% 
% dataH(dataH>pi) = pi - width^2;
% dataH(dataH<-pi) = -pi + width^2;
exPhaseZx = complex(zeros(nCh,nSj,bumpModl));   % [channels, subjects, model]
exPhaseZxm = complex(zeros(nCh,nSj,bumpModl));   % [channels, subjects, model]
binVal = edges(1:end-1) - diff(edges)/2;
    % assign to each sample a phase bin
for sj =1:1
    binIdx = zeros(size(dataH(find(subjectsF == sj)),1),nCh);       % [All samples x channels]
    for ch = 1:nCh
        binIdx(:,ch) = discretize(dataH(find(subjectsF == sj),ch), edges);
    end
    
    phasePDF = zeros(nCh,nBins,bumpModl);
    phasePDFm = zeros(nCh,nBins,bumpModl);
    
    for bu = 1:bumpModl
        for ch = 1:nCh
            % sums the prob of a bump at each phase bin for all saples
            phasePDF(ch,:,bu) = accumarray(binIdx(:,ch),probump(find(subjectsF == sj),bu));
            phasePDFm(ch,:,bu) = accumarray(binIdx(:,ch),probump(find(subjectsF == sj),bu),[],@mean);
        end
    end
    % normalization by number of trials
    phasePDF = phasePDF ./ length(y20{sj});
    %phasePDFm = phasePDFm ./ length(y20{sj});
    
    for bu =1:bumpModl
        % Expected phase
        exPhaseZx(:,sj,bu) = exp(1i*binVal) * squeeze(phasePDF(:,:,bu))';
        exPhaseZxm(:,sj,bu) = exp(1i*binVal) * squeeze(phasePDFm(:,:,bu))';
    end
    
    freqBand = 'Theta';
    fig = figure(1);
    titlen = ['Phase distribution of: ' freqBand ', subject: ' num2str(sj), ' @sum'];
    colorB = {'b', 'r','c','m','g','k','y','b'};
    colorB(bumpModl+1:end) = [];
    legendB = {'Bump 1','Bump 2','Bump 3','Bump 4','Bump 5','Bump 6','Bump 7','Bump 8'};
    legendB(bumpModl+1:end) = [];
    ymax = max(max(max(phasePDF))) + 0.1 * max(max(max(phasePDF)));
    plottopo(phasePDF, 'chanlocs', locElect, 'ylim', [0,ymax], 'title',titlen, 'colors', colorB, ...
                        'ydir', 1,'legend', legendB, 'showleg','on')
                    
    fig = figure(2);
    titlen = ['Phase distribution of: ' freqBand ', subject: ' num2str(sj), ' @mean'];
    colorB = {'b', 'r','c','m','g','k','y','b'};
    colorB(bumpModl+1:end) = [];
    legendB = {'Bump 1','Bump 2','Bump 3','Bump 4','Bump 5','Bump 6','Bump 7','Bump 8'};
    legendB(bumpModl+1:end) = [];
    ymax = max(max(max(phasePDFm))) + 0.1 * max(max(max(phasePDFm)));
    plottopo(phasePDFm, 'chanlocs', locElect, 'ylim', [0,ymax], 'title',titlen, 'colors', colorB, ...
                        'ydir', 1,'legend', legendB, 'showleg','on')

    % save figure to a .fig
%     namesave = ['sj' num2str(sj) freqBand 'Modl' num2str(bumpModl) 'HSMMallQ.fig' ];
%     savepath = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125bp05_58/events/ForBumps_OnRes_no20/histPhase/';
%     tosave = [savepath namesave ];
%     savefig(fig,tosave)
%     clear fig
%     close all
end
% 
% ISPC = squeeze(abs(mean(exp(1i*angle(exPhaseZx)),2)));       % Inter-Subject Phase Clustering
% ISPCz = nSj * ISPC.^2;                              % transformation for few subjects
% p = 0.05 / (nCh*bumpModl);
% pISPC = exp(sqrt(1 + 4*nSj + 4*(nSj^2 - (nSj*ISPC).^2))-(1+2*nSj));
% % Inter-Subject Phase Angle ISPA
% ISPA = squeeze(angle(mean(exp(1i*angle(exPhaseZx)),2)));
% ISPCri = sqrt(-log(p)/nSj);