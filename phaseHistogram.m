% makes a histogram of phases with the bump PDF and the pahse of the
% analytic signal at a frequency band. It works for one subject and one
% bump

clear all
close all
clc

nBump = 5;
nBins = 20;                        % number of bins in the phase histogram
edges = linspace(-pi, pi, nBins+1);  % Does nBins from -pi to pi Bin edges
width = abs(edges(1) - edges(2));        % width of bins
binVal = edges(1:end-1) - diff(edges)/2;    % value of each bin
nCh = 32;                           % number of channels
freqBand = 'theta';                 % frequency band to be analyzed

% prepare data paths to load
pathdata = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/';
varElec = [pathdata 'electrodes.mat'];
codePDF = ['*DS100HSMM' num2str(nBump) 'Bout235.mat'];
namesPDF = dir([pathdata codePDF]);
codeHilb = ['*Bads100.mat'];
namesHilb = dir([pathdata codeHilb]);
Nsj = length(namesPDF);

if length(namesPDF) ~= length(namesHilb), error('Unbalancend number of files per subject'), end


for sj = 1:Nsj
    % load datasets
    varPDF = [pathdata namesPDF(sj).name];
    varHilb =[pathdata namesHilb(sj).name];
    load(varElec)
    load(varPDF, 'eventprobs')
    data = load(varHilb, freqBand, 'x', 'y');
    x = data.x;
    y = data.y;
    % do pahse angle from analytic signal
    datAll = angle(data.(freqBand));
    clear data

    probump = [];
        % Concatenate trials vertically to remove zeros out of trial
    for tr = 1:length(y)
        yTr = y(tr) - x(tr) + 1;
        probump = vertcat(probump, squeeze(eventprobs(1:yTr,tr,:)));
    end
    clear eventprobs
    BinIdx = zeros(size(datAll));       % [All samples x channels]
        % assign to each sample a phase bin
    for ch = 1:nCh
        BinIdx(:,ch) = discretize(datAll(:,ch), edges);
    end

    phasePDF = zeros(nCh,nBins,nBump);
    for bu = 1:nBump
        for ch = 1:nCh
            % sums the prob of a bump at each phase bin for all saples
            phasePDF(ch,:,bu) = accumarray(BinIdx(:,ch),probump(:,bu));
        end
    end       
    % Normalization
    phasePDF = phasePDF ./ length(y);
    %clear dataPhase BinIdx
    
    % Expected phase from the phase PDF
    exPhaseZx = complex(zeros(nCh,nBump));   % [channels, subjects, model]
    for bu =1:nBump
        exPhaseZx(:,sj,bu) = exp(1i*binVal) * squeeze(phasePDF(:,:,bu))';
    end

%     % visual representation
%     fig = figure(sj);
%     titlen = ['Phase distribution of: ' freqBand ', subject: ' num2str(sj)];
%     colorB = {'b', 'r','c','m','g','k','y','b'};
%     colorB(nBump+1:end) = [];
%     legendB = {'Bump 1','Bump 2','Bump 3','Bump 4','Bump 5','Bump 6','Bump 7','Bump 8'};
%     legendB(nBump+1:end) = [];
%     ymax = max(max(max(phasePDF))) + 0.1 * max(max(max(phasePDF)));
%     plottopo(phasePDF, 'chanlocs', locElect, 'ylim', [0,ymax], 'title',titlen, 'colors', colorB, ...
%                         'ydir', 1,'legend', legendB, 'showleg','on')


        % Phase Histogram on polar form for all Channels
        side = 0.095;
        [~, chNames, Th, Rd] = readlocs(locElect);
        Th = pi/180*Th;                 % convert degrees to radians
        [yOr,xOr] = pol2cart(Th,Rd); % translate from polar to cart. coordinates
        % yOr = yOr + 0.5;
        % xOr = xOr + 0.5;
        yOr = yOr / (abs(max(yOr) - min(yOr)) + 3*side) + 0.5 + side/3;
        xOr = xOr / (abs(max(xOr) - min(xOr)) + 3*side) + 0.5 + side/3;
        position(:,1) = xOr - side;
        position(:,2) = yOr - side;
        position(:,3) = xOr + side;
        position(:,4) = yOr + side;
        % Test plot
        %plot(xOr,yOr,'o',position(:,1), position(:,2), '*',position(:,3), position(:,4), '+')
        colorbu = {'b','r','c','m','g'};
        fig = figure(sj)
        for i = 1:32
            subplot('Position', [position(i,1),position(i,2),side,side])
            set(0,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);
            %P = polaRad(edges, 0.15 * ones(size(edges)));
            %set(P, 'Visible', 'off')
            %hold on
            for bu = 1:5
                %plot([1,2])
                h = polaRad(edges,[phasePDF(i,:,bu),phasePDF(i,1,bu)]);
                title(chNames(i),'HorizontalAlignment', 'left')
                h.LineWidth = 2.5;
                h.Color = colorbu{bu};
                hold on
            end
        end

%                     % save figure to a .fig
    namesave = ['sj' num2str(sj) freqBand 'Modl' num2str(nBump) '.fig' ];
    savepath = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/histPhase/';
    tosave = [savepath namesave ];
    savefig(fig,tosave)
    clear fig
end
% close all