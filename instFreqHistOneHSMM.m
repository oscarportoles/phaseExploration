% makes a histogram of phases with the bump PDF and the pahse of the
% analytic signal at a frequency band. It uses the results from a HSMM
% computed from the PCs of all subjects concatenated

clear all
close all
clc

nBump = 5;                          % number of bumps in HSMM
plotHist = 0;                       % plot head histogram per subject 
limFr = [2, 11];                    % extreme bins [Hz]   
nBins = 36;                         % number of bins in the phase histogram
edges = linspace(limFr(1), limFr(2), nBins+1);  % Does nBins, Theta band was filtered from 4 to 9 Hz
width = abs(edges(1) - edges(2));        % width of bins
binVal = edges(1:end-1) + diff(edges)/2;    % value of each bin
nCh = 32;                           % number of channels
freqBand = 'theta';                 % frequency band to be analyzed
nRep = 400;                         % number of repetition to build a random PDF
fs = 100;                           % sampling frequency [Hz]              
n_order = 10;                       % median filter parameters (Cohen X et.al. FReq sliding...)
orders = round(linspace(0.02*fs,0.4*fs,n_order)); % recommended: 10 steps between 10 (here 20ms becoause of fs) and 400 ms
orders = floor((orders-1)/2); 


% prepare data paths to load
pathdata = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/';
%varElec = [pathdata 'electrodes.mat'];
namePDF = ['AllSj100HSMM' num2str(nBump) 'Bout235.mat'];
codeHilb = ['*Bads100.mat'];
namesHilb = dir([pathdata codeHilb]);
nSj = length(namesHilb);

% test that the subject are sorted equaly so loaded datasets agree
pathdataTest = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/';
snamesTest = dir([pathdataTest '*epochs235.mat']);
if length(snamesTest) ~= length(namesHilb), error('Unbalancend number of files per subject'), end
for sj = 1:nSj
    idTest = strfind(snamesTest(sj).name ,'_');
    idHil = strfind(namesHilb(sj).name ,'_');
    if namesHilb(sj).name(1:idHil) ~= snamesTest(sj).name(1:idTest)
        error('Error: Subjects are not the same')
    end
end
clear pathdataTest snamesTest idHil idTest codeHilb

% load variables
load([pathdata 'electrodes.mat'])
load([pathdata namePDF], 'eventprobs')

info = ['Instantaneous frequency histogram, Bump Prefered frequency of a HSMM' ...
        ' calculated for all subject concatenated.Fs:100Hz,Band: ' freqBand];

pValMi = zeros(nCh,nBump,nSj);          % p-values of Modulation index
miRa = zeros(nCh,nBump,nRep,nSj);       % MI of randomized trials
iFrePDF = zeros(nCh,nBins,nBump,nSj);    % Phase probability distribution
exIfe = zeros(nCh,nBump,nSj);         % expected complex value of pahse PD
mi = zeros(nCh,nBump,nSj);              % MI of emprical data

for sj = 1:nSj
    % load datasets
    %varHilb =[pathdata namesHilb(sj).name];
    data = load([pathdata namesHilb(sj).name], freqBand, 'x', 'y');
    x = data.x;
    y = data.y;
    % compute instantaneous requency from analytic signal
    datAll = angle(data.(freqBand)); 
    clear data
    dataIfre = [];
    probump = [];
    for tr = 1:length(y)    % iterate trials
        yTr = y(tr) - x(tr); % + 1; % One sample less to compensate for one-sample drop on 'diff()'
        % Concatenate trials vertically to remove zeros out of trial
        probump = vertcat(probump, squeeze(eventprobs(1:yTr,tr,:)));
        % Instantaneous frequency
        instFr = diff(unwrap(datAll(x(tr):y(tr),:)));
        % Add mirrowed tails at the edges of the trial to avoid edge artifacts
        sMirr = orders(end) + 1;        % number of samples to mirrow
        instFr = [flipud(instFr(1:sMirr,:)); instFr; flipud(instFr(end-sMirr:end,:))];
        % Compute median filtered instantaneous frequency (Cohen X 2014 Freq sliding... 2014)
        phasedmed = zeros(length(orders),size(instFr,1),nCh);
        for oi=1:n_order
            for ti=1:size(instFr,1)
                for ch =1:nCh
                    temp = sort(instFr(max(ti-orders(oi),1):min(ti+orders(oi),size(instFr,1)),ch));
                    phasedmed(oi,ti,ch) = temp(floor(numel(temp)/2)+1);
                end
            end
        end
        instFr = fs .* squeeze(median(phasedmed(:,sMirr+1:end-sMirr-1,:),1)) ./ (2*pi);
        dataIfre = vertcat(dataIfre, instFr);
    end
    clear phasedmed instFr sMirr datAll
    % Remove trials from current subject
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --------------- ATTENTION!!!
    eventprobs(:,1:length(y),:) = [];
    %clear eventprobs % Only for the case of ONE Subject
    % do phase Histogram
    [iFrePDF(:,:,:,sj), exIfe(:,:,sj), mi(:,:,sj)] = doInstFreHist(dataIfre,edges,y,probump,binVal,limFr);
    
%     % Difference in phase between channels
%     chPairs = combnk2(1:nCh,2);      % all combination of channels. Use my own 'combnk(2)' because original need statistical toolbox license
%     chDiffArg = zeros(size(dataIfre,1),size(chPairs,1));
%     for chi = 1:size(chPairs,1)
%         chDiffArg(:,chi) = dataIfre(:,chPairs(chi,1)) - dataIfre(:,chPairs(chi,2));
%     end
%     % do pahse Histograms of phase differences between pairs of channels
%     [argPDFpair, exZxPDFPair, miPair] = doPhaseHist(chDiffArg,edges,y,probump,binVal);
    %clear %chDiffArg
    
    % %%%%%%%%% -------- Randomization ------
    % Shuffle time points at random locations for each trial
    
    for r =1:nRep,
        dataRand = [];   
        for tr = 1:length(y),
            idr = randi([y(tr)-x(tr)-2]);
            dataRand = vertcat(dataRand, dataIfre(x(tr)+idr+1-tr:y(tr)-tr,:));
            dataRand = vertcat(dataRand, dataIfre(x(tr)+1-tr:x(tr)+idr-tr,:));
        end

        % do phase Histogram
        [~, ~, miRa(:,:,r,sj)] = doInstFreHist(dataRand,edges,y,probump,binVal,limFr);

        % Difference in phase between channels
%         chDiffArg = zeros(size(dataRand,1),size(chPairs,1));
%         for chi = 1:size(chPairs,1)
%             chDiffArg(:,chi) = dataRand(:,chPair(chi,1)) - dataRand(:,chPair(chi,1)) 
%         end
%         % do pahse Histograms of phase differences between pairs of channels
%         [~, ~, miPairRa(:,:,r)] = doPhaseHist(chDiffArg,edges,y,probump,binVal);
         clear dataRand
    end
    
    %%%%%%%% ---- Test that phase histograms are significatnt respect null
    %%%%%%%% ---- hypothesis (Phase histogram has uniform distribution).
    %meanMIra = mean(miRa,3);
    %stdMIra = std(miRa,[],3);
    %meanMIpaiRa = mean(miPairRa,3);
    %stdMIpaiRa = std(miPairRa,[],3);
    
    for bu = 1:nBump,
        for ch = 1:nCh,
            pValMi(ch,bu,sj) = length(find(miRa(ch,bu,:,sj) >= mi(ch,bu,sj))) / nRep;
        end
    end
    
    %%%%%%%%%%%%%% ---- visual representation -----
    if plotHist
        % Phase Histogram on polar form for all Channels
        side = 0.095;
        [~, chNames, Th, Rd] = readlocs(locElect);
        Th = pi/180*Th;                 % convert degrees to radians
        [yOr,xOr] = pol2cart(Th,Rd); % translate from polar to cart. coordinates
        yOr = yOr / (abs(max(yOr) - min(yOr)) + 3*side) + 0.5 + side/3;
        xOr = xOr / (abs(max(xOr) - min(xOr)) + 3*side) + 0.5 + side/3;
        position(:,1) = xOr - side;
        position(:,2) = yOr - side;
        position(:,3) = xOr + side;
        position(:,4) = yOr + side;
        % Test plot
        %plot(xOr,yOr,'o',position(:,1), position(:,2), '*',position(:,3), position(:,4), '+')
        colorbu = {'b','r','c','m','g'};
        fig = figure(sj);
        for i = 1:32
            subplot('Position', [position(i,1),position(i,2),side,side])
            set(0,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);
            for bu = 1:nBump
                h = plot(binVal,iFrePDF(i,:,bu));
                title(chNames(i),'HorizontalAlignment', 'left')
                h.LineWidth = 2.5;
                h.Color = colorbu{bu};
                xlim([limFr(1) limFr(2)])
                ylim([0 0.1])
                marksF = 1:4:nBins;
                set(gca,'Xtick',edges(marksF),'XTickLabel',edges(marksF))
                xlabel(['Frequency [Hz]'])
                line([4 4], [0 0.2], 'Color', 'k'); % 4 Hz lowedge
                line([9 9], [0 0.2], 'Color', 'k'); % 9 Hz highedge
                line([exIfe(ch,bu,sj) exIfe(ch,bu,sj)], [0.01 0.05], 'Color', colorbu{bu})
                hold on
            end
        end
        
%         %                     % save figure to a .fig
        namesave = ['sj' num2str(sj) freqBand 'InstFreModl' num2str(nBump) '.jpg' ];
        savepath = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/histPhase/';
        saveas(fig,[savepath namesave])
        clear fig
    end
    % test only with one subject
    %break
    
end

% Test differences on Instantaneous frequency between bumps
totSj = [1:nSj];
pairB = combnk2(1:nBump,2);     % all bump pairs
pValBuDif = zeros(nCh,length(pairB));
permT = zeros(nRep);
for ch = 1:nCh,
    for pr = 1:length(pairB),
        for re = 1:nRep,
            % permutation of a random number of labels
            sp = randi(nSj-1);
            permi = randperm(nSj,sp);       
            permNo = totSj(totSj~=permi);
            % do permuted variables
            var1 = [squeeze(exIfe(ch,pairB(pr,1),permi)) squeeze(exIfe(ch,pairB(pr,2),permNo))];
            var2 = [squeeze(exIfe(ch,pairB(pr,1),permNo)) squeeze(exIfe(ch,pairB(pr,2),permi))];
            % T-statistic
            [~,~,~,stats] = ttest2(var1, var2);
            permT(re) = stats.tstat;
        end
        % observed test statistic
        [~,~,~,stats] = ttest2(exIfe(ch,pairB(pr,1),:), exIfe(ch,pairB(pr,2),:));
        obsT = stats.tstat;
        % z-transform observed test statistic vs. null-hypothesis
        zval = (obsT - mean(permT)) / std(permT);
        % p-value
        pValBuDif(ch,pr) = 1 - normcdf(zval);
    end
end

save([savepath 'instFreqPDF' num2str(nBump) 'Bu.mat'], ...
   'exIfe','iFrePDF','pValBuDif','pValMi','','')

