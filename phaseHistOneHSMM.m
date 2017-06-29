% makes a histogram of phases with the bump PDF and the pahse of the
% analytic signal at a frequency band. It uses the results from a HSMM
% computed from the PCs of all subjects concatenated

clear all
close all
clc

nBump = 5;                          % number of bumps in HSMM
plotHist = 0;                       % plot head histogram per subject 
nBins = 18;                         % number of bins in the phase histogram
edges = linspace(-pi, pi, nBins+1);  % Does nBins from -pi to pi Bin edges
width = abs(edges(1) - edges(2));        % width of bins
binVal = edges(1:end-1) - diff(edges)/2;    % value of each bin
nCh = 32;                           % number of channels
freqBand = 'theta';                 % frequency band to be analyzed
nRep = 400;                         % number of repetition to build a random PDF

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

info = ['Phase histograms, Modulation Index, Bump Prefered pahse of a HSMM' ...
        ' calculated for all subject concatenated.Fs:100Hz,Band: ' freqBand];

pValMi = zeros(nCh,nBump,nSj);          % p-values of Modulation index
miRa = zeros(nCh,nBump,nRep,nSj);       % MI of randomized trials
argPDF = zeros(nCh,nBins,nBump,nSj);    % Phase probability distribution
exZxPDF = zeros(nCh,nBump,nSj);         % expected complex value of pahse PD
mi = zeros(nCh,nBump,nSj);              % MI of emprical data

for sj = 1:nSj
    % load datasets
    %varHilb =[pathdata namesHilb(sj).name];
    data = load([pathdata namesHilb(sj).name], freqBand, 'x', 'y');
    x = data.x;
    y = data.y;
    % do pahse angle from analytic signal
    datAll = angle(data.(freqBand)); % argment of analytic signal
    clear data
    probump = [];
    % Concatenate trials vertically to remove zeros out of trial
    for tr = 1:length(y)
        yTr = y(tr) - x(tr) + 1;
        probump = vertcat(probump, squeeze(eventprobs(1:yTr,tr,:)));
    end
    % Remove trials from current subject
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --------------- ATTENTION!!!
    eventprobs(:,1:length(y),:) = [];
    %clear eventprobs % Only for the case of ONE Subject
    % do phase Histogram
    [argPDF(:,:,:,sj), exZxPDF(:,:,sj), mi(:,:,sj)] = doPhaseHist(datAll,edges,y,probump,binVal);
    
%     % Difference in phase between channels
%     chPairs = combnk2(1:nCh,2);      % all combination of channels. Use my own 'combnk(2)' because original need statistical toolbox license
%     chDiffArg = zeros(size(datAll,1),size(chPairs,1));
%     for chi = 1:size(chPairs,1)
%         chDiffArg(:,chi) = datAll(:,chPairs(chi,1)) - datAll(:,chPairs(chi,2));
%     end
%     % do pahse Histograms of phase differences between pairs of channels
%     [argPDFpair, exZxPDFPair, miPair] = doPhaseHist(chDiffArg,edges,y,probump,binVal);
    %clear %chDiffArg
    
    % %%%%%%%%% -------- Randomization ------
    % Shuffle time points at random locations for each trial
    
    for r =1:nRep,
        dataRand = [];   
        for tr = 1:length(y),
            idr = randi([y(tr)-x(tr)]);
            dataRand = vertcat(dataRand, datAll(x(tr)+idr:y(tr),:));
            dataRand = vertcat(dataRand, datAll(x(tr):x(tr)+idr-1,:));
        end

        % do phase Histogram
        [~, ~, miRa(:,:,r,sj)] = doPhaseHist(dataRand,edges,y,probump,binVal);

        % Difference in phase between channels
%         chDiffArg = zeros(size(dataRand,1),size(chPairs,1));
%         for chi = 1:size(chPairs,1)
%             chDiffArg(:,chi) = dataRand(:,chPair(chi,1)) - dataRand(:,chPair(chi,1)) 
%         end
%         % do pahse Histograms of phase differences between pairs of channels
%         [~, ~, miPairRa(:,:,r)] = doPhaseHist(chDiffArg,edges,y,probump,binVal);
%         clear chDiffArg
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
    
%     %%%%%%%%%%%%%% ---- visual representation -----
    if plotHist
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
        fig = figure(sj);
        for i = 1:32
            subplot('Position', [position(i,1),position(i,2),side,side])
            set(0,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);
            %P = polaRad(edges, 0.15 * ones(size(edges)));
            %set(P, 'Visible', 'off')
            %hold on
            for bu = 1:nBump
                %plot([1,2])
                h = polaRad(edges,[argPDF(i,:,bu),argPDF(i,1,bu)]);
                title(chNames(i),'HorizontalAlignment', 'left')
                h.LineWidth = 2.5;
                h.Color = colorbu{bu};
                hold on
            end
        end
        
        %                     % save figure to a .fig
        namesave = ['sj' num2str(sj) freqBand 'Modl' num2str(nBump) '.jpg' ];
        savepath = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/histPhase/';
        saveas(fig,[savepath namesave])
        clear fig
    end
end

bpcMag = abs(mean(exp(1i.*angle(exZxPDF)),3));       % Bump Phase Clustering Magnitude
bpcMagZ = nSj .* bpcMag.^2;                              % transformation for few subjects
p = 0.05/nCh;
bpcCri = sqrt(-log(p)/nSj);
pValBpcZ = exp(sqrt(1 + 4*nSj + 4*(nSj^2 - (nSj*bpcMagZ).^2))-(1+2*nSj));
pValBpc = exp(sqrt(1 + 4*nSj + 4*(nSj^2 - (nSj*bpcMag).^2))-(1+2*nSj));

bpcPhs = angle(mean(exp(1i*angle(exZxPDF)),3));       % Bump Phase Clustering Phase
% Test if the diffrences between bump angles are significant
buPair = combnk2(1:nBump,2);
gvTest = zeros(nCh, size(buPair,1));
for pr = 1:size(buPair,1)       % Iterates backwards for clarity of resluts from 'buPair'
    diffPha = (bpcPhs(:,buPair(pr,1)) - bpcPhs(:,buPair(pr,2)));
    gvTest(:,pr) = nSj .* bpcPhs(:,buPair(pr,1)) .* exp((-(diffPha).^2) ./ (4*pi/nSj)) .*(sqrt(2./nSj));
end
pValBPCpha = 1 - normcdf(gvTest);
savepath = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/histPhase/';
nameSave = [savepath 'phasePDFmi' num2str(nBump) 'Bu36bins.mat'];
save(nameSave,'pValMi','miRa','argPDF','exZxPDF','mi','bpcMag','bpcPhs','bpcCri','pValBpcZ','pValBpc','gvTest','bpcPhs','pValBPCpha')

