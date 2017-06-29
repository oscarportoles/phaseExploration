% Connectivity during stages. connectivity is evaluated on time for the
% duration of an stage and accros trials. Stages are defined as the window
% 'aPnts' around the middle point between the expected locataion of
% consecutive bumps

clear all
close all
clc

nBump = 5;                          % number of bumps in HSMM
nCh = 32;                           % number of channels
freqBand = 'theta';                 % frequency band to be analyzed
fs = 100;                           % sampling frequency [Hz]
pairCh = combnk2([1:nCh],2);        % all pairs of electrodes combinations [nComb, 2]
aPnts = 2;                        % number of points around the middle points of a stage to compute connectivty

% limits of each stage [ini, end]
stagLim = zeros(nBump+1,2);
stagLim(:,1) = 1:2*aPnts:2*aPnts*nBump+1;
stagLim(:,2) = 2*aPnts:2*aPnts:2*aPnts*(nBump+1);

info = ['fs = 100, Band: theta. connectivity on the phase space. Connectivity is assesed in the' ... 
       ' middle point btween expected locations plus a range of values aroun it (see stagLim)'];

% prepare data paths to load
pathdata = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/';
%varElec = [pathdata 'electrodes.mat'];
namePDF = ['AllSj100HSMM' num2str(nBump) 'Bout235.mat'];
codeHilb = ['*Bads100.mat'];
namesHilb = dir([pathdata codeHilb]);
nSj = length(namesHilb);

% test that the subjects are sorted equaly so loaded datasets agree
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
ispc = {};
pli = {};
wpli = {};
dwpli = {};
ispcCri = {};
for sj = 1:nSj
    % load datasets
    data = load([pathdata namesHilb(sj).name], freqBand, 'x', 'y');
    x = data.x;
    y = data.y;
    % do pahse angle from analytic signal
    dataZx = data.(freqBand);
    clear data
    %difPha = zeros((2*aPnts+1)*(nBump+1),size(pairCh,1),length(y));     % phase diff [stages concatenated, nPairs, nTrials]
    %csd = zeros((2*aPnts+1)*(nBump+1),size(pairCh,1),length(y));        % cross-spectral density

%     difPha = zeros(nCh,nCh,(2*aPnts+1)*(nBump+1),length(y));     % phase diff [nCh, nCh, nSamples per stage, nTrials]
    csd = zeros(nCh,nCh,(2*aPnts+1)*(nBump+1),length(y));     % phase diff [nCh, nCh, nSamples per stage, nTrials]
    for tr = 1:length(y)
        yTr = y(tr) - x(tr) + 1;
        % define stage bondaries:
            % bump expected locations
        expLo = [1:yTr] * squeeze(eventprobs(1:yTr,tr,:));
        expLo = [1 expLo yTr];
            % mean point between two bumps --> stage
        %meanStag = movmean(expLo, 2);  % matlab 2016a
        meanStag = zeros(1,nBump+1);
        for i=1:nBump+1, meanStag(i) = round(mean(expLo(i:i+1))); end
            % Indexes of neigboring points around the mean location of each stage
        stagIdx = [];
        for bu=1:nBump+1, 
            stagIdx = vertcat(stagIdx, meanStag(bu)-aPnts:meanStag(bu)+aPnts); 
        end  
        trPha = angle(dataZx(x(tr):y(tr),:));
        trZx = dataZx(x(tr):y(tr),:);
        for pr = 1:size(pairCh,1)
            % difference in phase
            %difPha(:,pr,tr) = trPha(stagIdx,pairCh(pr,1)) - trPha(stagIdx,pairCh(pr,2));
            
%             difPha(pairCh(pr,1),pairCh(pr,2),:,tr) = trPha(stagIdx,pairCh(pr,1)) - trPha(stagIdx,pairCh(pr,2));
            % cross-spectral density
            %csd(:,pr,tr) = trZx(stagIdx,pairCh(pr,1)) .* conj(trZx(stagIdx,pairCh(pr,2)));
            csd(pairCh(pr,1),pairCh(pr,2),:,tr) = trZx(stagIdx,pairCh(pr,1)) .* conj(trZx(stagIdx,pairCh(pr,2)));
        end
    end
    clear dataZx
    % Remove trials from current subject
    eventprobs(:,1:length(y),:) = [];
    % Intersite Phase Clustering
%     ispcD{sj} = abs(mean(exp(1i*difPha),4));
    ispcS = abs(mean(exp(1i*angle(csd)),4));
    ispcArgS = angle(mean(exp(1i*angle(csd)),4));
    % Phase-lag Index
%     pliD{sj}  = abs(mean(sign(imag(exp(1i*difPha))),4));
    pliS = abs(mean(sign(imag(csd)),4));
    pliArgS = angle(mean(sign(imag(csd)),4));
    % weighted phase-lag index
    wpliS = abs(mean(abs(imag(csd)).*sign(imag(csd)) ,4) )./mean(abs(imag(csd)),4);
%     wpliArg{sj} = angle(mean(abs(imag(csd)).*sign(imag(csd)) ,4) )./mean(abs(imag(csd)),4);
    % debiased weighted phase-lag index
    imagsum      = sum(imag(csd),4);
    imagsumW     = sum(abs(imag(csd)),4);
    debiasfactor = sum(imag(csd).^2,4);
    dwpliS  = (imagsum.^2 - debiasfactor)./(imagsumW.^2 - debiasfactor); 
    % Do symetric matrices with zeros in the diagonal
    for ni = 1:size(ispcS,3)  % iterate all samples at which synchrony is computed at stages
        % ISPC (it could be a function)
        symm = triu(squeeze(ispcS(:,:,ni)))' + triu(squeeze(ispcS(:,:,ni)));
        symm(logical(eye(size(symm)))) = 0;
        ispcS(:,:,ni) = symm;
        % PLI
        symm = triu(squeeze(pliS(:,:,ni)))' + triu(squeeze(pliS(:,:,ni)));
        symm(logical(eye(size(symm)))) = 0;
        pliS(:,:,ni) = symm;
        % wPLI
        symm = triu(squeeze(wpliS(:,:,ni)))' + triu(squeeze(wpliS(:,:,ni)));
        symm(logical(eye(size(symm)))) = 0;
        wpliS(:,:,ni) = symm;
        % dwpli
        symm = triu(squeeze(dwpliS(:,:,ni)))' + triu(squeeze(dwpliS(:,:,ni)));
        symm(logical(eye(size(symm)))) = 0;
        dwpliS(:,:,ni) = symm;
    end
    ispc{sj} = ispcS;
    pli{sj} = pliS;
    wpli{sj} = wpliS;
    dwpli{sj} = dwpliS;
    % critical ISPC value
    p = 0.01;
    ispcCri{sj} = sqrt(-log(p)/length(y));
    % gv-test, angle significance from 0 radians -> volumen conductance
    for ch1=1:nCh
        for ch2=ch1:nCh
            for ni = 1:size(ispcS,3)
                % ISPC
                expArg = exp((-ispcArgS(ch1,ch2,ni).^2) ./ (4*pi/length(y)));
                gvIspc(ch1,ch2,ni) = length(y) .* ispcS(ch1,ch2,ni) .* expArg .* (sqrt(2.0/length(y)));
                % PLI
                expArg = exp((-pliArgS(ch1,ch2,ni).^2) ./ (4*pi/length(y)));
                gvPli(ch1,ch2,:) = length(y) .* pliS(ch1,ch2,ni) .* expArg .* (sqrt(2.0/length(y))); 
            end
        end
    end
    gvTestIspc{sj} = gvIspc;
    gvTestPli{sj} = gvPli;
    clear gvIspc gvPli debiasfactor imagsumW imagsum expArg csd pliS wpliS dwpliS ispcS
end % End of subject

save([pathdata 'ISPC&PLIsatgSeg.mat'],'pli','ispc','dwpli', 'wpli',...
    'gvTestIspc','gvTestPli', 'info', 'stagLim', 'ispcCri')

        
    