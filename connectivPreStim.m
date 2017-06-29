% Load band pass filtered baseline data, do the hilbert transform per
% trial, computes the phase clustering accross trials on the baseline and
% the median pahse connectivity of the baseline.

clear all
close all
clc

blPnt = 30;                     % prestimulus has 3 samples
nCh = 32;                       % number of channels
pairCh = combnk2([1:nCh],2);    % all pairs of electrodes combinations [nComb, 2]

savename = 'ISPC&PLIbaseline.mat';
pathdata = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/';
pathsave = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/';
snames = dir([pathdata '*PreStimEpochBands100.mat']);
nSj = length(snames);
varLoad = 'thetaEEGbl';
info = {'Baseline connectivty as the mean connectivity accross trial of a' ...
        '300 ms pre-stimulis time.'};
    
ispcBL = zeros(nCh,nCh,nSj);
pliBL = zeros(nCh,nCh,nSj);
wpliBL = zeros(nCh,nCh,nSj);
dwpliBL = zeros(nCh,nCh,nSj);

for sj = 1:nSj
    fileload = [pathdata snames(sj).name];
    load(fileload, varLoad);     % band specific data   
    data = reshape(thetaEEGbl,blPnt,[],nCh);    % new shape [blPnt, number trials, nCh]
    clear thetaEEGlb 
    
    nTr = size(data, 2);                        % number of trials
    csd = zeros(nCh,nCh,blPnt,nTr);      % phase diff [nCh,nCh,bl,nTr]
    for ch = 1:nCh,
        % Analytic signal by hilber transform
        data(:,:,ch) = hilbert(squeeze(data(:,:,ch)));
    end
    for pr = 1:size(pairCh,1),
        % pahse difference between pairs of channels,
        csd(pairCh(pr,1),pairCh(pr,2),:,:) = data(:,:,pairCh(pr,1)) .* conj(data(:,:,pairCh(pr,2)));
    end
    clear data
    % Intersite Phase Clustering
    ispcS = squeeze(abs(mean(exp(1i*angle(csd)),4)));    % across trials
    % Phase-lag Index
    pliS = squeeze(abs(mean(sign(imag(csd)),4)));
    % weighted phase-lag index
    wpliS = squeeze(abs(mean(abs(imag(csd)).*sign(imag(csd)) ,4) )./mean(abs(imag(csd)),4));
    % debiased weighted phase-lag index
    imagsum      = sum(imag(csd),4);
    imagsumW     = sum(abs(imag(csd)),4);
    debiasfactor = sum(imag(csd).^2,4);
    dwpliS  = squeeze((imagsum.^2 - debiasfactor)./(imagsumW.^2 - debiasfactor));
    
    % median connectivity value accross the pre-stimulus
    ispcS = median(ispcS,3);
    pliS = median(pliS,3);
    wpliS = median(wpliS,3);
    dwpliS = median(dwpliS,3);
    % Do symetric matrices with zeros in the diagonal
    % ISPC (it could be a function)
    symm = triu(ispcS)' + triu(ispcS);
    symm(logical(eye(size(symm)))) = 0;
    ispcS = symm;
    % PLI
    symm = triu(pliS)' + triu(pliS);
    symm(logical(eye(size(symm)))) = 0;
    pliS = symm;
    % wPLI
    symm = triu(wpliS)' + triu(wpliS);
    symm(logical(eye(size(symm)))) = 0;
    wpliS = symm;
    % dwpli
    symm = triu(dwpliS)' + triu(dwpliS);
    symm(logical(eye(size(symm)))) = 0;
    dwpliS = symm;
        
    ispcBL(:,:,sj) = ispcS;
    pliBL(:,:,sj)= pliS;
    wpliBL(:,:,sj) = wpliS;
    dwpliBL(:,:,sj) = dwpliS;
end

save([pathsave savename],'pliBL','ispcBL','dwpliBL', 'wpliBL', 'info')