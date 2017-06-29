% Does the connectivity matrices in nCh by nCh format after they are treated
% estimates the most cosistent connectivity across all samples of a
% network. It is donde after baseline correction, substraction of estimation noise (std) and removal of spurious
% connections

clear all
close all
clc

nCh = 32;                           % number of channels
nBump = 5;                          % number of transitions
nStg = nBump + 1;                   % number of stages
pairs = combnk2([1:nCh],2);        % all pairs of electrodes combinations [nComb, 2]
nPnt = 5;                           % number of samples on each network time window
thpli = 2;                          % dwPLI threshold, number of std above the mean
thispc = 1;                          % ISPC threshold, number of std above the mean
pVal = 0.05;                        % p-value threshold for angle GV-test
thCpl = 0.0;                          % consistency threshold phase lag index
thCis = 0.0;                       % consistency threshold inter-site pahse clustering

pathdata = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/PhaseConectiv/';
names = dir([pathdata '*RoTheta.mat']);
nSj = length(names);

load([pathdata 'electrodes.mat'])

cAllsj = struct([]);

for sj = 1:nSj
    load([pathdata names(sj).name])
    % variavility on the connectivity window
    dwpliTraStd    = squeeze(std(dwpliTra));
    dwpliStgStd    = squeeze(std(dwpliStg));
    % median connectivity across network window
    dwpliStg    = squeeze(median(dwpliStg,1));
    dwpliTra    = squeeze(median(dwpliTra,1));
    
    % remove baseline & unceratanity on connectivity network
    dwpliStg    = dwpliStg - repmat(dwpliPre,nBump+1,1)' - dwpliStgStd;
    dwpliTra    = dwpliTra - repmat(dwpliPre,nBump,1)' - dwpliTraStd;    
    % connectivity below confidence Interval
    ndwpliCIT = repmat(2*dwpliCItr,1, nBump,1);
    ndwpliCIS = repmat(2*dwpliCItr,1, nBump+1,1);
     
    noiseBL = dwpliStg <= ndwpliCIS;
    dwpliStg(noiseBL) = 0;
    noiseBL = dwpliTra <= ndwpliCIT;
    dwpliTra(noiseBL) = 0;
  
    % remove those below the median + mean/median absoloute deviation
%     medwpliT = repmat(median(dwpliTra)+mad(dwpliTra,1),size(dwpliTra,1),1);
%     medwpliS = repmat(median(dwpliStg)+mad(dwpliStg,1),size(dwpliStg,1),1);  
%     medwpliT = repmat(mean(dwpliTra)+mad(dwpliTra,0),size(dwpliTra,1),1);
%     medwpliS = repmat(mean(dwpliStg)+mad(dwpliStg,0),size(dwpliStg,1),1); 
%     thMed = dwpliStg <= medwpliS;
%     dwpliStg(thMed) = 0;
%     thMed = dwpliTra <= medwpliT;
%     dwpliTra(thMed) = 0;   
    for bu = 1:nBump+1,
        for pr = 1:size(pairs,1)
            xS(pairs(pr,1),pairs(pr,2)) = dwpliStg(pr,bu);
            xS(pairs(pr,2),pairs(pr,1)) = dwpliStg(pr,bu);
            if bu ~= nBump+1,
                xT(pairs(pr,1),pairs(pr,2)) = dwpliTra(pr,bu);
                xT(pairs(pr,2),pairs(pr,1)) = dwpliTra(pr,bu);
            end
        end
        eval(['csj(' num2str(sj) ').s' num2str(bu) '=xS;']);
        if bu ~= nBump+1,
            eval(['csj(' num2str(sj) ').t' num2str(bu) '=xT;']);
        end
    end
    
%     cAllsj(sj) = [csj];
%     clear csj
end  
info = 'Connectivity matrix via dwPLI, Fs=100Hz, Theta band, threshold ~ 0 (CI), 2017 March 247, s:stage, t:transition, 20 subjects; electrodes identifier';
save([pathdata 'connectivityD2703.mat'],'csj', 'info', 'locElect')
