% estimates the most cosistent connectivity across all samples of a
% network. It is donde after baseline correction and removal of spurious
% connections

clear all
% close all
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
namesRa = dir([pathdata '*RandTheta.mat']);
nSj = length(names);

load([pathdata 'electrodes.mat'])
% 
% ispcTraS = zeros(size(pairs,1),nBump,nSj);
% ispcStgS = zeros(size(pairs,1),nBump+1,nSj);
dwpliTraS = zeros(size(pairs,1),nBump,nSj);
dwpliStgS = zeros(size(pairs,1),nBump+1,nSj);

for sj = 14:14
    load([pathdata names(sj).name])
    load([pathdata namesRa(sj).name])
    % mean connectivity across network window
%     ispcStg     = squeeze(mean(ispcStg,1));
%     ispcTra     = squeeze(mean(ispcTra,1));
    dwpliTraStd    = squeeze(std(dwpliTra));
    dwpliStgStd    = squeeze(std(dwpliStg));
    dwpliStg    = squeeze(median(dwpliStg,1));
    dwpliTra    = squeeze(median(dwpliTra,1));
    
    % remove baseline
%     ispcStg     = ispcStg - repmat(abs(ispcPre),nBump+1,1)';
dwpliStg2    = dwpliStg - repmat(dwpliPre,nBump+1,1)';
    dwpliStg    = dwpliStg - repmat(dwpliPre,nBump+1,1)' - dwpliStgStd;
%     ispcTra     = ispcTra - repmat(abs(ispcPre),nBump,1)';
    dwpliTra    = dwpliTra - repmat(dwpliPre,nBump,1)' - dwpliTraStd;    
    % spurious volumn conductance connections
%     false0      = ispcGV0stg < pVal;
%     falsePi     = ispcGVpistg < pVal;
%     ispcStg(false0) = 0;
%     ispcStg(falsePi) = 0;
%     false0 = ispcGV0tra < pVal;
%     falsePi = ispcGVpitra < pVal;
%     ispcTra(false0) = 0;
%     ispcTra(falsePi) = 0;
    % threshold of noise
%      ndwpli  = ones(size(dwpliPre))*mean(dwpliPreStd);
%      ndwpliT = repmat(ndwpli,nBump,1)';
%      ndwpliS = repmat(ndwpli,nBump+1,1)';
     ndwpliCIT = repmat(dwpliCIt,1, nBump,1);
     ndwpliCIS = repmat(dwpliCIs,1, nBump+1,1);
     
    noiseBL = dwpliStg <= ndwpliCIS;
    dwpliStg(noiseBL) = 0;
    noiseBL = dwpliTra <= ndwpliCIT;
    dwpliTra(noiseBL) = 0;
    
%    ndwpli  = mean(dwpliPreStd);
%     nispc   = mean(ispcPreStd/sqrt(30));
%     noiseBL = ispcStg <= nispc;
%     ispcStg(noiseBL) = 0;
%     noiseBL = ispcTra <= nispc;
%     ispcTra(noiseBL) = 0;

%     noiseBL = dwpliStg <= ndwpliS;
%     dwpliStg(noiseBL) = 0;
%     noiseBL = dwpliTra <= ndwpliT;
%     dwpliTra(noiseBL) = 0;
    
    % remove those below the median + mean/median absoloute deviation
%     medwpliT = repmat(median(dwpliTra)+mad(dwpliTra,1),size(dwpliTra,1),1);
%     medwpliS = repmat(median(dwpliStg)+mad(dwpliStg,1),size(dwpliStg,1),1);  
%     medwpliT = repmat(mean(dwpliTra)+mad(dwpliTra,0),size(dwpliTra,1),1);
%     medwpliS = repmat(mean(dwpliStg)+mad(dwpliStg,0),size(dwpliStg,1),1); 
%     thMed = dwpliStg <= medwpliS;
%     dwpliStg(thMed) = 0;
%     thMed = dwpliTra <= medwpliT;
%     dwpliTra(thMed) = 0;    
%     ispcTraS(:,:,sj) = ispcTra;
%     ispcStgS(:,:,sj) = ispcStg;
    dwpliTraS(:,:,sj) = dwpliTra;
    dwpliStgS(:,:,sj) = dwpliStg;

%     % figures
    ispcTraSP = cell(1,nBump);
    ispcStgSP = cell(1,nBump+1);
    dwpliTraSP = cell(1,nBump); 
    dwpliStgSP = cell(1,nBump+1); 
    for bu = 1:nBump+1,
%         ispcStgStd   = std(nonzeros(ispcStg(:,bu)));
        dwpliStgStd  = std(nonzeros(dwpliStg(:,bu)));
%         ispcStgSm    = median(nonzeros(ispcStg(:,bu)));
        dwpliStgSm   = median(nonzeros(dwpliStg(:,bu)));

        if bu ~= nBump+1,
            dwpliTraStd  = std(nonzeros(dwpliTra(:,bu)));
%             ispcTraStd   = std(nonzeros(ispcTraS(:,bu)));
%             ispcTraSm    = median(nonzeros(ispcTraS(:,bu)));
            dwpliTraSm   = median(nonzeros(dwpliTra(:,bu)));
        end

        for pr = 1:size(pairs,1)
%             if ispcStg(pr,bu) > 0%ispcStgSm% + 1*ispcStgStd,
%                 ispcStgSP{bu} = vertcat(ispcStgSP{bu}, [pairs(pr,1),pairs(pr,2)]);
%             end
            if dwpliStg(pr,bu) > 0%dwpliStgSm + dwpliStgStd,
                dwpliStgSP{bu} = vertcat(dwpliStgSP{bu}, [pairs(pr,1),pairs(pr,2), dwpliStg(pr,bu)]);
            end
            if bu ~= nBump+1,
%                 if ispcTra(pr,bu) > 0%ispcTraSm% + 1*ispcTraStd,
%                     ispcTraSP{bu} = vertcat(ispcTraSP{bu}, [pairs(pr,1),pairs(pr,2)]);
%                 end
                if dwpliTra(pr,bu) > 0%dwpliTraSm + dwpliTraStd,
                    dwpliTraSP{bu} = vertcat(dwpliTraSP{bu}, [pairs(pr,1),pairs(pr,2),dwpliTra(pr,bu)]);
                end  
            end
        end
        figure(sj*10)
%         subplot(2,6,bu)
%         ds.chanPairs = ispcStgSP{bu};
%         topoplot_connect(ds, locElect)
%         title(['ISPC Stag. ' num2str(bu)])
        subplot(1,6,bu)
        ds.chanPairs = dwpliStgSP{bu}(:,[1,2]);
        ds.connectStrength = dwpliStgSP{bu}(:,3);
        topoplot_connect(ds, locElect)
        title(['dwpli Stag. ' num2str(bu)])
        set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',10);
        if bu ~= nBump+1,
            figure(sj*10+1)
%             subplot(2,5,bu)
%             ds.chanPairs = ispcTraSP{bu};
%             topoplot_connect(ds, locElect)
            subplot(1,5,bu)
            ds.chanPairs = dwpliTraSP{bu}(:,[1,2]);
            ds.connectStrength = dwpliTraSP{bu}(:,3);
            topoplot_connect(ds, locElect)
            title(['dwpli trans. ' num2str(bu)])
            set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',10);
        end
    end
    
    disp(['Subject:' num2str(sj)])
%     disp([length(nonzeros(dwpliStg(:,1))) length(nonzeros(dwpliStg(:,2))) length(nonzeros(dwpliStg(:,3))) ...
%         length(nonzeros(dwpliStg(:,4))) length(nonzeros(dwpliStg(:,5))) length(nonzeros(dwpliStg(:,6)))])
%     disp([length(nonzeros(dwpliTraSP(:,1))) length(nonzeros(dwpliTra(:,2))) length(nonzeros(dwpliTra(:,3))) ...
%         length(nonzeros(dwpliTra(:,4))) length(nonzeros(dwpliTra(:,5)))])
    disp([length(dwpliTraSP{1}),length(dwpliTraSP{2}),length(dwpliTraSP{3}),length(dwpliTraSP{4}),length(dwpliTraSP{5})])
    disp([length(dwpliStgSP{1}),length(dwpliStgSP{2}),length(dwpliStgSP{3}),length(dwpliStgSP{4}),length(dwpliStgSP{5}),length(dwpliStgSP{6})])
end  

% ispcTraS    = mean(ispcTraS,3);
% ispcStgS    = mean(ispcStgS,3);
dwpliTraSi   = sum(logical(dwpliTraS),3);
dwpliStgSi   = sum(logical(dwpliStgS),3);
% 
% ispcTraSP = cell(1,nBump);
% ispcStgSP = cell(1,nBump+1);
dwpliTraSP = cell(1,nBump); 
dwpliStgSP = cell(1,nBump+1); 


for bu = 1:nBump+1,
%     ispcStgStd   = std(nonzeros(ispcStgS(:,bu)));
    dwpliStgStd  = std(nonzeros(dwpliStgS(:,bu)));
%     ispcStgSm    = median(nonzeros(ispcStgS(:,bu)));
    dwpliStgSm   = median(nonzeros(dwpliStgS(:,bu)));
    
    if bu ~= nBump+1,
        dwpliTraStd  = std(nonzeros(dwpliTraS(:,bu)));
%         ispcTraStd   = std(nonzeros(ispcTraS(:,bu)));
%         ispcTraSm    = median(nonzeros(ispcTraS(:,bu)));
        dwpliTraSm   = median(nonzeros(dwpliTraS(:,bu)));
    end
    
    for pr = 1:size(pairs,1)
%         if ispcStgS(pr,bu) > 0%ispcStgSm% + 1*ispcStgStd,
%             ispcStgSP{bu} = vertcat(ispcStgSP{bu}, [pairs(pr,1),pairs(pr,2)]);
%         end
        if dwpliStgSi(pr,bu) > 12%dwpliStgSm + 2*dwpliStgStd,
            dwpliStgSP{bu} = vertcat(dwpliStgSP{bu}, [pairs(pr,1),pairs(pr,2)]);
        end
        if bu ~= nBump+1,
%             if ispcTraS(pr,bu) > 0%ispcTraSm% + 1*ispcTraStd,
%                 ispcTraSP{bu} = vertcat(ispcTraSP{bu}, [pairs(pr,1),pairs(pr,2)]);
%             end
            if dwpliTraSi(pr,bu) > 12%dwpliTraSm + 2*dwpliTraStd,
                dwpliTraSP{bu} = vertcat(dwpliTraSP{bu}, [pairs(pr,1),pairs(pr,2)]);
            end  
        end
    end

    figure(1)
    subplot(1,6,bu)
%     ds.chanPairs = ispcStgSP{bu};
%     topoplot_connect(ds, locElect)
%     title(['ISPC Stag. ' num2str(bu)])
%     subplot(2,6,bu+nBump+1)
    ds.chanPairs = dwpliStgSP{bu};
    topoplot_connect(ds, locElect)
    title(['Stage ' num2str(bu)])
    set(0,'DefaultAxesFontSize',20,'DefaultTextFontSize',20);
    if bu ~= nBump+1,
        figure(2)
%         subplot(2,5,bu)
%         ds.chanPairs = ispcTraSP{bu};
%         topoplot_connect(ds, locElect)
%         title(['ISPC trans. ' num2str(bu)])
        subplot(1,5,bu)
        ds.chanPairs = dwpliTraSP{bu};
        topoplot_connect(ds, locElect)
        title(['Transition ' num2str(bu)])
        set(0,'DefaultAxesFontSize',20,'DefaultTextFontSize',20);
    end
end