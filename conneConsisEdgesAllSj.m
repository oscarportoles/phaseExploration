% estimates the most cosistent connectivity across all samples of a
% network. It is donde after baseline correction and removal of spurious
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
names = dir([pathdata '*Theta.mat']);
nSj = length(names);

load([pathdata 'electrodes.mat'])

ispcTraS = zeros(size(pairs,1),nBump,nSj);
ispcStgS = zeros(size(pairs,1),nBump+1,nSj);
dwpliTraS = zeros(size(pairs,1),nBump,nSj);
dwpliStgS = zeros(size(pairs,1),nBump+1,nSj);

for sj = 1:nSj
    load([pathdata names(sj).name])
    
    ispcTraCu = zeros(size(pairs,1),nBump);
    ispcTraCw = zeros(size(pairs,1),nBump);
    ispcStgCu = zeros(size(pairs,1),nBump+1);
    ispcStgCw = zeros(size(pairs,1),nBump+1);
    dwpliTraCu = zeros(size(pairs,1),nBump);
    dwpliTraCw = zeros(size(pairs,1),nBump);
    dwpliStgCu = zeros(size(pairs,1),nBump+1);
    dwpliStgCw = zeros(size(pairs,1),nBump+1);
    
%     ispcTraP = cell(1,nBump);
%     ispcStgP = cell(1,nBump+1);
%     dwpliTraP = cell(1,nBump); 
%     dwpliStgP = cell(1,nBump+1); 
    for bu = 1:nBump+1, 
        for pnt = 1:nPnt,
            ndwpli = mean(dwpliPreStd); %mean(abs(dwpliCI));
            nispc = mean(ispcPreStd/sqrt(30));
            % ISPC - baseline correction
            ispcStg(pnt,:,bu) = (squeeze(ispcStg(pnt,:,bu)) - ispcPre);% ./ ispcPre;
            % spurious volumn conductance connections
            false0 = ispcGV0stg(pnt,:,bu) < pVal;
            falsePi = ispcGVpistg(pnt,:,bu) < pVal;
            ispcStg(pnt,false0,bu) = 0;
            ispcStg(pnt,falsePi,bu) = 0;
            % connections below threshold
            noiseC = ispcStg(pnt,:,bu) <= nispc; %mean(ispcCI);
            ispcStg(pnt,noiseC,bu) = 0;
%             ispcTh = mean(nonzeros(ispcStg(pnt,:,bu))) + thispc * std(nonzeros(ispcStg(pnt,:,bu)));
%             noiseC = find(ispcStg(pnt,:,bu) < ispcTh);
%             ispcStg(pnt,noiseC,bu) = 0;
            %%% dwPLI - baseline correction
            dwpliStg(pnt,:,bu) = (dwpliStg(pnt,:,bu) - dwpliPre);% ./ dwpliPre;
            % connections below threshold
            noiseC = dwpliStg(pnt,:,bu) <=  ndwpli;% mean(dwpliCI);
            dwpliStg(pnt,noiseC,bu) = 0;
%             dwpliTh = mean(nonzeros(dwpliStg(pnt,:,bu))) + thpli * std(nonzeros(dwpliStg(pnt,:,bu)));
%             noiseC = find(dwpliStg(pnt,:,bu) < dwpliTh);
%             dwpliStg(pnt,noiseC,bu) = 0;
            if bu ~= nBump+1
                % ISPC - baseline correction
                ispcTra(pnt,:,bu) = (ispcTra(pnt,:,bu) - ispcPre);% ./ ispcPre;
                % spurious volumn conductance connections
                false0 = ispcGV0stg(pnt,:,bu) < pVal;
                falsePi = ispcGVpistg(pnt,:,bu) < pVal;
                ispcTra(pnt,false0,bu) = 0;
                ispcTra(pnt,falsePi,bu) = 0;
                % connections below threshold
                noiseC = ispcTra(pnt,:,bu) <= nispc; % mean(ispcCI);
                ispcTra(pnt,noiseC,bu) = 0;
%                 ispcTh = mean(nonzeros(ispcTra(pnt,:,bu))) + thispc * std(nonzeros(ispcTra(pnt,:,bu)));
%                 noiseC = find(ispcTra(pnt,:,bu) < ispcTh);
%                 ispcTra(pnt,noiseC,bu) = 0;
                %%% dwPLI - baseline correction
                dwpliTra(pnt,:,bu) = (dwpliTra(pnt,:,bu) - dwpliPre);% ./ dwpliPre;
                % connections below threshold
                noiseC = dwpliTra(pnt,:,bu) <= ndwpli; % mean(dwpliCI);
                dwpliTra(pnt,noiseC,bu) = 0;
%                 dwpliTh = mean(nonzeros(dwpliTra(pnt,:,bu))) + thpli * std(nonzeros(dwpliTra(pnt,:,bu)));
%                 noiseC = find(dwpliTra(pnt,:,bu) < dwpliTh);
%                 dwpliTra(pnt,noiseC,bu) = 0;
            end
        end
        % get consistent links
        for pr = 1:size(pairs,1)
            % ISPC Stages consistency of an edge
%             w = (1 - (wpdStg(:,bu,pr) / sum(wpdStg(:,bu,pr)))) / sum(1 - (wpdStg(:,bu,pr) / sum(wpdStg(:,bu,pr))));
%             edgeC = sum(w .* logical(ispcStg(:,pr,bu)));
            edgeC = mean(ispcStg(:,pr,bu));
            
            if edgeC > nispc*5%thCis
                ispcStgCu(pr,bu) = 1;
                ispcStgCw(pr,bu) = edgeC;
%                 ispcStgP{bu} = vertcat(ispcStgP{bu}, [pairs(pr,1),pairs(pr,2)]);
            else
                ispcStgCu(pr,bu) = 0;
                ispcStgCw(pr,bu) = 0;
            end
            % dwPLI Stages consistency of an edge
%             w = (1 - (wpdStg(:,bu,pr) / sum(wpdStg(:,bu,pr)))) / sum(1 - (wpdStg(:,bu,pr) / sum(wpdStg(:,bu,pr))));
%             edgeC = sum(w .* logical(dwpliStg(:,pr,bu)));
            edgeC = mean(dwpliStg(:,pr,bu));
            if edgeC > ndwpli*5%thCpl
                dwpliStgCu(pr,bu) = 1;
                dwpliStgCw(pr,bu) = edgeC;
%                 dwpliStgP{bu} = vertcat(dwpliStgP{bu}, [pairs(pr,1),pairs(pr,2)]);
            else
                dwpliStgCu(pr,bu) = 0;
                dwpliStgCw(pr,bu) = 0;
            end
            if bu ~= nBump+1
                % ISPC transitions consistency of an edge
%                 edgeC = sum((wpdTra(:,bu,pr) / sum(wpdTra(:,bu,pr))) .* logical(ispcTra(:,pr,bu)));
                edgeC = mean((ispcTra(:,pr,bu)));
                if edgeC > nispc*5%thCis
                    ispcTraCu(pr,bu) = 1;
                    ispcTraCw(pr,bu) = edgeC;
%                     ispcTraP{bu} = vertcat(ispcTraP{bu}, [pairs(pr,1),pairs(pr,2)]);
                else
                    ispcTraCu(pr,bu) = 0;
                    ispcTraCw(pr,bu) = 0;
                end
                % dwPLI transitions consistency of an edge
%                 edgeC = sum((wpdTra(:,bu,pr) / sum(wpdTra(:,bu,pr))) .* logical(dwpliTra(:,pr,bu)));
                edgeC = mean((dwpliTra(:,pr,bu)));
                if edgeC > ndwpli*5%thCpl
                    dwpliTraCu(pr,bu) = 1;
                    dwpliTraCw(pr,bu) = edgeC;
%                     dwpliTraP{bu} = vertcat(dwpliTraP{bu}, [pairs(pr,1),pairs(pr,2)]);
                else
                    dwpliTraCu(pr,bu) = 0;
                    dwpliTraCw(pr,bu) = 0;
                end
            end
        end
    end
    ispcTraS(:,:,sj) = ispcTraCw;
    ispcStgS(:,:,sj) = ispcStgCw;
    dwpliTraS(:,:,sj) = dwpliTraCw;
    dwpliStgS(:,:,sj) = dwpliStgCw;
end  

ispcTraS    = mean(ispcTraS,3);
ispcStgS    = mean(ispcStgS,3);
dwpliTraS   = mean(dwpliTraS,3);
dwpliStgS   = mean(dwpliStgS,3);

ispcTraSP = cell(1,nBump);
ispcStgSP = cell(1,nBump+1);
dwpliTraSP = cell(1,nBump); 
dwpliStgSP = cell(1,nBump+1); 


for bu = 1:nBump+1,
    ispcStgStd   = std(nonzeros(ispcStgS(:,bu)));
    dwpliStgStd  = std(nonzeros(dwpliStgS(:,bu)));
    ispcStgSm    = median(nonzeros(ispcStgS(:,bu)));
    dwpliStgSm   = median(nonzeros(dwpliStgS(:,bu)));
    
    if bu ~= nBump+1,
        dwpliTraStd  = std(nonzeros(dwpliTraS(:,bu)));
        ispcTraStd   = std(nonzeros(ispcTraS(:,bu)));
        ispcTraSm    = median(nonzeros(ispcTraS(:,bu)));
        dwpliTraSm   = median(nonzeros(dwpliTraS(:,bu)));
    end
    
    for pr = 1:size(pairs,1)
        if ispcStgS(pr,bu) > 0%ispcStgSm% + 1*ispcStgStd,
            ispcStgSP{bu} = vertcat(ispcStgSP{bu}, [pairs(pr,1),pairs(pr,2)]);
        end
        if dwpliStgS(pr,bu) > 0%dwpliStgSm + 2*dwpliStgStd,
            dwpliStgSP{bu} = vertcat(dwpliStgSP{bu}, [pairs(pr,1),pairs(pr,2)]);
        end
        if bu ~= nBump+1,
            if ispcTraS(pr,bu) > 0%ispcTraSm% + 1*ispcTraStd,
                ispcTraSP{bu} = vertcat(ispcTraSP{bu}, [pairs(pr,1),pairs(pr,2)]);
            end
            if dwpliTraS(pr,bu) > 0%dwpliTraSm + 2*dwpliTraStd,
                dwpliTraSP{bu} = vertcat(dwpliTraSP{bu}, [pairs(pr,1),pairs(pr,2)]);
            end  
        end
    end

    figure(1)
    subplot(2,6,bu)
    ds.chanPairs = ispcStgSP{bu};
    topoplot_connect(ds, locElect)
    title(['ISPC Stag. ' num2str(bu)])
    subplot(2,6,bu+nBump+1)
    ds.chanPairs = dwpliStgSP{bu};
    topoplot_connect(ds, locElect)
    title(['dwpli Stag. ' num2str(bu)])
    if bu ~= nBump+1,
        figure(2)
        subplot(2,5,bu)
        ds.chanPairs = ispcTraSP{bu};
        topoplot_connect(ds, locElect)
        title(['ISPC trans. ' num2str(bu)])
        subplot(2,5,bu+nBump)
        ds.chanPairs = dwpliTraSP{bu};
        topoplot_connect(ds, locElect)
        title(['dwpli trans. ' num2str(bu)])
    end
end