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
ThC = 0.66;                          % consistency threshold

pathdata = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/PhaseConectiv/';
names = dir([pathdata '*Theta.mat']);
nSj = length(names);

load([pathdata 'electrodes.mat'])

for sj = 4:4
    load([pathdata names(sj).name])
    ispcTraCu = zeros(size(pairs,1),nBump);
    ispcTraCw = zeros(size(pairs,1),nBump);
    ispcStgCu = zeros(size(pairs,1),nBump+1);
    ispcStgCw = zeros(size(pairs,1),nBump+1);
    dwpliTraCu = zeros(size(pairs,1),nBump);
    dwpliTraCw = zeros(size(pairs,1),nBump);
    dwpliStgCu = zeros(size(pairs,1),nBump+1);
    dwplistgCw = zeros(size(pairs,1),nBump+1);
    
    ispcTraP = cell(1,nBump);
    ispcStgP = cell(1,nBump+1);
    dwpliTraP = cell(1,nBump); 
    dwpliStgP = cell(1,nBump+1); 
    for bu = 1:nBump+1, 
        for pnt = 1:nPnt,
            % ISPC - baseline correction
            ispcStg(pnt,:,bu) = ispcStg(pnt,:,bu) - ispcPre;
            % spurious volumn conductance connections
            false0 = find(ispcGV0stg(pnt,:,bu) < pVal);
            falsePi = find(ispcGVpistg(pnt,:,bu) < pVal);
            ispcStg(pnt,false0,bu) = 0;
            ispcStg(pnt,falsePi,bu) = 0;
            % connections below threshold
            noiseC = find(ispcStg(pnt,:,bu) < mean(ispcCI));
            ispcStg(pnt,noiseC,bu) = 0;
            ispcTh = mean(nonzeros(ispcStg(pnt,:,bu))) + thispc * std(nonzeros(ispcStg(pnt,:,bu)));
            noiseC = find(ispcStg(pnt,:,bu) < ispcTh);
            ispcStg(pnt,noiseC,bu) = 0;
            %%% dwPLI - baseline correction
            dwpliStg(pnt,:,bu) = dwpliStg(pnt,:,bu) - dwpliPre;
            % connections below threshold
            noiseC = find(dwpliStg(pnt,:,bu) < mean(dwpliCI));
            dwpliStg(pnt,noiseC,bu) = 0;
            dwpliTh = mean(nonzeros(dwpliStg(pnt,:,bu))) + thpli * std(nonzeros(dwpliStg(pnt,:,bu)));
            noiseC = find(dwpliStg(pnt,:,bu) < dwpliTh);
            dwpliStg(pnt,noiseC,bu) = 0;
            if bu ~= nBump+1
                % ISPC - baseline correction
                ispcTra(pnt,:,bu) = ispcTra(pnt,:,bu) - ispcPre;
                % spurious volumn conductance connections
                false0 = find(ispcGV0stg(pnt,:,bu) < pVal);
                falsePi = find(ispcGVpistg(pnt,:,bu) < pVal);
                ispcTra(pnt,false0,bu) = 0;
                ispcTra(pnt,falsePi,bu) = 0;
                % connections below threshold
                noiseC = find(ispcTra(pnt,:,bu) < mean(ispcCI));
                ispcTra(pnt,noiseC,bu) = 0;
                ispcTh = mean(nonzeros(ispcTra(pnt,:,bu))) + thispc * std(nonzeros(ispcTra(pnt,:,bu)));
                noiseC = find(ispcTra(pnt,:,bu) < ispcTh);
                ispcTra(pnt,noiseC,bu) = 0;
                %%% dwPLI - baseline correction
                dwpliTra(pnt,:,bu) = dwpliTra(pnt,:,bu) - dwpliPre;
                % connections below threshold
                noiseC = find(dwpliTra(pnt,:,bu) < mean(dwpliCI));
                dwpliTra(pnt,noiseC,bu) = 0;
                dwpliTh = mean(nonzeros(dwpliTra(pnt,:,bu))) + thpli * std(nonzeros(dwpliTra(pnt,:,bu)));
                noiseC = find(dwpliTra(pnt,:,bu) < dwpliTh);
                dwpliTra(pnt,noiseC,bu) = 0;
            end
        end
        % get consistent links
        for pr = 1:size(pairs,1)
            % ISPC Stages consistency of an edge
            w = (1 - (wpdStg(:,bu,pr) / sum(wpdStg(:,bu,pr)))) / sum((1 - (wpdStg(:,bu,pr) / sum(wpdStg(:,bu,pr)))));
            edgeC = sum(w .* logical(ispcStg(:,pr,bu)));
            if edgeC > ThC
                ispcStgCu(pr,bu) = 1;
                ispcStgCw(pr,bu) = edgeC;
                ispcStgP{bu} = vertcat(ispcStgP{bu}, [pairs(pr,1),pairs(pr,2)]);
            else
                ispcStgCu(pr,bu) = 0;
                ispcStgCw(pr,bu) = 0;
            end
            % dwPLI Stages consistency of an edge
            w = (1 - (wpdStg(:,bu,pr) / sum(wpdStg(:,bu,pr)))) / sum((1 - (wpdStg(:,bu,pr) / sum(wpdStg(:,bu,pr)))));
            edgeC = sum(w .* logical(dwpliStg(:,pr,bu)));
            if edgeC > ThC
                dwpliStgCu(pr,bu) = 1;
                dwpliStgCw(pr,bu) = edgeC;
                dwpliStgP{bu} = vertcat(dwpliStgP{bu}, [pairs(pr,1),pairs(pr,2)]);
            else
                dwpliStgCu(pr,bu) = 0;
                dwpliStgCw(pr,bu) = 0;
            end
            if bu ~= nBump+1
                % ISPC transitions consistency of an edge
                edgeC = sum((wpdTra(:,bu,pr) / sum(wpdTra(:,bu,pr))) .* logical(ispcTra(:,pr,bu)));
                if edgeC > ThC
                    ispcTraCu(pr,bu) = 1;
                    ispcTraCw(pr,bu) = edgeC;
                    ispcTraP{bu} = vertcat(ispcTraP{bu}, [pairs(pr,1),pairs(pr,2)]);
                else
                    ispcTraCu(pr,bu) = 0;
                    ispcTraCw(pr,bu) = 0;
                end
                % dwPLI transitions consistency of an edge
                 edgeC = sum((wpdTra(:,bu,pr) / sum(wpdTra(:,bu,pr))) .* logical(dwpliTra(:,pr,bu)));
                if edgeC > ThC
                    dwpliTraCu(pr,bu) = 1;
                    dwpliTraCw(pr,bu) = edgeC;
                    dwpliTraP{bu} = vertcat(dwpliTraP{bu}, [pairs(pr,1),pairs(pr,2)]);
                else
                    dwpliTraCu(pr,bu) = 0;
                    dwpliTraCw(pr,bu) = 0;
                end
            end
        end
        if bu ~= nBump+1,
            figure(10*sj)
            subplot(2,5,bu)
            ds.chanPairs = ispcTraP{bu};
            topoplot_connect(ds, locElect)
            title(['ISPC trans. ' num2str(bu)])
            subplot(2,5,bu+nBump)
            ds.chanPairs = dwpliTraP{bu};
            topoplot_connect(ds, locElect)
            title(['dwpli trans. ' num2str(bu)])
        end
        figure(10*sj+1)
        subplot(2,6,bu)
        ds.chanPairs = ispcStgP{bu};
        topoplot_connect(ds, locElect)
        title(['ISPC Stag. ' num2str(bu)])
        subplot(2,6,bu+nBump+1)
        ds.chanPairs = dwpliStgP{bu};
        topoplot_connect(ds, locElect)
        title(['dwpli Stag. ' num2str(bu)])        
    end
end        

