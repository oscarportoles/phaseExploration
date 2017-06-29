% processing of connectivity indexes to plot them out

nCh = 32;
pair = combnk2([1:nCh],2);        % all pairs of electrodes combinations [nComb, 2]
% bu = 5;
% pnt = 3;
pVal = 0.05;

ct = 1;
figure(3)
for pnt = 1:5
    for bu = 1:5,
        ispc = ispcStg(pnt,:,bu) - ispcPre;
        % spurious connections
        false0 = find(ispcGV0stg(pnt,:,bu) < pVal);
        falsePi = find(ispcGVpistg(pnt,:,bu) < pVal);

        ispc(false0) = 0;
        ispc(falsePi) = 0;

        % connections below threshold
        noiseC = find(ispc < mean(ispcCI));
        ispc(noiseC) = 0;

        connec = [];
        ispcw = [];
        for i = 1:length(ispc),
            if ispc(i) > 0,
                connec = vertcat(connec,[pair(i,1),pair(i,2)]);
                ispcw = vertcat(ispcw,ispc(i));
            end
        end
        subplot(5,5,ct)
        ct = ct + 1;
        %     figure(1)
        ds = struct();
        ds.chanPairs = connec;
        ds.connectStrength = ispcw;
        topoplot_connect(ds, locElect)
        title(['Stgnsition: ', num2str(bu),' Point: ',num2str(pnt)])
    end
end
ct = 1;
figure(4)
for pnt = 1:5
    for bu = 1:5,
        % debiased weighted phase lag index
        dwpli = dwpliStg(pnt,:,bu) - dwpliPre;
        dwpliTh = median(dwpli) + 2*std(dwpli);

        % connections below threshold
        noiseC = find(dwpli < mean(dwpliCI));
        dwpli(noiseC) = 0;
        noiseC = find(dwpli < dwpliTh);
        dwpli(noiseC) = 0;

        connec = [];
        dwpliw = [];
        for i = 1:length(ispc),
            if dwpli(i) > 0,
                connec = vertcat(connec,[pair(i,1),pair(i,2)]);
                dwpliw = vertcat(dwpliw,dwpli(i));
            end
        end
        ds = struct();
        ds.chanPairs = connec;
        ds.connectStrength = dwpliw;
        subplot(5,5, ct)
        ct = ct+1;
        topoplot_connect(ds, locElect)
        title(['Stgnsition: ', num2str(bu),' Point: ',num2str(pnt)])

    end
end