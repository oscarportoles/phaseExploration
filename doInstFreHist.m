function [iFrePDF, exIfre, mi] = doInstFreHist(data,edges,y,probump,binVal,lim)
    % assignt instantaneous frequency of analytic signals to frequency bins scaled by
    % the bump location probability distribution. It computes Modulation
    % Index (mi) (Tort et.al. 2010 Measuring Phase-Amplitud coupling betwen..)
    % MI is the distance of the computed phase histogram to a uniformly 
    % distributed histogram
    % @ data: eeg data [samples -1 , channels] (-1 because 'diff()')
    % @ edges: frequency bin edges
    % @ y: end of trials
    % @ probump: bumps probability distribution [samples, bumps]
    % @ binVal: frequency value of each bin [nBins]
    % @ lim: edges of the histogram to be buit [min, max] Hz
    nBump = size(probump,2);
    nBins = length(binVal);
    nCh = size(data, 2);
    fill = 1e-50;     % fill up empty bins to avoid dimension mismatchs. If it is NaN the IM = NaN
    % assign outlier frequencies to the closes extreme of the frequenies
    data(data<=lim(1)) = lim(1) + 0.1;
    data(data>=lim(2)) = lim(2) - 0.1;
    % assign to each sample a phase bin    
    BinIdx = zeros(size(data));       % [All samples x channels]
    for ch = 1:nCh
        BinIdx(:,ch) = discretize(data(:,ch), edges);
    end
    iFrePDF = zeros(nCh,nBins,nBump);
    for bu = 1:nBump
        for ch = 1:nCh
            % sums the prob of a bump at each phase bin for all saples
            iFrePDF(ch,:,bu) = accumarray(BinIdx(:,ch),probump(:,bu),[nBins,1]);
        end
    end       
    % Normalization
    iFrePDF = iFrePDF ./ length(y);
    % Expected Instantaneous frequency from the Inst. freq. PDF
    exIfre = zeros(nCh,nBump);   % [channels, subjects, model]
    for bu =1:nBump
        exIfre(:,bu) = binVal * squeeze(iFrePDF(:,:,bu))';
    end
    % Modulation Index
    mi = zeros(nCh,nBump);
    for bu = 1:nBump,
        for ch = 1:nCh,
            mi(ch,bu) = (log(nBins) + sum([iFrePDF(ch,:,bu) .* log(iFrePDF(ch,:,bu) )])) ./ log(nBins);
        end
    end
end