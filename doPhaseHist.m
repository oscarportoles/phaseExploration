function [phasePDF, exPhaseZx, mi] = doPhaseHist(data,edges,y,probump,binVal)
    % assignt imaginary data of analytic signals to pahse bins scaled by
    % the bump location probability distribution. It computes Modulation
    % Index (mi) (Tort et.al. 2010 Measuring Phase-Amplitud coupling betwen..)
    % MI is the distance of the computed phase histogram to a uniformly 
    % distributed histogram
    % @ data: eeg data [samples, channels]
    % @ edges: pahse bin edges
    % @ y: end of trials
    % @ probump: bumps probability distribution [samples, bumps]
    % @ binVal: phase value of each bin [nBins]
    nBump = size(probump,2);
    nBins = length(binVal);
    nCh = size(data, 2);
    % assign to each sample a phase bin    
    BinIdx = zeros(size(data));       % [All samples x channels]
    for ch = 1:nCh
        BinIdx(:,ch) = discretize(data(:,ch), edges);
    end
    phasePDF = zeros(nCh,nBins,nBump);
    for bu = 1:nBump
        for ch = 1:nCh
            % sums the prob of a bump at each phase bin for all saples
            phasePDF(ch,:,bu) = accumarray(BinIdx(:,ch),probump(:,bu));
        end
    end       
    % Normalization
    phasePDF = phasePDF ./ length(y);
    % Expected phase from the phase PDF
    exPhaseZx = complex(zeros(nCh,nBump));   % [channels, subjects, model]
    for bu =1:nBump
        exPhaseZx(:,bu) = exp(1i*binVal) * squeeze(phasePDF(:,:,bu))';
    end
    % Modulation Index
    mi = zeros(nCh,nBump);
    for bu = 1:nBump,
        for ch = 1:nCh,
            mi(ch,bu) = (log(nBins) + sum([phasePDF(ch,:,bu) .* log(phasePDF(ch,:,bu) )])) ./ log(nBins);
        end
    end
end