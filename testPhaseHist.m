sj = 1;
binIdxF = zeros(size(dataH(find(subjectsF == sj)),1),nCh);       % [All samples x channels]

dataF = [];
for tr = 1:length(y20{sj})
    dataF = vertcat(dataF, dataH(x20{sj}(tr)+31:y20{sj}(tr),:));
    dataF = vertcat(dataF, dataH(x20{sj}(tr):x20{sj}(tr)+30,:));
end

for ch = 1:nCh
    binIdxF(:,ch) = discretize(dataF(find(subjectsF == sj),ch), edges);
end

phasePDFf = zeros(nCh,nBins,bumpModl);

for bu = 1:bumpModl
    for ch = 1:nCh
        % sums the prob of a bump at each phase bin for all saples
        phasePDFf(ch,:,bu) = accumarray(binIdxF(:,ch),probump(find(subjectsF == sj),bu));
    end
end
% normalization by number of trials
phasePDFf = phasePDFf ./ length(y20{sj});

freqBand = 'Theta';
fig = figure(1);
titlen = ['Phase distribution of: ' freqBand ', subject: ' num2str(sj), ' @sum'];
colorB = {'b', 'r','c','m','g','k','y','b'};
colorB(bumpModl+1:end) = [];
legendB = {'Bump 1','Bump 2','Bump 3','Bump 4','Bump 5','Bump 6','Bump 7','Bump 8'};
legendB(bumpModl+1:end) = [];
ymax = max(max(max(phasePDFf))) + 0.1 * max(max(max(phasePDFf)));
plottopo(phasePDFf, 'chanlocs', locElect, 'ylim', [0,ymax], 'title',titlen, 'colors', colorB, ...
                        'ydir', 1,'legend', legendB, 'showleg','on')