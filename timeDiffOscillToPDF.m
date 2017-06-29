% time delay of transition pdf to the mean angle of all subjects at the
% transitions. You need to lead mauanaly the variable 'bpcAllArg' from
% wPhaseClust5Bu\freqband\.mat

nCh = 32;
nBump = 5;
% dT = 154;   % mean duration of a theta cycle in ms 
% dT = 87;   % mean duration of a alpha cycle in ms 
dT = 45;   % mean duration of a beta cycle in ms 
tDif = zeros(nCh,nBump);
plLim = [-11, 11];

bpcArg = angle(mean(exp(1i*bpcAllArg),3));

for ch =1:nCh,
    for bu = 1:nBump,
        if bpcArg(ch,bu) >= -pi/2 & bpcArg(ch,bu) < pi/2 % positive amplitude (Re)
            tDif(ch,bu) = bpcArg(ch,bu) * dT / (2*pi); 
        else
            if bpcArg(ch,bu) > 0,
                tDif(ch,bu) = - ( pi - bpcArg(ch,bu)) * dT / (2*pi); 
            elseif bpcArg(ch,bu) < 0,
                tDif(ch,bu) = ( pi - abs(bpcArg(ch,bu))) * dT / (2*pi);
            end
        end
    end
end

figure(3)
im = imagesc(tDif,plLim);
c = colorbar;
cc = colormap(jet);
xticks = {'1st Trans.','2nd Trans.','3rd Trans.','4th Trans.', '5th Trans.'};
title('\beta-band; time difference from transition PD and peak of oscillation [ms]')
set(im,'DefaultAxesFontSize',40,'DefaultTextFontSize',40);
set(gca,'Ytick',[1:nCh],'YTickLabel',{locElect(:).labels})
%set(gca,'Xtick',[1:nBump],'XTickLabel',xticks)
set(gca,'Xtick',[1:nBump],'XTickLabel',xticks,'XTickLabelRotation',45)
ylabel('Channels')
%xlabel(['Critical Phase Consistency: ' num2str(bpcCri)])
% xlabel('0.6 or higher consistency --> above-chance leve')
c.Label.String = ' Negative: Oscillation behind // Positive: Oscillation leads';
c.Limits = plLim;