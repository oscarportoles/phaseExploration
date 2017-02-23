
tr = 13;

% Probability distributions of bumps
figure(tr)
yTr = y(tr) - x(tr) + 1;
plot(zeros(yTr,1),'LineWidth', 0.1)
hold on
plot([squeeze(eventprobs(1:yTr,tr,:))], 'LineWidth', 4)
title('Probability Distribution of Transition Locations')
ylabel('Probability')
xlabel('Samples [n] // One Trial //')
legend({'1st trans.','2nd trans.','3rd trans.','4th trans.','5th trans.'})
set(0,'DefaultAxesFontSize',40,'DefaultTextFontSize',40);


for bu = 1:bumpModl
    for tr = 1:length(y20{sj})
        ym(tr) = y20{sj}(tr) - x20{sj}(tr) + 1;
        sumProbtr(tr,bu) = sum(eventprob20i{sj}(:,tr,bu));
        nonZeroProb(tr,bu) = nnz(eventprob20i{sj}(:,tr,bu)); 
    end
end
%ym = ym';
figure(1)
plot(1:tr, sumProbtr(:,1))
figure(2)
plot(1:tr, sumProbtr(:,2))
figure(3)
plot(1:tr, sumProbtr(:,3))
figure(4)
plot(1:tr, sumProbtr(:,4))
figure(5)
plot(1:tr, sumProbtr(:,5))



for bu = 1:bumpModl
    for tr = 1:length(y20{sj})
        yq(tr) = y20{sj}(tr) - x20{sj}(tr) + 1;
        %sumProbtr(tr,bu) = sum(eventprobSj(:,tr,bu));
        nonZeroProbQ(tr,bu) = nnz(eventprobs5(:,tr,bu)); 
    end
end

title('Theta band phase and bump 1 PDF')
yyaxis left
plot(dataH(x(1):y(1),1))
ylabel('Phase')
ylim([-pi pi])
yyaxis right
plot(probump(x(1):y(1),1))
ylim([0.2 0])
ylabel('probability [0, 1]')
xlabel('Samples [n]')

for i =1:y(1)
    sampleBin(i) = binVal(binIdx(i));
end


% Real part of EEG trial
tr = 11;
figure(tr)
yTr = y(tr) - x(tr) + 1;
plot(squeeze(real(dataZx(x(tr):y(tr),tr,1))), 'LineWidth', 3)
title('Amplitude & Phase of EEG Theta Band')
ylabel('Amplitude')
xlabel('Samples [n] // One Trial //')
%legend({'1st transition','2nd transition','3rd transition','4th transition','5th transition'})
set(0,'DefaultAxesFontSize',40,'DefaultTextFontSize',40);

% Real part and pahse
tr = 11;
yTr = y(tr) - x(tr) + 1;
t = 1:yTr;
[hAx,hphase,hprob] =  plotyy(t, angle(dataZx(x(tr):y(tr),1)),t,real(dataZx(x(tr):y(tr),1)));
% hphase.LineStyleOrder = {'o'};
% hprob.LineStyleOrder = '*';
hphase.LineWidth = 4;
hprob.LineWidth = 4;
xlabel('Samples [n] // One Trial //')
ylabel(hAx(1),'Amplitud') % left y-axis
hAx(1).YTick = [-pi,-pi/2,0,pi/2,pi];
hAx(1).YTickLabel = {'\pi','\pi/2','0','-\pi/2','-\pi'};
hAx(2).FontSize = 40;
hAx(1).FontSize = 40;
ylabel(hAx(1),'Phase [rad]') % right y-axis
ylabel(hAx(2),'Amplitud theta band') % right y-axis
title('Amplitude & Phase of EEG Theta Band')
legend({'EEG Phase', 'EEG Ampli.'})
set(0,'DefaultAxesFontSize',40,'DefaultTextFontSize',40);



% Phase and probability distributions
tr = 11;
yTr = y(tr) - x(tr) + 1;
t = 1:yTr;
nBins = 18;
edges = linspace(-pi, pi, nBins+1);  % Does nBins from -pi to pi Bin edges
[hAx,hphase,hprob] =  plotyy(t, angle(dataZx(x(tr):y(tr))),t,[squeeze(eventprobs(1:yTr,tr,:))]);
% hphase.LineStyleOrder = {'o'};
% hprob.LineStyleOrder = '*';
hphase.LineWidth = 4;
hprob(1).LineWidth = 4;
hprob(2).LineWidth = 4;
hprob(3).LineWidth = 4;
hprob(4).LineWidth = 4;
hprob(5).LineWidth = 4;
xlabel('Samples [n] // One Trial //')
ylabel(hAx(1),'Phase [rad]') % left y-axis
hAx(1).YLim = 1.05 .* [-pi, pi];
%hAx(1).YTick = [-pi,-pi/2,0,pi/2,pi];
%hAx(1).YTickLabel = {'pi','pi/2','0','-pi/2','-pi'};
hAx(1).YTick = edges;
hAx(1).YTickLabel = {'-\pi','','','-2\pi/3','','','-\pi/3','','','0','','','\pi/3','','','2\pi/3','','','\pi'};
hAx(1).YGrid = 'on';
hAx(2).YLim = [0, 0.3];
hAx(2).YTick = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35];
hAx(2).YTickLabel = {'0','0.05','0.1','0.15','0.2','0.25',' '};
hAx(1).FontSize = 40;
hAx(2).FontSize = 40;
ylabel(hAx(2),'Probability [0, 1]') % right y-axis
title('EEG Theta Band Phase & Transition Probability Distributions')
legend({'EEG phase','1st trans.','2nd trans.','3rd trans.','4th trans.','5th trans.'})
set(0,'DefaultAxesFontSize',40,'DefaultTextFontSize',40);



%
figure(2);

set(0,'DefaultAxesFontSize',16,'DefaultTextFontSize',16);
h = polaRad(edges,[phasePDF(1,:,1),phasePDF(1,1,1)]);
title('Histogram of theta band phase at bump 1')
h.LineWidth = 2.5;
hold on
set(0,'DefaultAxesFontSize',16,'DefaultTextFontSize',16);
h = polaRad(edges,[phasePDF(1,:,2),phasePDF(1,1,2)]);
title('Histogram of theta band phase at bump 1')
h.LineWidth = 2.5;
h.Color = 'r';

% all head histograms
side = 0.095;
[~, chNames, Th, Rd] = readlocs(locElect);
Th = pi/180*Th;                 % convert degrees to radians
[yOr,xOr] = pol2cart(Th,Rd); % translate from polar to cart. coordinates
% yOr = yOr + 0.5;
% xOr = xOr + 0.5;
yOr = yOr / (abs(max(yOr) - min(yOr)) + 3*side) + 0.5 + side/3;
xOr = xOr / (abs(max(xOr) - min(xOr)) + 3*side) + 0.5 + side/3;
position(:,1) = xOr - side;
position(:,2) = yOr - side;
position(:,3) = xOr + side;
position(:,4) = yOr + side;
% Test plot
%plot(xOr,yOr,'o',position(:,1), position(:,2), '*',position(:,3), position(:,4), '+')
colorbu = {'b','r','c','m','g'};
%side2 = 0.03;
sj = 12;
figure(4)
for i = 1:32
    subplot('Position', [position(i,1),position(i,2),side,side])
    set(0,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);
    %P = polaRad(edges, 0.15 * ones(size(edges)));
    %set(P, 'Visible', 'off')
    %hold on
    for bu = 1:5
        %plot([1,2])
        h = polaRad(edges,[phasePDF(i,:,bu),phasePDF(i,1,bu)]);
        title(chNames(i),'HorizontalAlignment', 'left')
        h.LineWidth = 2.5;
        h.Color = colorbu{bu};
        hold on
    end
end


% Example of bump histogram, uniform distribution and random distribution
% Matlab default colors (from 2 to 6) after 2014 
colorbu = [ 0 0.4470 0.7410;
            0.8500 0.3250 0.0980;
            0.2500 0.2500 0.2500];
uniforRa = unifrnd(1./nBins-0.007, 1./nBins+0.007,1,nBins );
uniforRa = uniforRa ./ sum(uniforRa);
unifor = ones(1,nBins)./ nBins;
ch = 2;
figure(1)
set(0,'DefaultAxesFontSize',40,'DefaultTextFontSize',40);
h = polaRad(edges,[argPDF(ch,:,bu,sj),argPDF(ch,1,bu,sj)]);
h.LineWidth = 4;
h.Color = colorbu(2,:);
hold on
h = polaRad(edges,[uniforRa,uniforRa(1)]);
h.LineWidth = 4;
h.Color = colorbu(1,:);
h = polaRad(edges,[unifor,unifor(1)]);
h.LineWidth = 4;
h.Color = colorbu(3,:);
legend({'empirical', 'random', 'uniform'})
title('e.g. significance uniform distribution')


% One channel plots
% Matlab default colors (from 2 to 6) after 2014 
colorbu = [0.8500 0.3250 0.0980;
           0.9290 0.6940 0.1250;
           0.4940 0.1840 0.5560;
           0.4660 0.6740 0.1880;
           0.3010 0.7450 0.9330];

ch = 2;
figure(1)
set(0,'DefaultAxesFontSize',40,'DefaultTextFontSize',40);
for bu = 1:5
    %plot([1,2])
    h = polaRad(edges,[argPDF(ch,:,bu,sj),argPDF(ch,1,bu,sj)]);
    title('Phase polar histogram')
    h.LineWidth = 4;
    h.Color = colorbu(bu,:);
    hold on
end
legend({'1st trans.','2nd trans.','3rd trans.','4th trans.','5th trans.'},'Location', 'eastoutside')


% ISPC all head
[carX, carY] = pol2cart(bpcPhs, bpcMag);  % pol2cart([angle, radious])
figure(5)
for ch = 1:32
    subplot('Position', [position(ch,1),position(ch,2),side,side])
    set(0,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);
    for bu = 1:5
        %plot([1,2])
        h = compassRad(carX(ch,bu), carY(ch,bu));
        title(chNames(ch))
        h.LineWidth = 2.5;
        h.Color = colorbu{bu};
        hold on
    end
    P = polaRad(edges, bpcCri * ones(size(edges)));
    P.Color = 'k';
    hold on
end
legend({'1st bump', '2nd bump','3th bump','4th bump','5th bump'})

% One channel ISPC
[carX, carY] = pol2cart(ISPA,ISPC);
ch = 31;
figure(1)
set(0,'DefaultAxesFontSize',16,'DefaultTextFontSize',16);
for bu = 1:5
    %plot([1,2])
    h = compassRad(carX(ch,bu), carY(ch,bu));
    %title(chNames(ch))
    title('OZ: ISPC & ISPA')
    h.LineWidth = 2.5;
    h.Color = colorbu{bu};
    hold on
end
P = polaRad(edges, ISPCri * ones(size(edges)));
P.Color = 'k';

% plots needs data from phase histograms analysis (file: phaseHistOnseHSMM.m)

% Test significancy on the prefrec phase angle between bumps at subject
% level
figure(1)
xticks = {'4th vs 5th','3rd vs 5th','3rd vs 4th','2nd vs 5th','2nd vs 4th','2nd vs 3rd','1st vs 5th','1st vs 4th','1st vs 3th',' 1st vs 2nd'};
alpha = 0.01;
im = imagesc(pValBPCpha,[0,alpha]);
c = colorbar;
colormap(autumn)
title('p-values testing differences between bumps phase angle')
set(im,'DefaultAxesFontSize',14,'DefaultTextFontSize',14);
set(gca,'Ytick',[1:nCh],'YTickLabel',{locElect(:).labels})
set(gca,'Xtick',[1:size(buPair,1)],'XTickLabel',xticks,'XTickLabelRotation',45)
xlabel('Bump "X" vs Bump "Y"')
ylabel('Channels')
c.Label.String = '';

% Test Bump Phase Consistency BPC is higher than a Critical BPC. BPC higher
% thant Critical BPC indicates that the Phase ocnsistency is ibove chance
figure(22)
p = 0.05/(nCh*nBump);
bpcCriH = sqrt(-log(p)/nSj);
im = imagesc(bpcMag,[bpcCriH,1]);
c = colorbar;
cc = colormap(autumn);
cc(1,:) = [0,0,0];
colormap(cc)
xticks = {'1st Trans.','2nd Trans.','3rd Trans.','4th Trans.', '5th Trans.'};
title('Phase clustering across all subjects')
set(im,'DefaultAxesFontSize',40,'DefaultTextFontSize',40);
set(gca,'Ytick',[1:nCh],'YTickLabel',{locElect(:).labels})
%set(gca,'Xtick',[1:nBump],'XTickLabel',xticks)
set(gca,'Xtick',[1:nBump],'XTickLabel',xticks,'XTickLabelRotation',45)
ylabel('Channels')
%xlabel(['Critical Phase Consistency: ' num2str(bpcCri)])
% xlabel('0.6 or higher consistency --> above-chance leve')
c.Label.String = 'Clustering [0, 1]';
c.Ticks = [7:10]./10;
% c.TickLabels = {'0.6'};
% c.TickLength = 0.05;



% Plot consistency of phase angles between channels in a subjecta and a
% single bump
p = 0.05/(nSj*nBump);
bpcCriCh = sqrt(-log(p)/nCh);
figure(33)
im = imagesc(squeeze(abs(mean(exp(1i.*angle(exZxPDF)),1)))',[bpcCriCh,1]);
c = colorbar;
cc = colormap(autumn);
cc(1,:) = [0,0,0];
colormap(cc)
xticks = {'1st Trans.','2nd Trans.','3rd Trans.','4th Trans.', '5th Trans.'};
set(gca,'Xtick',[1:nBump],'XTickLabel',xticks,'XTickLabelRotation',45)
%set(gca,'CLimMode', 'manual', 'AmbientLightColor', 'r')
yticks = {'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12','S13','S14','S15','S16','S17','S18','S19','S20'};
set(gca,'Ytick',[1:nSj],'YTickLabel',yticks)
ylabel('Subjects')
title('Phase clustering across all channels')
set(im,'DefaultAxesFontSize',40,'DefaultTextFontSize',40);
c.Label.String = 'Clustering [0, 1]';
c.Ticks = [5:10]./10;


% Plot and do a table of the number of p-values avobe alpha per bump and
% channel
countPval = zeros(nCh,nBump);
alpha = 0.05/(nCh*nSj*nBump); % Bomferroni correction
for bu =1: nBump,
    for ch = 1:nCh,
        countPval(ch,bu) = length(find(pValMi(ch,bu,:) < alpha));
    end
end


figure(4)
im = imagesc(countPval);
c = colorbar;
colormap(gray(max(max(countPval))+1))
xticks = {'1st Trans.','2nd Trans.','3rd Trans.','4th Trans.', '5th Trans.'};
title('Significance analysis ( \alpha = 1.5x10^{-5})')
set(im,'DefaultAxesFontSize',40,'DefaultTextFontSize',40);
set(gca,'Ytick',[1:nCh],'YTickLabel',{locElect(:).labels})
%set(gca,'Xtick',[1:nBump],'XTickLabel',xticks)
set(gca,'Xtick',[1:nBump],'XTickLabel',xticks,'XTickLabelRotation',45)
ylabel('Channels')
c.Label.String = 'Number of subjects with significant Phase clustering';
% c.Ticks = [0:max(max(countPval))+1];
% c.TickLabels = {'0', '1','2','3','4','5','6','7'};\

% Golbal field variability. gloabal mean of VAriability accross channels
% All trials are time locked to the bumps expected location. The signals
% were contracted or expanded to be time locked
figure(1)
plot(mean(gfEnv,2))
xlabel('samples')
ylabel('Envelop Variability [std]')
title('Global mean of between-channel Envelop variability')
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
for i = 1:6, line([duration(i) duration(i)], [1.7 2.1], 'Color', 'k'); end

figure(2)
plot(mean(gfInstF,2))
xlabel('samples')
ylabel('Instantaneous Frequency Variability [std]')
title('Global mean of between-channel Instantaneous Frequency variability')
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
for i = 1:6, line([duration(i) duration(i)], [0.9 1.3], 'Color', 'k'); end


figure(3)
plot(mean(gfRe,2))
xlabel('samples')
ylabel('Amplitudey Variability [std]')
title('Global mean of between-channel Amplitude variability')
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
for i = 1:6, line([duration(i) duration(i)], [2.0 2.5], 'Color', 'k'); end

figure(4)
plot(mean(gfPot,2))
xlabel('samples')
ylabel('Power Variability [std]')
title('Global mean of between-channel Power variability')
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
for i = 1:6, line([duration(i) duration(i)], [12.0 16], 'Color', 'k'); end

figure(5)
plot(mean(meanGFpha,2))
xlabel('samples')
ylabel('Phase Angle Variability [std]')
title('Global mean of between-channel Phase angle variability')
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
for i = 1:6, line([duration(i) duration(i)], [1.1 1.5], 'Color', 'k'); end


figure(5)
plot(mean(gfUpha,2))
xlabel('samples')
ylabel('Corrected Phase Angle Variability [std]')
title('Global mean of between-channel Corrected Phase angle variability')
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
for i = 1:6, line([duration(i) duration(i)], [0.7 0.82], 'Color', 'k'); end

% Global means

figure(6)
plot(squeeze(mean(meanEnv,3)))
xlabel('samples')
ylabel('Envelop amplitude')
title('Global mean Envelop amplitude')
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
for i = 1:6, line([duration(i) duration(i)], [3 7], 'Color', 'k'); end

figure(7)
for ch = 1:32
    subplot('Position', [position(ch,1),position(ch,2),side,side])
    set(0,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);
    h = plot(squeeze(mean(meanEnv(:,ch,:),3)));
    title(chNames(ch))
    h.LineWidth = 1.5;
    for i = 1:6, line([duration(i) duration(i)], [4 6], 'Color', 'k'); end
end


figure(8)
plot(squeeze(mean(meanIfre,3)))
xlabel('samples')
ylabel('Instantaneous frequency [Hz]')
title('Global Instantaneous Frequency')
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
for i = 1:6, line([duration(i) duration(i)], [5.8 6.5], 'Color', 'k'); end


figure(9)
plot(squeeze(mean(meanPot,3)))
xlabel('samples')
ylabel('Power [V^2]')
title('Global mean Power')
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
for i = 1:6, line([duration(i) duration(i)], [6 35], 'Color', 'k'); end

figure(9)
for ch = 1:32
    subplot('Position', [position(ch,1),position(ch,2),side,side])
    set(0,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);
    h = plot(squeeze(mean(meanPot(:,ch,:),3)));
    title(chNames(ch))
    ylabel('Power [V^2]')
    h.LineWidth = 1.5;
    for i = 1:6, line([duration(i) duration(i)], [10 15], 'Color', 'k'); end
end




figure(10)
plot(squeeze(mean(meanRe,3)))
xlabel('samples')
ylabel('Amplitude [V]')
title('Global mean Amplitude')
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
for i = 1:6, line([duration(i) duration(i)], [-3 3], 'Color', 'k'); end

figure(10)
for ch = 1:32
    subplot('Position', [position(ch,1),position(ch,2),side,side])
    set(0,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);
    h = plot(squeeze(mean(meanRe(:,ch,:),3)));
    title(chNames(ch))
    ylabel('Amplitude [V]')
    h.LineWidth = 1.5;
    for i = 1:6, line([duration(i) duration(i)], [-1 1], 'Color', 'k'); end
end



figure(11)
plot(squeeze(mean(meanPha,3)))
xlabel('samples')
ylabel('Phase [Rad]')
title('Global mean Phase')
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
for i = 1:6, line([duration(i) duration(i)], [-1 1], 'Color', 'k'); end

figure(11)
for ch = 1:32
    subplot('Position', [position(ch,1),position(ch,2),side,side])
    set(0,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);
    h = plot(squeeze(mean(meanPha(:,ch,:),3)));
    title(chNames(ch))
    ylabel('Phase [Rad]')
    h.LineWidth = 1.5;
    for i = 1:6, line([duration(i) duration(i)], [-1 1], 'Color', 'k'); end
end

figure(12)
plot(squeeze(mean(meanUpha,3)))
xlabel('samples')
ylabel('Corrected Phase [Rad]')
title('Global mean Corrected Phase')
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
for i = 1:6, line([duration(i) duration(i)], [1.3 1.8], 'Color', 'k'); end

figure(12)
for ch = 1:32
    subplot('Position', [position(ch,1),position(ch,2),side,side])
    set(0,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);
    h = plot(squeeze(mean(meanUpha(:,ch,:),3)));
    title(chNames(ch))
    ylabel('C. Phase [Rad]')
    h.LineWidth = 1.5;
    for i = 1:6, line([duration(i) duration(i)], [1.5 1.6], 'Color', 'k'); end
end

% Mean Inter all Sites Phase clusering, ISPC across all the scalp at each
% sample
figure(13)
plot(mean(meanPhaCons,2))
title('Mean Inter all Sites Phase clusering')
ylabel('mean IallSPC [0,1]')
xlabel('transition time-locked [samples]')
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
for i = 1:6, line([duration(i) duration(i)], [0.65, 0.75], 'Color', 'k'); end
