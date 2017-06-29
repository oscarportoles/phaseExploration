%% Figure 2
% clear all
% %close all
% clc
% 
% load onsetStgTimeInfo5Bu.mat
% 
% % Mean duration and Standard deviation between subjects
% %meanStags = round([dur(2).stag1, dur(2).stag2, dur(2).stag3, dur(2).stag4, dur(2).stag5, dur(2).stag6]);
% stdOnStags = round([dur(4).stag1, dur(4).stag2, dur(4).stag3, dur(4).stag4, dur(4).stag5]); 
% stdOnCorStags = round([dur(5).stag1, dur(5).stag2, dur(5).stag3, dur(5).stag4, dur(5).stag5]);
% cumOnStags = round([dur(1).stag1, dur(1).stag2, dur(1).stag3, dur(1).stag4, dur(1).stag5]);
% 
% figure(3)
% y = [0.5, 1, 1.5, 2, 2.5];
% %y = ones(1,5);
% ax = gca;
% ax.YLim = [0,10];
% ax.XLim = [0, 1800];
% e = errorbar(ax,cumOnStags,y, stdOnCorStags,'horizontal','o','MarkerSize',6,...
%     'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',15, ...
%     'LineWidth', 2);
% %e.Color = 'k';
% %hold on
% % errorbar(ax,cumOnStags,y, stdOnStags,'horizontal','o','MarkerSize',6,...
% %     'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',15, ...
% %     'LineWidth', 2);
% ax.YLim = [0,3];
% ax.XLim = [0, 1400];
% ax.YTick = y;
% ax.YTickLabel = {'onS1','onS2','onS3','onS4','onS5'};
% ax.XTick = [0, cumOnStags, dur(1).avDur];
% ax.XTickLabel = {'0', num2str(cumOnStags(1)),num2str(cumOnStags(2)),num2str(cumOnStags(3)),...
%                 num2str(cumOnStags(4)),num2str(cumOnStags(5)), ...
%                 'End'};
%                 %num2str(cumStags(6))};
% %ax.XTickLabelRotation = 45;
% %ylim([0,10])
% %title('Mean stage onset & stage duration standard deviation')
% xlabel('Time [milliseconds]')

%% Fgure 4

% Phase, probability distributions and windows of cognitive events at a random trial
eventprobs = load('S8_DS100HSMM5Bout235.mat', 'eventprobs');
eventprobs = eventprobs.eventprobs;
data = load('S8_HilBads100.mat','theta', 'x','y');
dataZx = data.theta;
x = data.x;
y = data.y;
clear data

tr = 37;

yTr = y(tr) - x(tr) + 1;
t = 1:yTr;
nBins = 4;
tOn = round(t * squeeze(eventprobs(1:yTr,tr,:)));
tStg = round(movmean([0,tOn,yTr],2,'Endpoints','discard'));
%yOn = diag(squeeze(eventprobs(tOn,tr,:)))';
yOn = 0.13 * ones(1,5);
%yStg = movmean([0,yOn,0],2,'Endpoints','discard');
yStg = 0.07 * ones(1,6);
erOn = 2 * ones(size(yOn));
erStg = 2 * ones(size(yStg));

figure(tr)
%edges = linspace(-pi, pi, nBins+1);  % Does nBins from -pi to pi Bin edges
[hAx,hprob,hphase] =  plotyy(t,[squeeze(eventprobs(1:yTr,tr,:))],t,angle(dataZx(x(tr):y(tr))));
% hphase.LineStyleOrder = {'o'};
% hprob.LineStyleOrder = '*';
hphase.LineWidth = 1.4;
hphase.LineStyle = ':';
hprob(1).LineWidth = 2;
hprob(2).LineWidth = 2;
hprob(3).LineWidth = 2;
hprob(4).LineWidth = 2;
%hprob(4).LineStyle = ':';
hprob(5).LineWidth = 2;
xlabel('time [milliseconds]')
ylabel(hAx(1),'Phase [rad]') % left y-axis
hAx(2).YLim = 1.05 .* [-pi, pi];
hAx(2).YTick = [-pi,-pi/2,0,pi/2,pi];
hAx(2).YTickLabel = {'\pi','\pi/2','0','-\pi/2','-\pi'};
%hAx(1).YTick = edges;
%hAx(1).YTickLabel = {'-\pi','','','-2\pi/3','','','-\pi/3','','','0','','','\pi/3','','','2\pi/3','','','\pi'};
%hAx(2).YGrid = 'on';
hAx(1).YLim = [0, 0.35];
hAx(1).YTick = [0,0.05,0.1,0.15,0.2,0.25,0.3, 0.35];
hAx(1).YTickLabel = {'0','0.05','0.1','0.15','0.2','0.25','0.3','0.35'};
hAx(2).FontSize = 10;
hAx(1).FontSize = 10;
ylabel(hAx(1),'Probability [0, 1]') % right y-axis
%title('\theta-band phase & onset probability distributions')
%legend({'\theta phase','onS1','onS2','onS3','onS4','onS5'},'Location', 'northeast')
set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',10);

hold on
e1 = errorbar(tOn,yOn,erOn,'o','horizontal');
e1.MarkerSize = 0.1;
e1.Color = 'k';
e1.CapSize = 10;
e1.MarkerFaceColor = 'k';
e1.LineWidth = 2;

ax = gca;
gr = [0.5, 0.5, 0.5];
e2 = errorbar(ax,tStg,yStg,erStg,'o','horizontal');
e2.MarkerSize = 0.1;
e2.Color = gr;
e2.CapSize = 10;
e2.MarkerFaceColor = gr;
e2.LineWidth = 2;
ticks =  [0:10:yTr];
tickVal = {};
for i = 1:length(ticks), tickVal{i} = num2str(10*ticks(i)); end 

ax.XTick = ticks;
ax.XTickLabel = tickVal;


legend({'P(onS2)','P(onS3)','P(onS4)','P(onS5)','P(onS6)','wOnset','wStage','\theta phase'},'Location', 'northeast')
