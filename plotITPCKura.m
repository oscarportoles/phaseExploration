% plot and save results from ITCP, Critical ITPC and Kuramoto order
% parameter for bump-lock data as well as onset lock data
clear all 
close all
clc

load electrodes.mat
load SubjItcpMax.mat itpcBu itpcOn itpCriBu itpCriOn kuOrderBu kuOrderOn lockTo
load bumpLocationsMax.mat bumpLoc

nModl = 8;                  % number of models
nFreq = 5;                  % number of frequency bands
nSj = 20;                   % number of subjects
freqName = {'Delta','Theta','Alpha','Beta','Gamma'};     % frequency bands names

for sj = 1:nSj
    for mo = 1:nModl,
        for fq = 1:nFreq 
            for bu = 1:mo
                data = [];
                % plot ITPC and Critical ITPC, time-lock to Bumps
                fig(1) = figure(1);
                data(:,:,1) = itpcBu{sj,mo,fq,bu}';
                data(:,:,2) = itpCriBu{sj,mo,fq,bu};
                titlen = ['Subject: ' num2str(sj) ' ,Bump-lock ITCP, Rythm: ' freqName{fq} ', Model: ' num2str(mo) ', Bump: ' num2str(bu)];
                plottopo(data, 'chanlocs', locElect, 'ylim', [0 1],'title',titlen, 'colors', {'b','r'}, ...
                    'ydir', 1, 'vert', lockTo{mo}(bu), 'legend', {'ITCP', 'Critic ITCP', 'Bump'}, 'showleg','on')
                % plot ITPC and Critical ITPC, time-lock to Onset
                fig(2) = figure(2);
                data = [];
                data(:,:,1) = itpcOn{sj,fq}';
                data(:,:,2) = itpCriOn{sj,fq};
                lockIx = round(mean(bumpLoc{mo}(:,bu,1)));      % mean bumps location for nodel and bump {model}(N,bump,[lockIx,t_pos, sj, condition])
                titlen = ['Subject: ' num2str(sj) ' ,Onset-lock ITCP, Rythm: ' freqName{fq} ', Model: ' num2str(mo) ', Bump: ' num2str(bu)];
                plottopo(data, 'chanlocs', locElect, 'ylim', [0 1],'title',titlen, 'colors', {'b','r'}, ...
                    'ydir', 1, 'vert', lockIx, 'legend', {'ITCP', 'Critic ITCP', 'Bump'}, 'showleg','on')
                % plot Kuramoto Order parameter,  time-lock to Bumps
                fig(3) = figure(3);
                data = [];
                title1 = ['Subject: ' num2str(sj) ' ,Bump-lock Kuramoto Order, Rythm: ' freqName{fq} ', Model: ' num2str(mo) ', Bump: ' num2str(bu)];
                subplot(2,1,1)
                plot(kuOrderBu{sj,mo,fq,bu})
                title(title1)
                line([lockTo{mo}(bu)  lockTo{mo}(bu)],[0 1],'Color','k')
                % plot Kuramoto Order parameter,  time-lock to Onset           
                title2 = ['Subject: ' num2str(sj) ' ,Onset-lock Kuramoto Order, Rythm: ' freqName{fq} ', Model: ' num2str(mo) ', Bump: ' num2str(bu)];
                subplot(2,1,2)
                plot(kuOrderOn{sj,fq})
                title(title2)
                line([lockIx lockIx],[0 1],'Color','k')
                % save figure to a .fig
                namesave = ['sj' num2str(sj) freqName{fq} 'Modl' num2str(mo) 'Bu' num2str(bu) '.fig' ];
                savepath = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125bp05_58/events/ForBumps_OnRes_no20/LocBandsHilbert/figurITCPkuraBySubjmax/';
                tosave = [savepath namesave ];
                savefig(fig,tosave)
                clear fig
                close all
            end
        end
    end
end