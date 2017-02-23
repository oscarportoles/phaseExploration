% plot and save results from ITCP, Critical ITPC and Kuramoto order
% parameter for bump-lock data as well as onset lock data
clear all 
close all
clc

load electrodes.mat
load SubjItcpMax35.mat itpcBu itpcOn itpCriBu itpCriOn lockTo
load bumpLocationsMax35.mat bumpLoc 

nModl = 8;                  % number of models
nFreq = 5;                  % number of frequency bands
nSj = 20;                   % number of subjects
nCh = 32;                   % number of channels
table = {};                 % information about significant ITPC
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
                % kep track of cases with significant ITPC
                chCT = 0;
                chName = {};
                chArea = {};
                for ch = 1:nCh,
                    area = 0;
                    positive = 0;
                    for i = 4:size(data,1)-4,
                        criStd = std(data(ch,i-3:i+3,2));
                        if criStd < 0.02 & data(i,ch,1) > data(i,ch,2)     % avobe threshold to be cosidered
                            area = area + (data(i,ch,1) - data(i,ch,2));
                            positive = 1;
                        end
                    end
                    if positive,
                        chCT = chCT + 1;
                        chName{1,chCT} = locElect(ch).labels; 
                        chArea{1,chCT} = area;
                    end
                end
                if chCT >= 1,
                    infotable = {sj, freqName{fq}, mo, bu, chCT, chName, chArea};
                    table = vertcat(table, infotable);
                end
                % plot ITPC and Critical ITPC, time-lock to Onset
                if bu == 1,
                    fig(2) = figure(2);
                    data = [];
                    data(:,:,1) = itpcOn{sj,fq}';
                    data(:,:,2) = itpCriOn{sj,fq};
                    lockIx = round(mean(bumpLoc{mo}(:,bu,1)));      % mean bumps location for nodel and bump {model}(N,bump,[lockIx,t_pos, sj, condition])
                    titlen = ['Subject: ' num2str(sj) ' ,Onset-lock ITCP, Rythm: ' freqName{fq} ', Model: ' num2str(mo) ', Bump: ' num2str(bu)];
                    plottopo(data, 'chanlocs', locElect, 'ylim', [0 1],'title',titlen, 'colors', {'b','r'}, ...
                        'ydir', 1, 'vert', lockIx, 'legend', {'ITCP', 'Critic ITCP', 'Bump'}, 'showleg','on')
                end
                % save figure to a .fig
                namesave = ['sj' num2str(sj) freqName{fq} 'Modl' num2str(mo) 'Bu' num2str(bu) '.fig' ];
                savepath = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125bp05_58/events/ForBumps_OnRes_no20/LocBandsHilbert/figurITCPkuraBySubjmax35/';
                tosave = [savepath namesave ];
                savefig(fig,tosave)
                clear fig
                close all
            end
        end
    end
end
saveTable = [savepath, 'summaryITPC.mat'];
save(saveTable, 'table')