% It opens the text file with the events code and matches the events to the
% events in the EEG.event structure to epoch the data in EEGlab.
% Each trial has two events: trial onset and subject's response
% Event's file has seven columns (categories):
% subject; trial number; event code; Reaction Time; accuracy; RT outlier; word pair 
% 
% BUG!! A 'boundary' event can lay between Onset and Response events on the
% same trial. If epoching is done with EEGLAB pop_epoch it takes care of
% it. But if epoching is done manually you need to take care off it. The
% scrpit 'epochForBups' deals with it.

clear all
close all
clc

pathevents = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/Corrected/';
fnames = dir([pathevents '*.txt']);
Nsj = size(fnames,1);

pathdata = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125bp05_58/';
snames = dir([pathdata '*.set']);

pathsave = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125bp05_58/events/';

for sj=1:Nsj,
%sj = 16;
    % Open subject dataset
    EEG = pop_loadset('filename',snames(sj).name ,'filepath',pathdata);
    EEG = eeg_checkset( EEG );
    % Copy event.type to a event.urtype
    [EEG.event(:).urtype] = deal(EEG.event(:).type);
    % read TXT docuemnt with events info
    name = [pathevents, fnames(sj).name];
    fileID = fopen(name);
    C = textscan(fileID,'%d %d %d %d %d %d %d');
    fclose(fileID);
    %celldisp(C);
    % columns of the document:   
    % 2:trialN, 3:code, 4:RT, 5:accuracy (T(1)/F(0)), 
    % 6:outlierRT(T(1)/F(0)), 7:wordPair
    
    for i=1:size(EEG.event, 2) % search first non bounday event
        if not(strcmp((EEG.event(i).type), 'boundary')), 
            if mod(EEG.event(i).urevent,2)  % odd
                % First trial is always an Onset marker. Onsets are odd
                tr_o=ceil(EEG.event(i).urevent/2); 
            else % even
                % Discars first non boundary trial if it is even
                tr_o=ceil(EEG.event(i).urevent/2) + 1; 
                EEG.event(i).type = 'NotFull';
            end
            break
        else
            EEG.event(i).type = 'NotFull'; 
        end
        
    end
    
%     % modify EEG.event.type for those that won't iterete later
%     if tr_o ~= 1,
%         for tr = 1:tr_o - 1,
%             EEG.event(tr*2).type = 'NotFull';
%             EEG.event(tr*2-1).type = 'NotFull';
%         end
%     end
    
    for tr = tr_o:size(C{2},1), % iterates the trials     
        if find(strcmp(int2str(C{2}(tr)*2), {EEG.event.type})) & find(strcmp(int2str(C{2}(tr)*2-1),{EEG.event.type}))
            % find indexes 
            idx1 = find(strcmp(int2str(C{2}(tr)*2-1),{EEG.event.type}));
            idx2 = find(strcmp(int2str(C{2}(tr)*2), {EEG.event.type}));

            [EEG.event(idx1).trialN, EEG.event(idx2).trialN] = deal(C{2}(tr));
            [EEG.event(idx1).code, EEG.event(idx2).code] = deal(C{3}(tr));
            [EEG.event(idx1).RT, EEG.event(idx2).RT] = deal(C{4}(tr));
            [EEG.event(idx1).accur, EEG.event(idx2).accur] = deal(C{5}(tr));
            [EEG.event(idx1).outRT, EEG.event(idx2).outRT] = deal(C{6}(tr));
            [EEG.event(idx1).wordp, EEG.event(idx2).wordp] = deal(C{7}(tr));
            EEG.event(idx1).OnRes = 1; % Onset
            EEG.event(idx2).OnRes = 0; % Response
            % Second trial
            if EEG.event(idx2).OnRes, 
                EEG.event(idx2).type = 'On';
            else
                EEG.event(idx2).type = 'Resp';
            end
            switch EEG.event(idx2).code,
                case 2
                    EEG.event(idx2).type = [EEG.event(idx2).type 'TargetFan1S'];
                case 3
                    EEG.event(idx2).type = [EEG.event(idx2).type 'TargetFan1L'];
                case 4
                    EEG.event(idx2).type = [EEG.event(idx2).type 'TargetFan2S'];
                case 5
                    EEG.event(idx2).type = [EEG.event(idx2).type 'TargetFan2L'];
                case 6
                    EEG.event(idx2).type = [EEG.event(idx2).type 'FoilFan1S'];
                case 7
                    EEG.event(idx2).type = [EEG.event(idx2).type 'FoilFan1L'];
                case 8
                    EEG.event(idx2).type = [EEG.event(idx2).type 'FoilFan2S'];
                case 9
                    EEG.event(idx2).type = [EEG.event(idx2).type 'FoilFan2L'];
                case 10
                    EEG.event(idx2).type = [EEG.event(idx2).type 'NewFan1S'];
                case 11
                    EEG.event(idx2).type = [EEG.event(idx2).type 'NewFan1L'];
            end
            % First trial
            if EEG.event(idx1).OnRes
                EEG.event(idx1).type = 'On';
            else
                EEG.event(idx1).type = 'Resp';
            end
            switch EEG.event(idx1).code,
                case 2
                    EEG.event(idx1).type = [EEG.event(idx1).type 'TargetFan1S'];
                case 3
                    EEG.event(idx1).type = [EEG.event(idx1).type 'TargetFan1L'];
                case 4
                    EEG.event(idx1).type = [EEG.event(idx1).type 'TargetFan2S'];
                case 5
                    EEG.event(idx1).type = [EEG.event(idx1).type 'TargetFan2L'];
                case 6
                    EEG.event(idx1).type = [EEG.event(idx1).type 'FoilFan1S'];
                case 7
                    EEG.event(idx1).type = [EEG.event(idx1).type 'FoilFan1L'];
                case 8
                    EEG.event(idx1).type = [EEG.event(idx1).type 'FoilFan2S'];
                case 9
                    EEG.event(idx1).type = [EEG.event(idx1).type 'FoilFan2L'];
                case 10
                    EEG.event(idx1).type = [EEG.event(idx1).type 'NewFan1S'];
                case 11
                    EEG.event(idx1).type = [EEG.event(idx1).type 'NewFan1L'];
            end 
        else
            if find(strcmp(int2str(C{2}(tr)*2), {EEG.event.type})),
                idx2 = find(strcmp(int2str(C{2}(tr)*2), {EEG.event.type}));
                EEG.event(idx2).type = 'NotFull'; 
            end
            if find(strcmp(int2str(C{2}(tr)*2-1), {EEG.event.type})), 
                idx1 = find(strcmp(int2str(C{2}(tr)*2-1),{EEG.event.type}));
                EEG.event(idx1).type = 'NotFull'; 
            end
%             if not(find(strcmp(int2str(C{2}(tr)*2), {EEG.event.type}))) & ...
%                     not(find(strcmp(int2str(C{2}(tr)*2-1), {EEG.event.type}))),
%                 
%                 EEG.event(tr*2).type = 'NotFull';
%                 EEG.event(tr*2-1).type = 'NotFull';
%            end
        end
    end
    name = strsplit(snames(sj).name,'_');
    newname = [name{1,1} '_DS125Events.set'];
    EEG = pop_saveset( EEG, 'filename',newname,'filepath',pathsave);
    EEG = eeg_checkset( EEG );
    clear EEG
end
        
    