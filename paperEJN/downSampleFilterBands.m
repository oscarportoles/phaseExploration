% It downsamples raw EEG data to 100 Hz and filter the data into EEG freqency bands
% it requires EEGLAB.
% You need to declare the file path where you have the raw data and the path to save
% the resuls. 
%%%%

clear all
close all
clc
 
fs_DS = 100; % new sampling frequency
savepath = '/Data_JB_AssoRec/DS100_EEGbands/';

fnames = dir('/Data_JB_AssoRec/Corrected/*.set');
Nsj = size(fnames,1); % number of subjects



%%%%%%%%%%%%%--------------------------------------%%%%%%%%%%%%%%%%%%%%%%%
% Band pass filtering
% It downsamples and bandpass filters the data
%%%%%%%%%%%%%--------------------------------------%%%%%%%%%%%%%%%%%%%%%%

for sj=1:Nsj,
    % Load raw data
    EEG = pop_loadset('filename',fnames(sj).name ,'filepath','/Data_JB_AssoRec/Corrected/'); 
    EEG = eeg_checkset( EEG );
    EEG = pop_select( EEG,'nochannel',{'HEOG' 'VEOG'});
    EEG = eeg_checkset( EEG );
    EEGo = pop_resample( EEG, fs_DS);       % Raw downsamped EEG data
    
%     % Delta band 2 - 4 Hz
%     savepath = '/Data_JB_AssoRec/DS100_EEGbands/Delta/'
%     EEG = pop_iirfilt( EEGo, 0, 4, [1], 0, 0);
%     EEG = pop_iirfilt( EEG, 2, 0,[1], 0, 0);
% %     EEG = pop_eegfiltnew(EEG, 2, [], [], 0, [], 0); % low cutoff
% %     EEG = pop_eegfiltnew(EEG, [], 35, [], 0, [], 0); % high cutoff
%     EEG.setname='DS100_Delta';
%     EEG = eeg_checkset( EEG );
%     name = strsplit(fnames(sj).name,'.');
%     newname = [name{1,1} 'DS100_Delta' '.set'];
%     EEG = pop_saveset( EEG, 'filename',newname,'filepath', savepath);
%     EEG = eeg_checkset( EEG );
%     clear EEG
    % Theta band 4 - 9 Hz
    savepath = '/Data_JB_AssoRec/DS100_EEGbands/Theta/'
    EEG = pop_iirfilt( EEGo, 0, 9, [1], 0, 0);
    EEG = pop_iirfilt( EEG, 4, 0,[1], 0, 0);
%     EEG = pop_eegfiltnew(EEG, 2, [], [], 0, [], 0); % low cutoff
%     EEG = pop_eegfiltnew(EEG, [], 35, [], 0, [], 0); % high cutoff
    EEG.setname='DS100_Theta';
    EEG = eeg_checkset( EEG );
    name = strsplit(fnames(sj).name,'.');
    newname = [name{1,1} 'DS100_Theta' '.set'];
    EEG = pop_saveset( EEG, 'filename',newname,'filepath', savepath);
    EEG = eeg_checkset( EEG );
    clear EEG
    % Alpha band 9 - 14 Hz
    savepath = '/Data_JB_AssoRec/DS100_EEGbands/Alpha/'
    EEG = pop_iirfilt( EEGo, 0, 14, [1], 0, 0);
    EEG = pop_iirfilt( EEG, 9, 0,[1], 0, 0);
%     EEG = pop_eegfiltnew(EEG, 2, [], [], 0, [], 0); % low cutoff
%     EEG = pop_eegfiltnew(EEG, [], 35, [], 0, [], 0); % high cutoff
    EEG.setname='DS100_Alpha';
    EEG = eeg_checkset( EEG );
    name = strsplit(fnames(sj).name,'.');
    newname = [name{1,1} 'DS100_Alpha' '.set'];
    EEG = pop_saveset( EEG, 'filename',newname,'filepath', savepath);
    EEG = eeg_checkset( EEG );
    clear EEG
    % Beta band 14 - 30 Hz
    savepath = '/Data_JB_AssoRec/DS100_EEGbands/Beta/'
    EEG = pop_iirfilt( EEGo, 0, 30, [1], 0, 0);
    EEG = pop_iirfilt( EEG, 14, 0,[1], 0, 0);
%     EEG = pop_eegfiltnew(EEG, 2, [], [], 0, [], 0); % low cutoff
%     EEG = pop_eegfiltnew(EEG, [], 35, [], 0, [], 0); % high cutoff
    EEG.setname='DS100_Beta';
    EEG = eeg_checkset( EEG );
    name = strsplit(fnames(sj).name,'.');
    newname = [name{1,1} 'DS100_Beta' '.set'];
    EEG = pop_saveset( EEG, 'filename',newname,'filepath', savepath);
    EEG = eeg_checkset( EEG );
    clear EEG
%     % Gamma band 30 - 58 Hz
%     savepath = '/Data_JB_AssoRec/DS100_EEGbands/Gamma'
%     EEG = pop_iirfilt( EEGo, 0, 58, [2], 0, 0);
%     EEG = pop_iirfilt( EEG, 30, 0,[1], 0, 0);
% %     EEG = pop_eegfiltnew(EEG, 2, [], [], 0, [], 0); % low cutoff
% %     EEG = pop_eegfiltnew(EEG, [], 35, [], 0, [], 0); % high cutoff
%     EEG.setname='DS100_Gamma';
%     EEG = eeg_checkset( EEG );
%     name = strsplit(fnames(sj).name,'.');
%     newname = [name{1,1} 'DS100_Gamma' '.set'];
%     EEG = pop_saveset( EEG, 'filename',newname,'filepath', savepath);
%     EEG = eeg_checkset( EEG );
%     clear EEG
end


