% It downsamples the data to 125Hz and filter the data at 58 Hz with 
% the same filter slope than the only low.
% After a visual inspection of the data, we saw some subjects had high
% amplitude 60 Hz line noise and a peak of power at 70 Hz
% I introduce a modification to replicate the same filtering parameters
% than in Qion's study but with 125 Hz Sampling frequency
%%%%

clear all
close all
clc

fs_DS = 125;

fnames = dir('/Users/lab/Documents/MATLAB/Data_JB_AssoRec/Corrected/*.set');
Nsj = size(fnames,1);

%eeglab

% for sj=1:Nsj,
%     EEG = pop_loadset('filename',fnames(sj).name ,'filepath','/Users/lab/Documents/MATLAB/Data_JB_AssoRec/Corrected/');
%     EEG = eeg_checkset( EEG );
%     EEG = pop_select( EEG,'nochannel',{'HEOG' 'VEOG'});
%     EEG = eeg_checkset( EEG );
%     EEG = pop_resample( EEG, fs_DS,0.928,0.032);
%     EEG = pop_iirfilt( EEG, 0, 58, [2], 0, 0);
%     %EEG = pop_iirfilt( EEG, 0, 30,[2], 0, 0);
%     EEG = pop_iirfilt( EEG, 0.5, 0,[0.4], 0, 0);
%     EEG.setname='DS125_0558';
%     EEG = eeg_checkset( EEG );
%     name = strsplit(fnames(sj).name,'.');
%     newname = [name{1,1} 'DS125_0558' '.set'];
%     EEG = pop_saveset( EEG, 'filename',newname,'filepath','/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125bp05_58/');
%     EEG = eeg_checkset( EEG );
%     clear EEG
% end

%%%%%%%%%%%%%--------------------------------------%%%%%%%%%%%%%%%%%%%%%%%
% Band pass filtering
% It downsamples and bandpass filters the data
%%%%%%%%%%%%%--------------------------------------%%%%%%%%%%%%%%%%%%%%%%

for sj=1:Nsj,
    EEG = pop_loadset('filename',fnames(sj).name ,'filepath','/Users/lab/Documents/MATLAB/Data_JB_AssoRec/Corrected/');
    EEG = eeg_checkset( EEG );
    EEG = pop_select( EEG,'nochannel',{'HEOG' 'VEOG'});
    EEG = eeg_checkset( EEG );
    EEGo = pop_resample( EEG, fs_DS,0.928,0.032);       % Raw downsamped EEG data
    
    % Delta band 2 - 4 Hz
    savepath = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125_EEGbands/Delta/'
    EEG = pop_iirfilt( EEGo, 0, 4, [1], 0, 0);
    EEG = pop_iirfilt( EEG, 2, 0,[1], 0, 0);
    EEG.setname='DS125_Delta';
    EEG = eeg_checkset( EEG );
    name = strsplit(fnames(sj).name,'.');
    newname = [name{1,1} 'DS125_Delta' '.set'];
    EEG = pop_saveset( EEG, 'filename',newname,'filepath', savepath);
    EEG = eeg_checkset( EEG );
    clear EEG
    % Theta band 4 - 9 Hz
    savepath = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125_EEGbands/Theta/'
    EEG = pop_iirfilt( EEGo, 0, 9, [1], 0, 0);
    EEG = pop_iirfilt( EEG, 4, 0,[1], 0, 0);
    EEG.setname='DS125_Theta';
    EEG = eeg_checkset( EEG );
    name = strsplit(fnames(sj).name,'.');
    newname = [name{1,1} 'DS125_Theta' '.set'];
    EEG = pop_saveset( EEG, 'filename',newname,'filepath', savepath);
    EEG = eeg_checkset( EEG );
    clear EEG
    % Alpha band 9 - 14 Hz
    savepath = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125_EEGbands/Alpha/'
    EEG = pop_iirfilt( EEGo, 0, 14, [1], 0, 0);
    EEG = pop_iirfilt( EEG, 9, 0,[1], 0, 0);
    EEG.setname='DS125_Alpha';
    EEG = eeg_checkset( EEG );
    name = strsplit(fnames(sj).name,'.');
    newname = [name{1,1} 'DS125_Alpha' '.set'];
    EEG = pop_saveset( EEG, 'filename',newname,'filepath', savepath);
    EEG = eeg_checkset( EEG );
    clear EEG
    % Beta band 14 - 30 Hz
    savepath = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125_EEGbands/Beta/'
    EEG = pop_iirfilt( EEGo, 0, 30, [1], 0, 0);
    EEG = pop_iirfilt( EEG, 14, 0,[1], 0, 0);
    EEG.setname='DS125_Beta';
    EEG = eeg_checkset( EEG );
    name = strsplit(fnames(sj).name,'.');
    newname = [name{1,1} 'DS125_Beta' '.set'];
    EEG = pop_saveset( EEG, 'filename',newname,'filepath', savepath);
    EEG = eeg_checkset( EEG );
    clear EEG
    % Gamma band 30 - 58 Hz
    savepath = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125_EEGbands/Gamma'
    EEG = pop_iirfilt( EEGo, 0, 58, [2], 0, 0);
    EEG = pop_iirfilt( EEG, 30, 0,[1], 0, 0);
    EEG.setname='DS125_Gamma';
    EEG = eeg_checkset( EEG );
    name = strsplit(fnames(sj).name,'.');
    newname = [name{1,1} 'DS125_Gamma' '.set'];
    EEG = pop_saveset( EEG, 'filename',newname,'filepath', savepath);
    EEG = eeg_checkset( EEG );
    clear EEG
end


