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

for sj=1:Nsj,
    EEG = pop_loadset('filename',fnames(sj).name ,'filepath','/Users/lab/Documents/MATLAB/Data_JB_AssoRec/Corrected/');
    EEG = eeg_checkset( EEG );
    EEG = pop_select( EEG,'nochannel',{'HEOG' 'VEOG'});
    EEG = eeg_checkset( EEG );
    EEG = pop_resample( EEG, fs_DS,0.928,0.032);
    %EEG = pop_iirfilt( EEG, 0, 58, [2], 0, 0);
    EEG = pop_iirfilt( EEG, 0, 30,[2], 0, 0);
    EEG = pop_iirfilt( EEG, 0.5, 0,[0.4], 0, 0);
    EEG.setname='CorrDS125';
    EEG = eeg_checkset( EEG );
    name = strsplit(fnames(sj).name,'.');
    newname = [name{1,1} 'DS125' '.set'];
    EEG = pop_saveset( EEG, 'filename',newname,'filepath','/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS125bp05_30/');
    EEG = eeg_checkset( EEG );
    clear EEG
end