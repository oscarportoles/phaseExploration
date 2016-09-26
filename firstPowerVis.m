% Open a subject and compute the power spectrum with visaulization
% purpuouses


fnames = dir('Data_JB_AssoRec/Corrected/*.set');
Nsj = size(fnames,1);

for sj=1:Nsj,
    EEG = pop_loadset('filename',fnames(sj).name ,'filepath','/Users/lab/Documents/MATLAB/Data_JB_AssoRec/Corrected/');
    EEG = eeg_checkset( EEG );
    figure; pop_spectopo(EEG, 1, [0      3050000], 'EEG' , 'percent', 15, 'freq', [60 10 70], 'freqrange',[0 125],'electrodes','off');
end