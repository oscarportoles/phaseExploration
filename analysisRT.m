% Computes and visualizes the mean reaction time and 3*STD for all subjects

close all
clear all
clc

pathdata = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/Downsamp125Lowpass58/events/';
snames = dir([pathdata '*.set']);
Nsj = size(snames,1);

m_rt = [];
sd_rt = [];
var_rt = [];

for sj=1:Nsj,
%sj = 16;
    % Open subject dataset
    EEG = pop_loadset('filename',snames(sj).name ,'filepath',pathdata);
    EEG = eeg_checkset( EEG );
    m_rt(end+1) = mean(single([EEG.event.RT]));
    sd_rt(end+1) = std(single([EEG.event.RT]));
    var_rt(end+1) = var(single([EEG.event.RT]));
end
s = [1:Nsj];

figure(1)
title('mean')
plot(m_rt)

figure(2)
title('mean + 3*SD')
plot(m_rt+3*sd_rt)

figure(3)
title('mean + 3*var')
plot(m_rt+3*var_rt)

figure(4)
plot(s,m_rt,s,sd_rt,s, m_rt+3*sd_rt, s,m_rt+3*sd_rt+500)
title('Reaction time analysis to choose epoch lenght')
legend('m RT','sd RT','m RT+3*sd RT','m RT+3*sd RT+500')
