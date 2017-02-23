% Load band pass filtered data, do the hilbert transform per trial and save
% it

clear all
close all
clc

pathdata = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/';
snames = dir([pathdata '*Bands100.mat']);
Nsj = length(snames);
varLoad = {'alphaEEG','betaEEG','thetaEEG'};
varsave = {'alpha', 'beta', 'theta'};
varLoadExtra = {'x','x5','y','y5','conds','condsB','subject'};
info = {'Hilber transformed data (analytic signal) by EEG frequency bands, Fs = 100Hz'};
% varLoad = {'alphaEEG','betaEEG','deltaEEG','thetaEEG','gammaEEG'};
% savenames = {'alphaH','betaH','deltaH','thetaH','gammaH'};
savepath = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/';

for sj = 1:Nsj,
    for vr = 1:length(varLoad)
        dataout = [];                           % put output data
        fileload = [pathdata snames(sj).name];
        data = load(fileload, varLoad{vr});     % band specific data
        data = data.(varLoad{vr});
        load(fileload, varLoadExtra{:});
        for tr = 1:length(y)
            dataH = hilbert(data(x(tr):y(tr),:));
            dataout = vertcat(dataout, dataH);
        end
        eval([varsave{vr} '= dataout;'])
        clear('dataout', 'dataH', 'data') 
    end
    idSj = regexp(snames(sj).name, '__');
    namesave = [snames(sj).name(1:idSj) 'HilBads100.mat'];
    savefile = [savepath namesave];
    save(savefile, varsave{1}, varsave{2}, varsave{3}, 'x', 'y', 'x5', 'y5', 'conds', 'condsB', 'subject')
    clear x y x5 y5 conds condsB subject
end