% Load band pass filtered data, do the hilbert transform per trial and save
% it

clear all
close all
clc

fileLoad = 'SegbyBand125_0558.mat';
varLoad = {'alphaEEG','betaEEG','deltaEEG','thetaEEG','gammaEEG'};
savenames = {'alphaH','betaH','deltaH','thetaH','gammaH'};
savefile = 'HilbertBands.mat';

load varForBumpsOn_Res12575.mat x y 

for vr = 1:length(varLoad)
    dataout = [];                           % put output data
    data = load(fileLoad, varLoad{vr});     % band specific data
    data = data.(varLoad{vr});
    for tr = 1:length(y)
        dataH = hilbert(data(x(tr):y(tr),:));
        dataout = vertcat(dataout, dataH);
    end
    save(savefile, 'dataout')
    nameS = savenames{vr};
    assignin('base', nameS, dataout);
    clear(varLoad{vr}, 'nameS', 'dataout', 'dataH', 'data') 
end