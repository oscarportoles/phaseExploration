% Concatenate the normed PCAs from all subjects into one data set. It also
% concatenates variables 'x' and 'y'. It also adds cells with the vaiables
% 'x5', 'y5','condsB' and 'conds'. The porpous is to use this new data
% set to run HSMM over all subjects together.

clear all
close all
clc

pathdata = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/';
snames = dir([pathdata '*epochs235.mat']);
Nsj = size(snames,1);
info = ['10 normed PCA from all subjects concatened. PCA computed individually per each subject' ...
        ' Fs = 100Hz, flitering [2, 35] Hz. x and y are the begning and end of trial respectively'];
nCondi = 4; % number fo conditions of interest (to pass)
    
data = []; 
ny = [];
nx = [];
nyEnd = 0;
nx5 = cell(1,nCondi); % New condition is not in the dataset
ny5 = cell(1,nCondi);
nconds = [];
for sj = 1:Nsj,
    filename = [pathdata snames(sj).name];
    load(filename, 'normedscore10', 'x', 'y','x5', 'y5', 'conds')
    data = vertcat(data, normedscore10);
    clear normedscore10 
    x = x + nyEnd;
    nx = vertcat(nx, x);
    y = y + nyEnd;
    ny = vertcat(ny, y);
    for c = 1:nCondi,
        x5{c} = x5{c} + nyEnd;
        y5{c} = y5{c} + nyEnd;
        nx5{1,c} = vertcat(nx5{1,c}, x5{c}); 
        ny5{1,c} = vertcat(ny5{1,c}, y5{c}); 
    end
    nconds = vertcat(nconds, conds);
    nyEnd = ny(end);
    clear x y 
end

normedscore10 = data;
clear data
y = ny;
x = nx;
conds = nconds;
x5 = nx5;
y5 = ny5;
savefile = [pathdata 'DS100epoch235allSj.mat'];
save(savefile, 'normedscore10','x','y', 'x5', 'y5', 'conds')