
% First part does the .mat file needed for the LOOCV.
% Second part does the LOOCV
% comment the part of the code that you do not want
% The second part should run in Peregrine or any parallel computer
close all
clear all
clc


Nbumps = [1:8];
mags8 = cell(1,8);
param8 = cell(1,8);


load varForBumpsOn_Res75BL.mat
clear y5 y20 x5 x20 subjects score10 latent10 info data condsB conds coeff10

for i =1:8
    filename = ['iniCondOut10075bl_' num2str(Nbumps(i)) '_Bu.mat'];
    load(filename)
    mags8{i} = bumpMag2;
    param8{i} = gammPara2;
    
    
    clear likehood2 bumpMag2 gammPara2
end

save forLOOCV3075BL.mat mags8 x y normedscore10 subjectsF param8

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second part: Do LOOCV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%