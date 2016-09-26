
% Calculate the gamma paramter as well as the Bumps magnitides. The proces
% is repetated multiple time with different initial conditions. The median
% of multiple initializtions is taken as the best parameter to avoud local
% minimum. It makes use of the parallel computing peregirne cluster.

load 'HSMMsegments.mat'

% initial condition gamma parameter
%gammaIni40 = [2, 40;2, 40;2, 40;2, 40;2, 40;2, 40;2, 40;2, 40;2, 40];

rep = 50;
Nbumps = 4;

halfT = 188; % half of data points in the longest trial (375)

%gammaIni = gammaIni40(1:Nbumps+1,:);
gammaIni=repmat([2 halfT/Nbumps],Nbumps+1,1);

% Initial conditions bumps magnitudes
magIni = zeros(10,Nbumps,rep);

for i=1:rep,
    magIni(:,:,i) = -1 + 2*rand(10,Nbumps); 
end

% output parameters gamma = 40
%likehood = zeros(1,rep);
bumpMag = zeros(size(magIni));
gammPara = zeros(Nbumps+1,2,rep);
%eventprobs = zeros(375, 19297, Nbumps, rep);

% number of tasks/cores --> ntask in sbatch file
parpool(24);

parfor r=1:rep,
    % Gamma 40
    [~,bumpMag(:,:,r),gammPara(:,:,r), ~]=hsmmEEG(normedscore10,magIni(:,:,r),gammaIni,1,x,y);
end

% gamma 40
[likehood2,bumpMag2,gammPara2, eventprobs2]= ...
        hsmmEEG(normedscore10,median(bumpMag,3),median(gammPara,3),1,x,y);


    
%save iniCondOutG12_8Bu.mat likehood bumpMag gammPara
filename = ['iniCondOutG40_' num2str(Nbumps) '_Bu.mat'];
save(filename, 'likehood2', 'bumpMag2', 'gammPara2', 'eventprobs2')
