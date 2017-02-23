
% Calculate the gamma paramter as well as the Bumps magnitides. The proces
% is repetated multiple time with different initial conditions. The median
% of multiple initializtions is taken as the best parameter to avoud local
% minimum. It makes use of the parallel computing peregirne cluster.

load 'HSMMsegments125_0558.mat'

% initial condition gamma parameter
rep = 500;
Nbumps = 1;
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
eventprob20i = zeros(375, 19297, Nbumps, 20);
eventprob20iGF = zeros(375, 19297, Nbumps, 20);

% number of tasks/cores --> ntask in sbatch file
parpool(24);

parfor r=1:rep,
    % Gamma 40
    [~,bumpMag(:,:,r),gammPara(:,:,r), ~]=hsmmEEG(normedscore10,magIni(:,:,r),gammaIni,1,x,y);
end

poolobj = gcp('nocreate');
delete(poolobj);

% Best inital conditions
[likehood2,bumpMag2,gammPara2, eventprobs2]= ...
        hsmmEEG(normedscore10,median(bumpMag,3),median(gammPara,3),1,x,y);

% Correct bump location for individual differences by subject
parpool(24);
parfor i = 1:20 
    [lk20i(i,1),params20i(:,:,i), eventprob20i(:,:,:,i)] = ...
        hsmmEEGfixMags(normedscore10(subjectsF==i,:),bumpMag2,gammPara2,1,x20{i},y20{i});
    
    [lk20iGF(i,1),params20iGF(:,:,i), eventprob20iGF(:,:,:,i)] = ...
        hsmmEEGfixMags(normedscore10(subjectsF==i,:),bumpMag2,gammaIni,1,x20{i},y20{i}); 
end

%save 
filename = ['stages125_0558_' num2str(Nbumps) '_Bu.mat'];
save(filename, 'likehood2', 'bumpMag2', 'gammPara2', 'eventprobs2', ...
    'lk20i', 'lk20iGF', 'params20i', 'params20iGF', 'eventprob20i','eventprob20iGF')
