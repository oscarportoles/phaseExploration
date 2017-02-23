
% Calculate the gamma paramter as well as the Bumps magnitides. The proces
% is repetated mfor each subject. The median
% of multiple initializtions is taken as the best parameter to avoud local
% minimum. It makes use of the parallel computing peregirne cluster. It
% needs a the function parsave to save the datasets

snames = dir('*135.mat');
Nsj = size(snames,1);
rep = 500;
Nbumps = 5;
maxDur = 300;


% initial condition gamma parameter
gammaIni=repmat([2 maxDur/Nbumps],Nbumps+1,1);
% Initial conditions bumps magnitudes
magIni = zeros(10,Nbumps,rep);
for i=1:rep,
    magIni(:,:,i) = -1 + 2*rand(10,Nbumps); 
end
% number of tasks/cores --> ntask in sbatch file
parpool(Nsj);

parfor sj= 1:Nsj
    [data,x,y] = parload(snames(sj).name);
    % output parameters 
    bumpMag = zeros(size(magIni));
    gammPara = zeros(Nbumps+1,2,rep);
    %eventprobs = zeros(maxDur, 19297, Nbumps); % eventprobs[Max Scamples a trial, N of trials, N of bump]

    for r=1:rep,
        [~,bumpMag(:,:,r),gammPara(:,:,r), ~] = hsmmEEG_o(data,magIni(:,:,r),gammaIni,1,x,y);
    end

    % Best inital conditions
    [likehood,bumpMag,gammPara, eventprobs]= ...
            hsmmEEG_o(data,median(bumpMag,3),median(gammPara,3),1,x,y);

    %save 
    filename = strrep(snames(sj).name,'epochs','HSMMout');
    parsave(filename, likehood, bumpMag, gammPara, eventprobs)
end
