% EStimates bumps magnitudes and gamma parameter per subjects. It will
% attemts to remove part of the variability between subjects.

clear all
close all
clc

nsamp = 375;        % max number of samples in a trial 
ntrial = 10400;     % max number of trials seen by a subject
nsubj = 20;         % number of subjects
nbump = 8;          % number of bumps
%eventprobs = zeros(nsamp, ntrial, 1, nsubj,nbump);  % prob. bump location
lkP20 = zeros(nsubj, 1,nbump);          % likelikood
mags20 = zeros(10,nbump,nsubj,nbump);   % bump's maginitudes
params20 = zeros(10,2,nsubj,nbump);     % gamma parameters

load('HSMMsegments125_0158.mat')
load('Mags8Params8.mat')
%parpool(nsubj)
%parfor i = 1:nsubj
for i = 1:nsubj
    data = normedscore10(subjectsF==i,:);
    for b = 1:nbump
        [lkP20(i,b),mags20(:,:,i,b),params20(:,:,i,b),eventprob20(:,:,:,i,b)] = ...
            hsmmEEG(data,mags8{b},params8{b},1,x20{i},y20{i});
    end
end
save subjectOutput.mat lkP20 mags20 params20 eventprob20 -v7.3

% parfor i = 1:20 
%     [lkP20(i,1),params20(:,:,i)]= hsmmEEGfixMags(normedscore10(subjectsF==i,:),mags8{5},params8{5},1,x20{i},y20{i}); 
% end
% 
% parfor i = 1:20
%     [lkM20(i,1),mags20(:,:,i)] = hsmmEEGfixParams(normedscore10(subjectsF==i,:),mags8{5},params8{5},1,x20{i},y20{i}); 
% end

%%%%%
% HEre it does two variables with the best magnitudes and gamma scale parameters
% mags8 = cell(8,1);
% params8 = cell(8,1);
% 
% for f =1:8
%     fileload = ['iniCondOut125_0158_' num2str(f) '_Bu.mat'];
%     load(fileload, 'bumpMag2', 'gammPara2');
%     mags8{f} = bumpMag2;
%     params8{f} = gammPara2;
% end
% save Mags8Params8.mat mags8 params8
%     
    
    