% It creates a matrix with bump lock data. Then it computes ITPC for all
% subjects, p-value and critical ITPC
 clear all
 close all
 clc

load bumpLocationsMax35.mat
load bumpLocations.mat x y

Nch = 32;                       % number of channels
nBand = 5;                      % frequency bands to analyze
nModl = 8;                      % number o4f modeld to analyze
nSj = 20;                       % number of subjects
nTr = length(x);                % number of trials
lockTo = cell(nModl);         % sample time locked data
itpcBu = cell(nSj,nModl,nBand,nModl);    	 % ITCP to bumps
itpcOn = cell(nSj,nBand);    		 % ITCP to Onset
pValBu = cell(nSj,nModl,nBand,nModl);    	 % p-value ITCP to bumps
pValOn = cell(nSj,nBand);    		 % p-value ITCP to Onset
itpCriBu = cell(nSj,nModl,nBand,nModl);  	 % critical ITCP to bumps 
itpCriOn = cell(nSj,nBand);  		 % critical ITCP to Onset 
%kuOrderBu = cell(nSj,nModl,nBand,nModl); 	 % mean accross trials of Kuramoto order parameter, to bumps
%kuOrderOn = cell(nSj,nBand); 		 % mean accross trials of Kuramoto order parameter, to Onset

for fq = 1:nBand
    if fq == 1, load('deltaH.mat'); 
    elseif fq == 2, load('thetaH.mat'); 
    elseif fq == 3, load('alphaH.mat'); 
    elseif fq == 4, load('betaH.mat'); 
    elseif fq == 5, load('gammaH.mat');
    end
    for sj = 1:nSj,
        done = 0;                       % Trials lock to the onset and ITCP vlaues calculated
        for mo = 1:nModl
            for bu = 1:mo
                trS = 0;        % counter of number of trials per subject
                % Select  Indexes
                lockIx = bumpLoc{mo}(:,bu,1);
                ixPost = bumpLoc{mo}(:,bu,2);
                lockTo{mo}(bu) = max(lockIx);
                dataLk = NaN(noNaNtr(mo+1), 750, Nch);     % matrix with the time locked data [nTr, samples, Nch]
                if not(done),
                    dataOn = NaN(noNaNtr(mo+1), 375, Nch);     % matrix with trials locked to onset
                end
                for tr = 1:nTr,
                    if (bumpLoc{mo}(tr,bu,3) == sj) & not(isnan(bumpLoc{mo}(tr,bu,1)))       
                        trS = trS + 1;        % counter of number of trials per subject
                        data = dataAll(x(tr):y(tr),:);
                        pre = data(1:lockIx(tr), :);
                        pos = data(lockIx(tr)+1:end,:);
                        dataLk(trS,lockTo{mo}(bu)-lockIx(tr)+1:lockTo{mo}(bu),:) = pre;
                        dataLk(trS,lockTo{mo}(bu)+1:lockTo{mo}(bu)+ixPost(tr),:) = pos;
                        if not(done)
                            dataOn(trS,1:size(data,1),:) = data;
                        end
                        clear data pre pos
                    end
                end
%                 if not(done)
%                     dataOn = dataOn(1:trS,:,:);
%                 end
%                 dataLk = dataLk(1:trS,:,:);

                % number of trials per data point
                N = sum(~isnan(dataLk(:,:,1)),1)';      % [nTr,1]
                % extract ITPC
                itpcBu{sj,mo,fq,bu} = squeeze(abs(nanmean(exp(1i.*angle(dataLk)),1)));
                % extract kuramoto order parameter
                %kuOrderBu{sj,mo,fq,bu} = squeeze(nanmean(abs(mean(exp(1i*angle(dataLk)),3)),1));
                if not(done)
                    % number of trials per data point
                    
                    Non = sum(~isnan(dataOn(:,:,1)),1)';
                    % extract ITPC
                    itpcOn{sj,fq} = squeeze(abs(nanmean(exp(1i.*angle(dataOn)),1)));
                    %kuOrderOn{sj,fq} = squeeze(nanmean(abs(mean(exp(1i*angle(dataOn)),3)),1));           
                end
                for ch = 1:Nch
                    % ITPC p-value
                    pValBu{sj,mo,fq,bu}(ch,:) = exp(sqrt(1 + 4.*N + 4.*(N.^2 - (N.*itpcBu{sj,mo,fq,bu}(:,ch)).^2))-(1+2.*N));
                    % ITPC p-value approximation, (It is the same for large N)
                    % pB = exp(-1*N.*itpcBu(:,1).^2);
                    % pO = exp(-1.*Non.*itpcOn(:,1).^2);
                    % Critical ITCP, p-value with Boferroni's correction
                    pCriB = 0.05/(300*Nch*nModl);%*factorial(nModl)); % Factorial does a extream correction
                    itpCriBu{sj,mo,fq,bu}(ch,:) = sqrt(-1*log(pCriB) ./ N);
                    if not(done)
                        % ITPC p-value
                        pValOn{sj,fq}(ch,:) = exp(sqrt(1 + 4.*Non + 4.*(Non.^2 - (Non.*itpcOn{sj,fq}(:,ch)).^2))-(1+2.*Non));
                        % critical ITPC
                        pCriO = 0.05/(300*Nch*nModl);%*factorial(nModl));
                        itpCriOn{sj,fq}(ch,:) = sqrt(-1*log(pCriO) ./ Non);
                    end
                end
                if bu == 1 && mo == 1, 
                    done = 1;
               	    clear dataOn
                end
            end
        end
    end
    clear dataAll
end % frequency bands end
save SubjItcpMax35.mat itpcBu pValBu itpCriBu itpcOn pValOn itpCriOn lockTo %kuOrderOn kuOrderBu
