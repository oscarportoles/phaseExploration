% It creates a matrix with bump lock data. The it computes ITPC fro all
% subjects, p-value and critical ITPC
 clear all
 close all
 clc

load bumpLocations.mat

Nch = 32;                       % number of channels
nBand = 5;                      % frequency bands to analyze
nModl = 8;                      % number of modeld to analyze
nSj = 20;                       % number of subjects
nTr = length(x);                % number of trials
lockTo = cell(nSj,nModl);         % sample time locked data
itpcBu = cell(nSj,nModl,nBand,nModl);    	 % ITCP to bumps
itpcOn = cell(nSj,nBand);    		 % ITCP to Onset
pValBu = cell(nSj,nModl,nBand,nModl);    	 % p-value ITCP to bumps
pValOn = cell(nSj,nBand);    		 % p-value ITCP to Onset
itpCriBu = cell(nSj,nModl,nBand,nModl);  	 % critical ITCP to bumps 
itpCriOn = cell(nSj,nBand);  		 % critical ITCP to Onset 
kuOrderBu = cell(nSj,nModl,nBand,nModl); 	 % mean accross trials of Kuramoto order parameter, to bumps
kuOrderOn = cell(nSj,nBand); 		 % mean accross trials of Kuramoto order parameter, to Onset
for sj = 1:nSj, nTrSj(sj) = length(find(subjects == sj)); end   % nummber of trials per subject



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
                for tr = 1:nTr,
                    if subjects(tr) == sj
                        dataLk = NaN(nTrSj(sj), 600, 32);     % matrix with the time locked data [nTr, samples, Nch]
                        if not(done),
                            dataOn = NaN(nTrSj(sj), 375, 32);     % matrix with trials locked to onset
                        end
                        % Select  Indexes
                        lockIx = bumpLoc{mo}(:,bu,1);
                        ixPost = bumpLoc{mo}(:,bu,2);
                        lockTo{sj,mo}(bu) = round(mean(lockIx)) + max(lockIx);
             
                        trS = trS + 1;        % counter of number of trials per subject
                        data = dataAll(x(tr):y(tr),:);
                        pre = data(1:lockIx(tr), :);
                        pos = data(lockIx(tr)+1:end,:);
                        dataLk(trS,lockTo{sj,mo}(bu)-lockIx(tr)+1:lockTo{sj,mo}(bu),:) = pre;
                        dataLk(trS,lockTo{sj,mo}(bu)+1:lockTo{sj,mo}(bu)+ixPost(tr),:) = pos;
                        if not(done)
                            dataOn(trS,1:size(data,1),:) = data;
                        end
                        clear data pre pos
                    end
                end
                % number of trials per data point
                N = sum(~isnan(dataLk(:,:,1)),1)';

                % extract ITPC
                itpcBu{nSj,mo,fq,bu} = squeeze(abs(nanmean(exp(1i*angle(dataLk)),1)));
                % extract kuramoto order parameter
                kuOrderBu{nSj,mo,fq,bu} = squeeze(nanmean(abs(mean(exp(1i*angle(dataLk)),3)),1));
                if not(done)
                    % number of trials per data point
                    Non = sum(~isnan(dataOn(:,:,1)),1)';
                    % extract ITPC
                    itpcOn{nSj,fq} = squeeze(abs(nanmean(exp(1i*angle(dataOn)),1)));
                    kuOrderOn{nSj,fq} = squeeze(nanmean(abs(mean(exp(1i*angle(dataOn)),3)),1));           
                end
                for ch = 1:Nch
                    % ITPC p-value
                    pValBu{nSj,mo,fq,bu}(ch,:) = exp(sqrt(1 + 4.*N + 4.*(N.^2 - (N.*itpcBu{mo,fq,bu}(:,ch)).^2))-(1+2.*N));
                    % ITPC p-value approximation, (It is the same for large N)
                    % pB = exp(-1*N.*itpcBu(:,1).^2);
                    % pO = exp(-1.*Non.*itpcOn(:,1).^2);
                    % Critical ITCP, p-value with Boferroni's correction
                    pCriB = 0.05/(nnz(N)*Nch*nModl);%*factorial(nModl)); % Factorial does a extream correction
                    itpCriBu{nSj,mo,fq,bu}(ch,:) = sqrt(-1*log(pCriB) ./ N);
                    if not(done)
                        % ITPC p-value
                        pValOn{nSj,fq}(ch,:) = exp(sqrt(1 + 4.*Non + 4.*(Non.^2 - (Non.*itpcOn{fq}(:,ch)).^2))-(1+2.*Non));
                        % critical ITPC
                        pCriO = 0.05/(nnz(Non)*Nch*nModl);%*factorial(nModl));
                        itpCriOn{nSj,fq}(ch,:) = sqrt(-1*log(pCriO) ./ Non);
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
save SubjItcp.mat itpcBu pValBu itpCriBu itpcOn pValOn itpCriOn lockTo kuOrderOn kuOrderBu
