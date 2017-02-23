% It creates a matrix with bump lock data. The it computes ITPC fro all
% subjects, p-value and critical ITPC
 clear all
 close all
 clc

load bumpLocations.mat

Nch = 32;                       % number of channels
nBand = 5;                      % frequency bands to analyze
nModl = 8;                      % number of modeld to analyze
nTr = length(x);                % number of trials
lockTo = cell(nModl,1);         % sample time locked data
itpcBu = cell(nModl,nBand,nModl);    	 % ITCP to bumps
itpcOn = cell(1,nBand);    		 % ITCP to Onset
pValBu = cell(nModl,nBand,nModl);    	 % p-value ITCP to bumps
pValOn = cell(1,nBand);    		 % p-value ITCP to Onset
itpCriBu = cell(nModl,nBand,nModl);  	 % critical ITCP to bumps 
itpCriOn = cell(1,nBand);  		 % critical ITCP to Onset 
kuOrderBu = cell(nModl,nBand,nModl); 	 % mean accross trials of Kuramoto order parameter, to bumps
kuOrderOn = cell(1,nBand); 		 % mean accross trials of Kuramoto order parameter, to Onset

for fq = 1:nBand
    if fq == 1, load('deltaH.mat'); 
    elseif fq == 2, load('thetaH.mat'); 
    elseif fq == 3, load('alphaH.mat'); 
    elseif fq == 4, load('betaH.mat'); 
    elseif fq == 5, load('gammaH.mat'); 
    end
    done = 0;                       % Trials lock to the onset and ITCP vlaues calculated
    for mo = 1:nModl
        for bu = 1:mo
            dataLk = NaN(nTr, 600, 32);     % matrix with the time locked data [nTr, samples, Nch]
            if not(done),
                dataOn = NaN(nTr, 375, 32);     % matriz with trials locked to onset
            end
	    % Select  Indexes
            lockIx = bumpLoc{mo}(:,bu,1);
            ixPost = bumpLoc{mo}(:,bu,2);
            lockTo{mo}(bu) = round(mean(lockIx)) + max(lockIx);
            for tr = 1:nTr,
                data = dataAll(x(tr):y(tr),:);
                pre = data(1:lockIx(tr), :);
                pos = data(lockIx(tr)+1:end,:);
                dataLk(tr,lockTo{mo}(bu)-lockIx(tr)+1:lockTo{mo}(bu),:) = pre;
                dataLk(tr,lockTo{mo}(bu)+1:lockTo{mo}(bu)+ixPost(tr),:) = pos;
                if not(done)
                    dataOn(tr,1:size(data,1),:) = data;
                end
            end
            clear data pre pos

            % number of trials per data point
            N = sum(~isnan(dataLk(:,:,1)),1)';
            
            % extract ITPC
            itpcBu{mo,fq,bu} = squeeze(abs(nanmean(exp(1i*angle(dataLk)),1)));
            % extract kuramoto order parameter
            kuOrderBu{mo,fq,bu} = squeeze(nanmean(abs(mean(exp(1i*angle(dataLk)),3)),1));
            if not(done)
                % number of trials per data point
                Non = sum(~isnan(dataOn(:,:,1)),1)';
                % extract ITPC
                itpcOn{fq} = squeeze(abs(nanmean(exp(1i*angle(dataOn)),1)));
                kuOrderOn{fq} = squeeze(nanmean(abs(mean(exp(1i*angle(dataOn)),3)),1));           
            end
            for ch = 1:Nch
                % ITPC p-value
                pValBu{mo,fq,bu}(ch,:) = exp(sqrt(1 + 4.*N + 4.*(N.^2 - (N.*itpcBu{mo,fq,bu}(:,ch)).^2))-(1+2.*N));
                % ITPC p-value approximation, (It is the same for large N)
                % pB = exp(-1*N.*itpcBu(:,1).^2);
                % pO = exp(-1.*Non.*itpcOn(:,1).^2);
                % Critical ITCP, p-value with Boferroni's correction
                pCriB = 0.05/(nnz(N)*Nch*nModl);%*factorial(nModl)); % Factorial does a extream correction
                itpCriBu{mo,fq,bu}(ch,:) = sqrt(-1*log(pCriB) ./ N);
                if not(done)
                    % ITPC p-value
                    pValOn{fq}(ch,:) = exp(sqrt(1 + 4.*Non + 4.*(Non.^2 - (Non.*itpcOn{fq}(:,ch)).^2))-(1+2.*Non));
                    % critical ITPC
                    pCriO = 0.05/(nnz(Non)*Nch*nModl);%*factorial(nModl));
                    itpCriOn{fq}(ch,:) = sqrt(-1*log(pCriO) ./ Non);
                end
            end
            if bu == 1 && mo == 1, 
                done = 1;
     		clear dataOn
            end
        end
    end
    clear dataAll
end % frequency bands end
save itcp.mat itpcBu pValBu itpCriBu itpcOn pValOn itpCriOn lockTo kuOrderOn kuOrderBu
