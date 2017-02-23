% instFrStd = squeeze(std(instFrN,[],2));
% plot(std(instFrStd,[],2))
% hold on
% plot(mean(instFrStd,2))
% for i = 1:6, line([duration(i) duration(i)], [5.8 6.2], 'Color', 'k'); end

fs = 100;                           % sampling frequency
tr = [6];
n_order = 10;                       % median filter parameters (Cohen X et.al. FReq sliding...)
orders = round(linspace(0.02*fs,0.4*fs,n_order)); % recommended: 10 steps between 10 (here 20ms becoause of fs) and 400 ms
orders = floor((orders-1)/2);         % pre/post halves

for tr = 6:10 
    instFr = diff(unwrap(dataPha(x(tr):y(tr),1)));
    instFrO  = instFr.* 100 ./(2*pi);
    yTr = y(tr) - x(tr) + 1;
    % Add mirrowed tails at the edges of the trial to avoid edge artifacts
    sMirr = orders(end) + 1;        % number of samples to mirrow
    instFrF = [flipud(instFr(1:sMirr)); instFr; flipud(instFr(end-sMirr:end))];
    % phasedmed = zeros(length(orders),yTr-1,nCh);
    phasedmed = zeros(length(orders),yTr-1); % Only one channel
    phasedmedF = zeros(length(orders),size(instFrF,1)); % Only one channel
    numelCt = [];
    for oi=1:n_order
        for ti=1:yTr-1
            for ch =1:1%nCh
                temp = sort(instFr(max(ti-orders(oi),1):min(ti+orders(oi),yTr-1),ch));
                phasedmed(oi,ti) = temp(floor(numel(temp)/2)+1);
                numelCt(end+1) = numel(temp);
            end
        end
        for ti=1:size(instFrF,1)
            for ch =1:1%nCh
                % this line could be simplified if it iterates from sMirr
                % to end - sMirr
                temp = sort(instFrF(max(ti-orders(oi),1):min(ti+orders(oi),size(instFrF,1)),ch));
                phasedmedF(oi,ti) = temp(floor(numel(temp)/2)+1);
                numelCt(end+1) = numel(temp);
            end
        end
    end

    instFr = fs .* squeeze(mean(phasedmed,1)) ./ (2*pi);
    instFrF = fs .* squeeze(mean(phasedmedF(:,sMirr+1:end-sMirr-1),1)) ./ (2*pi);
    
    figure(tr+10)
    plot(instFrO)
    hold on
    plot(instFr)
    plot(dataPha(x(tr):y(tr),1)+6)
    plot(instFrF)

end