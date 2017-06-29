% 
% T4=zeros(32,32);
% XX=zeros(20,1);
% for i=1:32
%     for j=1:32
%         counter=0;
%         XX=zeros(20,1);
%         for z=1:20
%             K=csj(z).t4;
%             if K(i,j)~=0
%                 counter=counter+1;
%             end
%         end
%         if counter>=14
%             for z=1:20
%             XX(z)=csj(z).t4(i,j);
%             end
%         XX=nonzeros(XX);
%         T4(i,j)=median(XX);
%        
%         end
%     end
% end

nCh = 32;
nSj = 20;
pairs = combnk2([1:nCh],2);
th = 14;

cAll = zeros(nCh,nCh,nSj);
cGlo = zeros(nCh,nCh);
gPair = [];

for sj = 1:nSj
    cAll(:,:,sj) = csj(sj).t1;
end

for pr = 1:size(pairs,1)
    if sum(logical(cAll(pairs(pr,1),pairs(pr,2),:))) >= th
        gPair(end+1,:) = [pairs(pr,1), pairs(pr,2)];
    end
end
ds.chanPairs = gPair;
topoplot_connect(ds, locElect)