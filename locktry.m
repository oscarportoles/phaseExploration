clear all
data = [1,1,2,1,1;
        1,2,1,1,1;
        1,1,1,2,1;
        2,1,1,1,1;
        1,1,1,1,2];

ix = [3,2,4,1,5];
loI = floor(mean(ix))+min(ix)+1;
 
lock = zeros(size(data,1),size(data,2)*2);

for tr=1:size(data,1)
    pre = data(tr,1:ix(tr));
    pos = data(tr,ix(tr)+1:end);
    lock(tr,loI-length(pre)+1:loI) = pre;
    lock(tr,loI+1:loI+length(pos)) = pos;
end