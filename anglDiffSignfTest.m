% Test if the diffrences between bump angles are significant

nBump = 5;
nCh = 32;
nSj = 20;
buPair = combnk2(1:nBump,2);
gvTest = zeros(nCh, nSj, size(buPair,1));
pValGVbu = zeros(nCh, nSj, size(buPair,1));

pathdata = '/Users/lab/Documents/MATLAB/Data_JB_AssoRec/DS100bp2_35/events/ForBumps_OnRes_no20/hilber&HSMM/';
codeHilb = ['*Bads100.mat'];
namesHilb = dir([pathdata codeHilb]);

for sj = 1:nSj,
    load([pathdata namesHilb(sj).name],'y')
    nTr = length(y);
    for ch = 1:nCh,
        for pr = 1:size(buPair,1)
            diffPha = abs(bpcTrArg(ch,buPair(pr,1),sj) - bpcTrArg(ch,buPair(pr,2),sj));
            gvTest(ch,sj,pr) = nTr .* bpcTrMod(ch,buPair(pr,1),sj) .* exp((-(diffPha).^2) ./ (4*pi/nTr)) .*(sqrt(2./nTr));
            pValGVbu(ch,sj,pr) = 1 - normcdf(gvTest(ch,sj,pr));
        end
    end
end

buPairLabels = cell(size(buPair,1),1);
for bp = 1:size(buPair,1),
    buPairLabels{bp} = ['trans. ' num2str(buPair(bp,1)), ' vs. trans. ' num2str(buPair(bp,2))];
end