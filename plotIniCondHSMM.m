


%load iniTestCondOut.mat bumpMag gammPara likehood timeLaps
t = 1:10;

for bp=1:8
    figure(bp)
    plot(t,bumpMag(:,bp,1),'o',t,bumpMag(:,bp,2),'*',t,bumpMag(:,bp,3),'*',t,bumpMag(:,bp,4),'*',t,bumpMag(:,bp,5),'x',t,bumpMag(:,bp,6),'x',t,bumpMag(:,bp,7),'x')
    title('bumps magnitudes gamma = 12')
end

for bp=1:8
    figure(bp+10)
    plot(t,bumpMag2(:,bp,1),'o',t,bumpMag2(:,bp,2),'*',t,bumpMag2(:,bp,3),'*',t,bumpMag2(:,bp,4),'*',t,bumpMag2(:,bp,5),'x',t,bumpMag2(:,bp,6),'x',t,bumpMag2(:,bp,7),'x')
    title('bumps magnitudes gamma = 40')
end


t = 1:9;
figure(100)
plot(t,gammPara(:,2,1),'*',t,gammPara(:,2,2),'*',t,gammPara(:,2,3),'*',t,gammPara(:,2,4),'*',t,gammPara(:,2,5),'*',t,gammPara(:,2,6),'*',t,gammPara(:,2,7),'*')
title('gamma paramters magnitudes')
hold on
%figure(101)
plot(t,gammPara2(:,2,1),'o',t,gammPara2(:,2,2),'o',t,gammPara2(:,2,3),'o',t,gammPara2(:,2,4),'o',t,gammPara2(:,2,5),'o',t,gammPara2(:,2,6),'o',t,gammPara2(:,2,7),'o')
%title('gamma paramters magnitudes')



t = 1:9;
gammaVar = var(gammPara,0,3);
gamma2Var = var(gammPara2,0,3);
plot(t, gammaVar(:,2), t, gamma2Var(:,2))

plot(likeh)

t = 1:10;
for bp = 1:8,
    figure(bp)
    for i=1:24,
        plot(t,bumpMag(:,bp,i), '.')
        hold on
    end
%     plot(t,median(bumpMag(:,bp,1:30),3), '.')
%     hold on
%     plot(t,mean(bumpMag(:,bp,1:30),3), '*')
%     hold on
    plot(t,bumpMag(:,bp,25), '*')
    hold on
    plot(t,bumpMag(:,bp,26), 'x')
    hold on
    plot(t,bumpMag(:,bp,27), 'o')
end

t = 1:9;

figure(2)
for i=1:24,
    plot(t,gammPara2(:,2,i), '.')
    hold on
end
plot(t,gammPara2(:,2,25), '*')
hold on
plot(t,gammPara2(:,2,26), 'x')
hold on
plot(t,gammPara2(:,2,27), 'o')
hold on
%     plot(t,bumpMag(:,bp,33), 'o')
