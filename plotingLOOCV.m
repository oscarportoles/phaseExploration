figure(1)
t = 1:8;
Ix = zeros(20,4);
for s =1:20
    [~,Ix(s,1)] = max(likesL_o(s,:));
    [~,Ix(s,2)] = max(likesL_oc(s,:));
    [~,Ix(s,3)] = max(likesL_bestG(s,:));
    [~,Ix(s,4)] = max(likesL_newG(s,:));
    plot(Ix(1),t,'bo',Ix(2),t,'ko',Ix(3),t,'r*',Ix(4),t,'m*')
end

plot(Ix(:,1),'bo')
hold on
plot(Ix(:,3),'r*')
plot(Ix(:,4),'gx')
plot(Ix(:,2),'k.')
grid on

histogram(Ix(:,1))
histogram(Ix(:,2))
histogram(Ix(:,3))
histogram(Ix(:,4))