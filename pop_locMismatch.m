function [answer] = pop_locMismatch(tr,lens,eventprob20i,locId)
% It pops up a user interactions window If two bumps have the same location.
% The user can specify the best location for both bumps. It is called from
% searchBumpLocation.m
wrong = 1;
while wrong
    t = 1:lens(tr);
    figure(1)
    plot(t,eventprob20i(:,tr,bu-1),'o-', t,eventprob20i(:,tr,bu),'o-')
    hold on
    plot(locId(tr,bu-1),eventprob20i(:,tr,bu-1),locId(tr,bu),eventprob20i(:,tr,bu))
    show = ['first bump (s): ', num2str(locId(tr,bu-1)),'second bump (s): ',num2str(locId(tr,bu))];
    title(show)
    text1 = ['first bump (s): ' num2str(locId(tr,bu-1))];
    text2 = ['second bump (s): ' num2str(locId(tr,bu))];
    prompt = {text1,text2};
    name = 'Bump locations';
    answer = inputdlg(prompt,name,[1 40]);
    close(1)
    if ~isempty(answer) & ~isempty(answer{1}) & ~isempty(answer{2})
        wrong = 0; % checks that all field are filled
    end
end
end