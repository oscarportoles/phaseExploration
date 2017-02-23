% figure(1)
% for sj = 1:20,
%     plot([dur(sj).stag1,dur(sj).stag2,dur(sj).stag3,dur(sj).stag4,dur(sj).stag5,dur(sj).stag6], ...
%         'LineWidth',1)
%     hold on
% end
% plot([dur(21).stag1,dur(21).stag2,dur(21).stag3,dur(21).stag4,dur(21).stag5,dur(21).stag6], ...
%     'LineWidth',5)
% plot([dur(22).stag1,dur(22).stag2,dur(22).stag3,dur(22).stag4,dur(22).stag5,dur(22).stag6], ...
%     '-o','LineWidth',5)
% title('All Conditions')
% 
% figure(2)
% for sj = 1:20,
%     plot([dur(sj).stag1Foil,dur(sj).stag2Foil,dur(sj).stag3Foil,dur(sj).stag4Foil, ...
%         dur(sj).stag5Foil,dur(sj).stag6Foil],'LineWidth',1)
%     hold on
% end
% plot([dur(21).stag1Foil,dur(21).stag2Foil,dur(21).stag3Foil,dur(21).stag4Foil,dur(21).stag5Foil,dur(21).stag6Foil], ...
%     'LineWidth',5)
% plot([dur(22).stag1Foil,dur(22).stag2Foil,dur(22).stag3Foil,dur(22).stag4Foil,dur(22).stag5Foil,dur(22).stag6Foil], ...
%     '-o','LineWidth',5)
% title('Foil Condition')
% 
% figure(3)
% for sj = 1:20,
%     plot([dur(sj).stag1Targ,dur(sj).stag2Targ,dur(sj).stag3Targ,dur(sj).stag4Targ, ...
%         dur(sj).stag5Targ,dur(sj).stag6Targ],'LineWidth',1)
%     hold on
% end
% plot([dur(21).stag1Targ,dur(21).stag2Targ,dur(21).stag3Targ,dur(21).stag4Targ,dur(21).stag5Targ,dur(21).stag6Targ], ...
%     'LineWidth',5)
% plot([dur(22).stag1Targ,dur(22).stag2Targ,dur(22).stag3Targ,dur(22).stag4Targ,dur(22).stag5Targ,dur(22).stag6Targ], ...
%     '-o','LineWidth',5)
% title('Target Condition')
% 
% 
% figure(4)
% for sj = 1:20,
%     plot([dur(sj).stag1Targ - dur(sj).stag1Foil,dur(sj).stag2Targ - dur(sj).stag2Foil, ...
%         dur(sj).stag3Targ - dur(sj).stag3Foil,dur(sj).stag4Targ - dur(sj).stag4Foil, ...
%         dur(sj).stag5Targ - dur(sj).stag5Foil,dur(sj).stag6Targ - dur(sj).stag6Foil],'LineWidth',1)
%     hold on
% end
% plot([dur(21).stag1Targ - dur(21).stag1Foil,dur(21).stag2Targ - dur(21).stag2Foil, ...
%     dur(21).stag3Targ - dur(21).stag3Foil,dur(21).stag4Targ - dur(21).stag4Foil, ...
%     dur(21).stag5Targ - dur(21).stag5Foil,dur(21).stag6Targ - dur(21).stag6Foil],'LineWidth',5)
% plot([dur(22).stag1Targ - dur(22).stag1Foil,dur(22).stag2Targ - dur(22).stag2Foil, ...
%     dur(22).stag3Targ - dur(22).stag3Foil,dur(22).stag4Targ - dur(22).stag4Foil, ...
%     dur(22).stag5Targ - dur(22).stag5Foil,dur(22).stag6Targ - dur(22).stag6Foil],'-o','LineWidth',5)
% title('Target - Foil Conditions')

% Standard deviation between subjects
figure(5)
stdStags = [std([dur(1:20).stag1]), std([dur(1:20).stag2]), std([dur(1:20).stag3]), ...
            std([dur(1:20).stag4]), std([dur(1:20).stag5]), std([dur(1:20).stag6])];
meanStags = [mean([dur(1:20).stag1]), mean([dur(1:20).stag2]), mean([dur(1:20).stag3]), ...
            mean([dur(1:20).stag4]), mean([dur(1:20).stag5]), mean([dur(1:20).stag6])]; 

errorbar(meanStags, stdStags)
hold on
plot([dur(22).stag1,dur(22).stag2,dur(22).stag3,dur(22).stag4,dur(22).stag5,dur(22).stag6], ...
    '-o','LineWidth',2)
legend('EachSj','AllConct')
title('Duration of each stage in a 5 bumps HSMM')
ylabel('Duration [milliseconds * 10]')
xlabel('Stages')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ---------------------- %%%%%%%%%%%%%%%%%%%%%%
%%% 6 BUMPS MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ---------------------- %%%%%%%%%%%%%%%%%%%%%%

figure(1)
for sj = 1:20,
    plot([dur(sj).stag1,dur(sj).stag2,dur(sj).stag3,dur(sj).stag4,dur(sj).stag5,dur(sj).stag6,dur(sj).stag7], ...
        'LineWidth',1)
    hold on
end
plot([dur(21).stag1,dur(21).stag2,dur(21).stag3,dur(21).stag4,dur(21).stag5,dur(21).stag6,dur(21).stag7], ...
    'LineWidth',5)
plot([dur(22).stag1,dur(22).stag2,dur(22).stag3,dur(22).stag4,dur(22).stag5,dur(22).stag6,dur(22).stag7], ...
    '-o','LineWidth',5)
title('All Conditions')

figure(2)
for sj = 1:20,
    plot([dur(sj).stag1Foil,dur(sj).stag2Foil,dur(sj).stag3Foil,dur(sj).stag4Foil, ...
        dur(sj).stag5Foil,dur(sj).stag6Foil,dur(sj).stag7Foil],'LineWidth',1)
    hold on
end
plot([dur(21).stag1Foil,dur(21).stag2Foil,dur(21).stag3Foil,dur(21).stag4Foil,dur(21).stag5Foil,dur(21).stag6Foil,dur(21).stag7Foil], ...
    'LineWidth',5)
plot([dur(22).stag1Foil,dur(22).stag2Foil,dur(22).stag3Foil,dur(22).stag4Foil,dur(22).stag5Foil,dur(22).stag6Foil,dur(22).stag7Foil], ...
    '-o','LineWidth',5)
title('Foil Condition')

figure(3)
for sj = 1:20,
    plot([dur(sj).stag1Targ,dur(sj).stag2Targ,dur(sj).stag3Targ,dur(sj).stag4Targ, ...
        dur(sj).stag5Targ,dur(sj).stag6Targ,dur(sj).stag7Targ],'LineWidth',1)
    hold on
end
plot([dur(21).stag1Targ,dur(21).stag2Targ,dur(21).stag3Targ,dur(21).stag4Targ,dur(21).stag5Targ,dur(21).stag6Targ,dur(21).stag7Targ], ...
    'LineWidth',5)
plot([dur(22).stag1Targ,dur(22).stag2Targ,dur(22).stag3Targ,dur(22).stag4Targ,dur(22).stag5Targ,dur(22).stag6Targ,dur(22).stag7Targ], ...
    '-o','LineWidth',5)
title('Target Condition')


figure(4)
for sj = 1:20,
    plot([dur(sj).stag1Targ - dur(sj).stag1Foil,dur(sj).stag2Targ - dur(sj).stag2Foil, ...
        dur(sj).stag3Targ - dur(sj).stag3Foil,dur(sj).stag4Targ - dur(sj).stag4Foil, ...
        dur(sj).stag5Targ - dur(sj).stag5Foil,dur(sj).stag6Targ - dur(sj).stag6Foil, ...
        dur(sj).stag7Targ - dur(sj).stag7Foil],'LineWidth',1)
    hold on
end
plot([dur(21).stag1Targ - dur(21).stag1Foil,dur(21).stag2Targ - dur(21).stag2Foil, ...
    dur(21).stag3Targ - dur(21).stag3Foil,dur(21).stag4Targ - dur(21).stag4Foil, ...
    dur(21).stag5Targ - dur(21).stag5Foil,dur(21).stag6Targ - dur(21).stag6Foil, ...
    dur(21).stag7Targ - dur(21).stag7Foil],'LineWidth',5)
plot([dur(22).stag1Targ - dur(22).stag1Foil,dur(22).stag2Targ - dur(22).stag2Foil, ...
    dur(22).stag3Targ - dur(22).stag3Foil,dur(22).stag4Targ - dur(22).stag4Foil, ...
    dur(22).stag5Targ - dur(22).stag5Foil,dur(22).stag6Targ - dur(22).stag6Foil, ...
    dur(22).stag7Targ - dur(22).stag7Foil],'-o','LineWidth',5)
title('Target - Foil Conditions')

% Standard deviation between subjects
figure(5)
stdStags = [std([dur(1:20).stag1]), std([dur(1:20).stag2]), std([dur(1:20).stag3]), ...
            std([dur(1:20).stag4]), std([dur(1:20).stag5]), std([dur(1:20).stag6]), std([dur(1:20).stag7])];
meanStags = [mean([dur(1:20).stag1]), mean([dur(1:20).stag2]), mean([dur(1:20).stag3]), ...
            mean([dur(1:20).stag4]), mean([dur(1:20).stag5]), mean([dur(1:20).stag6]), mean([dur(1:20).stag7])]; 

errorbar(meanStags, stdStags)
hold on
plot([dur(22).stag1,dur(22).stag2,dur(22).stag3,dur(22).stag4,dur(22).stag5,dur(22).stag6,dur(22).stag7], ...
    '-o','LineWidth',2)
legend('EachSj','AllConct')
title('Duration of each stage in a 6 bumps HSMM')
ylabel('Duration [milliseconds * 10]')
xlabel('Stages')
