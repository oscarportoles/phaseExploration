function plotPhaRePDF(prob, real, phase, cons, ch)
% plots the reeal and phase component of a trial overlap with the
% transition PDF

plot(prob,'LineWidth',2)
hold on
plot(phase(:,ch).*0.05)
plot(real(:,ch).*0.05)
plot(cons.*0.1)