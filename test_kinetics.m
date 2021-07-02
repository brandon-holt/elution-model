% kon = 3 * 10^4
% koff = 5.7 * 10^-5

% Method - Competing Antibody, Increasing Amount
%tstep = 1, tMax = 6000
figure;
kon = 3e4;
koff = 5.7e-5;
[time, noCompete] = ElutionModel(100, 100, 100, 0, 0, kon, koff);
plot(time, noCompete); hold on;
[time, compete10] = ElutionModel(100, 100, 100, 50, 0, kon, koff);
plot(time, compete10); hold on;
[time, compete100] = ElutionModel(100, 100, 100, 500, 0, kon, koff);
plot(time, compete100);
title('Competing Antibody Elution Method - Increasing Concentration');
xlabel('Time (s)');
ylabel('Number of Bound Signal Antibodies');
legend('No Competing Antibody', 'Competing Ab, .5X Concentration', 'Competing Ab, 5X Concentration');
