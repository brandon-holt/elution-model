% testing
clear;clc;close all;

% Method - Competing Antibody, Increasing Affinity
figure;
kon = [.1, .1, .1];
koff = [.1, .1, .1];
[time, noCompete] = ElutionModel(100, 100, 100, 0, 0, kon, koff);
plot(time, noCompete); hold on;
kon = [.1, .1, .2];
koff = [.1, .1, .05];
[time, compete10] = ElutionModel(100, 100, 100, 100, 0, kon, koff);
plot(time, compete10); hold on;
kon = [.1, .1, 10];
koff = [.1, .1, .001];
[time, compete100] = ElutionModel(100, 100, 100, 100, 0, kon, koff);
plot(time, compete100);
title('Competing Antibody Elution Method - Increasing Affinity');
xlabel('Time (s)');
ylabel('Number of Bound Signal Antibodies');
legend('No Competing Antibody', 'Competing Ab, 2X Kinetics', 'Competing Ab, 100X Kinetics');

% Method - Competing Antibody, Increasing Amount
figure;
kon = [.1, .1, .1];
koff = [.1, .1, .1];
[time, noCompete] = ElutionModel(100, 100, 100, 0, 0, kon, koff);
plot(time, noCompete); hold on;
[time, compete10] = ElutionModel(100, 100, 100, 50, 0, kon, koff);
plot(time, compete10); hold on;
[time, compete100] = ElutionModel(100, 100, 100, 200, 0, kon, koff);
plot(time, compete100);
title('Competing Antibody Elution Method - Increasing Concentration');
xlabel('Time (s)');
ylabel('Number of Bound Signal Antibodies');
legend('No Competing Antibody', 'Competing Ab, .5X Concentration', 'Competing Ab, 2X Concentration');

% Method - Competing Protein, Increasing Amount
figure;
kon = [.1, .1, .1];
koff = [.1, .1, .1];
[time, noCompete] = ElutionModel(100, 100, 100, 0, 0, kon, koff);
plot(time, noCompete); hold on;
[time, compete10] = ElutionModel(100, 100, 100, 0, 50, kon, koff);
plot(time, compete10); hold on;
[time, compete100] = ElutionModel(100, 100, 100, 0, 200, kon, koff);
plot(time, compete100);
title('Excess Protein Elution Method - Increasing Concentration');
xlabel('Time (s)');
ylabel('Number of Bound Signal Antibodies');
legend('No Excess Protein', 'Excess Protein, .5X Concentration', 'Excess Protein, 2X Concentration');
