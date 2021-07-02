% sars valency = 60
% if we do excessprotein method, original proteins are high valency, excess
% protein is 2

% % Method - Competing Antibody, Increasing Amount
[~, noCompete] = ElutionModel(60, 60, 60, 0, 0, .1, .1);
[~, compete10] = ElutionModel(60, 60, 60, 60, 0, .1, .1);
[time, compete100] = ElutionModel(60, 60, 60, 6000, 0, .1, .1);

figure;
plot(time, noCompete); hold on;
plot(time, compete10); hold on;
plot(time, compete100);
title('Competing Antibody Elution Method - Valency = 60');
xlabel('Time (s)');
ylabel('Number of Bound Signal Antibodies');
legend('No Competing Antibody', 'Competing Ab, 1X Concentration', 'Competing Ab, 100X Concentration');
ylim([0, 100]);

% Method - Competing Protein, Increasing Amount
figure;
kon = [.1, .1, .1];
koff = [.1, .1, .1];
[time, noCompete] = ElutionModel(60, 60, 60, 0, 0, kon, koff);
plot(time, noCompete); hold on;
[time, compete10] = ElutionModel(60, 60, 60, 0, 60, kon, koff);
plot(time, compete10); hold on;
[time, compete100] = ElutionModel(60, 60, 60, 0, 120, kon, koff);
plot(time, compete100);
title('Excess Protein Elution Method - Valency = 60');
xlabel('Time (s)');
ylabel('Number of Bound Signal Antibodies');
legend('No Excess Protein', 'Excess Protein, 1X Concentration', 'Excess Protein, 2X Concentration');