% sars valency = 60
% if we do excessprotein method, original proteins are high valency, excess
% protein is 2

% % Method - Competing Antibody, Increasing Amount
kon = 1; koff = 0.01;
[time1, noCompete] = ElutionModel(70, 70, 70, 0, 0, 60, kon, koff);
[time2, compete10] = ElutionModel(70, 70, 70, 70, 0, 60, kon, koff);
[time3, compete100] = ElutionModel(70, 70, 70, 7000, 0, 60, kon, koff);

figure;
plot(time1, noCompete); hold on;
plot(time2, compete10); hold on;
plot(time3, compete100);
title('Competing Antibody Elution Method - Valency = 60');
xlabel('Time (s)');
ylabel('Number of Bound Signal Antibodies');
legend('No Competing Antibody', 'Competing Ab, 1X Concentration', 'Competing Ab, 100X Concentration');
ylim([0, 100]);

% Method - Competing Protein, Increasing Amount
figure;
kon = 1; koff = 0.01;
[time1, noCompete] = ElutionModel(70, 70, 70, 0, 0, 60, kon, koff);
[time2, compete10] = ElutionModel(70, 70, 70, 0, 60, 60, kon, koff);
[time3, compete100] = ElutionModel(70, 70, 70, 0, 120, 60, kon, koff);
plot(time1, noCompete); hold on;
plot(time2, compete10); hold on;
plot(time3, compete100);
title('Excess Protein Elution Method - Valency = 60');
xlabel('Time (s)');
ylabel('Number of Bound Signal Antibodies');
legend('No Excess Protein', 'Excess Protein, 1X Concentration', 'Excess Protein, 2X Concentration');