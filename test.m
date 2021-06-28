% testing
clear;clc;close all;


method = Methods.ExcessProtein;

figure;
kon = [.1, .1];
koff = [.1, .1];
[time, noExcess] = ElutionModel(method, 1000, 1000, 0, kon, koff);
plot(time, noExcess); hold on;
[time, excess1] = ElutionModel(method, 1000, 1000, 1000, kon, koff);
plot(time, excess1); hold on;
[time, excess5] = ElutionModel(method, 1000, 1000, 5000, kon, koff);
plot(time, excess5);
title('Excess Protein Elution Method');
xlabel('Time (s)');
ylabel('Number of Bound Signal Antibodies');
legend('No Excess', '1X Excess Secto', '5X Excess Secto');


% 
% method = Methods.CompetingAntibody;
% 
% figure;
% kon = [.1, .1];
% koff = [.1, .1];
% [time, noCompete] = ElutionModel(method, 100, 100, 0, kon, koff);
% plot(time, noCompete); hold on;
% kon = [.1, 1];
% koff = [.1, .01];
% [time, compete10] = ElutionModel(method, 100, 100, 100, kon, koff);
% plot(time, compete10); hold on;
% kon = [.1, 10];
% koff = [.1, .001];
% [time, compete100] = ElutionModel(method, 100, 100, 100, kon, koff);
% plot(time, compete100);
% title('Competing Antibody Elution Method');
% xlabel('Time (s)');
% ylabel('Number of Bound Signal Antibodies');
% legend('No Competing Antibody', 'Competing Ab, 10X Kinetics', 'Competing Ab, 100X Kinetics');
