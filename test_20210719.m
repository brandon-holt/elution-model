clear; clc; close all;

% first test effect of separation betwee kon/koff
figure;
starting_power = 6;
orders_of_magnitude = 12;
repeats = 5;
kon = 1;
koff = logspace(starting_power, starting_power-orders_of_magnitude, orders_of_magnitude+1);
data = zeros(orders_of_magnitude+1, 2); % col-1 = kon/koff, % col-2 = #bound at endpoint
for i = 1:numel(koff)
    for j = 2:repeats+1
        [t, b] = ElutionModel(100, 100, 100, 100, 0, 2, kon, koff(i));
        data(i, 1) = kon / koff(i);
        data(i, j) = (b(1) - b(end)) / t(end);
    end
end

loglog(data(:,1), mean(data(:,2:end),2), 'Marker','o');
xlabel('Kon/Koff ratio');
ylabel('Protein Elution Velocity (s^{-1})');
title('Protein Elution at Varying Kd');

% test the effect of valency on the same graph
figure;
[t, b] = ElutionModel(1000, 100, 100, 1000, 0, 2, 1, .001);
plot(t, b); hold on;
[t, b] = ElutionModel(1000, 100, 100, 1000, 0, 60, 1, .001);
plot(t, b); hold on;
xlabel('Time (seconds)'); ylabel('Number of Proteins Bound');
title('Effect of Valency on Elution');
legend('Valency = 2', 'Valency = 60');