clearvars; close all; clc
% autor = Luca Moschen 
% date = 29.03.2021
%% Aufgabe 6.1
% Umsetzung von Ethylacetat (A) mit NaOH (B)
c_A_0 = 0.00486; % [mol/L]
c_B_0 = 0.0098; % [mol/L]
c_C_0 = 0; % [mol/L]
c_D_0 = 0; % [mol/L]
k_2AB = 0.11; %[L/(mol * s)]

t = linspace(0, 3600, 100);  %[s] Zeit t als Zeilenvektor

% A + B -> C + D
c_A = (c_A_0 * (c_A_0 - c_B_0) * exp((c_A_0 - c_B_0) * k_2AB * t)) ./...
      (c_A_0 * exp((c_A_0 - c_B_0) * k_2AB * t) - c_B_0); 

c_B = c_B_0 - (c_A_0 - c_A);

c_C = 0.0098 - c_B;
c_D = 0.0098 - c_B;

X_A = 1 - (c_A ./ c_A_0);

% Plot
figure; hold on; grid on
title('Aufgabe 6.1');
xlabel('Zeit [s]');
set(gca, 'XLim', [0 4000]);

yyaxis right;
plot(t, X_A);
ylabel('Umsetzung X_A [-]');
set(gca, 'YLim', [0 1]);

yyaxis left;
plot(t, c_A, t, c_B, t, c_C, t, c_D);
ylabel('Konzentration [mol/l]');
set(gca, 'YLim', [0 0.01]);

legend('$c_A$', '$c_B$', '$c_C$', '$c_D$', '$X_A$', 'interpreter', 'latex');