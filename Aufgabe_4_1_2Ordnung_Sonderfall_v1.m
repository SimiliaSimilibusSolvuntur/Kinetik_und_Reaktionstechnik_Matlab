% Aufgabe 4.1 Chemische Kinetik und Reaktionstechnik
% autor = Luca Moschen 
% date = 12.03.2021
% title = Reaktion 2. Ordnung Sonderfall

% Reaktion 2. Ordnung Sonderfall c_A = c_B: A + B -> Produkt

clearvars; close all; clc

c_A_0 = 2; % [mol/L] Anfangskonzentration
c_B_0 = 2; % [mol/L] Anfangskonzentration entspricht c_A_0
k_2 = 1e-4; % [1/s] Geschw. Konstante 

t = linspace(0, 60*120, 50);  %[s] Zeit t als Zeilenvektor

c_A = 1 ./ (k_2 * t + (1 ./ c_A_0)); % [mol/l] Formel zur Berechnung der Konzentration

% erstellen eines Graphen
figure; hold on; grid on
title('Aufgabe 4.1');
plot(t/60, c_A);
xlabel('Zeit [min]');
ylabel('Konzentration [mol/l]');
legend('$c_A$', 'interpreter', 'latex');
set(gca, 'YLim', [0 2]);

% erstellen eines Graphen (Linearisierung)
figure; hold on; grid on
title('Aufgabe 4.1 linearisiert 1 / c_A'); % Linearisierung 1 / c_A
plot(t/60, 1./c_A);
xlabel('Zeit [min]');
ylabel('[l/mol]');
legend('$c_A$', 'interpreter', 'latex');
