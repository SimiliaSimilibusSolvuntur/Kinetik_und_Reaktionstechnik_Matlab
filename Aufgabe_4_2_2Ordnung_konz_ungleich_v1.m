% Aufgabe 4.2 Chemische Kinetik und Reaktionstechnik
% autor = Luca Moschen 
% date = 12.03.2021
% title = Reaktion 2. Ordnung

% Reaktion 2. Ordnung c_A != c_B: A + B -> Produkt

clearvars; close all; clc

c_A_0 = 1; % [mol/L] Anfangskonzentration
c_B_0 = 1.5; % [mol/L] Anfangskonzentration
k_2AB = 1e-3; % [1/s] Geschw. Konstante 

t = linspace(0, 3600, 50);  %[s] Zeit t als Zeilenvektor

% [mol/l] Formel zur Berechnung der Konzentration
c_A = (c_A_0 * (c_A_0 - c_B_0) * exp((c_A_0 - c_B_0) * k_2AB * t)) ./...
      (c_A_0 * exp((c_A_0 - c_B_0) * k_2AB * t) - c_B_0); 
c_B = c_B_0 - (c_A_0 - c_A);

% erstellen eines Graphen
figure; hold on; grid on
title('Aufgabe 4.2');
plot(t/60, c_A, t/60, c_B);
xlabel('Zeit [min]');
ylabel('Konzentration [mol/l]');
legend('$c_A$', '$c_B$', 'interpreter', 'latex');
set(gca, 'YLim', [0 2]); % y-Achse bis 2
