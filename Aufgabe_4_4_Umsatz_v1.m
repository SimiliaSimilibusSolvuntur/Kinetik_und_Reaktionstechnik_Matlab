% Aufgabe 4.4 Chemische Kinetik und Reaktionstechnik
% autor = Luca Moschen 
% date = 12.03.2021

clearvars; close all; clc

t = linspace(0, 60*120, 50);  %[s] Zeit t als Zeilenvektor
c_A_0 = 1; % [mol/L] Anfangskonzentration
c_B_0 = 1.5; % [mol/L] Anfangskonzentration

% Reaktion 2. Ordnung Sonderfall c_A = c_B: A + B -> Produkt
k_1 = 1e-3; % [1/s] Geschw. Konstante 
% GL 4.10
c_A_1 = c_A_0 * exp(-k_1 * t); % [mol/l]

% Reaktion 2. Ordnung Sonderfall c_A = c_B: A + B -> Produkt
k_2 = 1e-3; % [1/s] Geschw. Konstante 
% GL 4.10
c_A_2 = 1 ./ (2 * k_2 * t + (1 ./ c_A_0)); % [mol/l]

% Reaktion 2. Ordnung c_A != c_B: A + B(Ueberschuss) -> Produkt
k_2AB_1 = 1e-3; % [l/(mol*s)] Geschw. Konstante 
% GL 4.20
c_A_3 = (c_A_0 * (c_A_0 - c_B_0) * exp((c_A_0 - c_B_0) * k_2AB_1 * t)) ./...
      (c_A_0 * exp((c_A_0 - c_B_0) * k_2AB_1 * t) - c_B_0); 

% erstellen eines Graphen
figure; hold on; grid on
title('Aufgabe 4.1');
plot(t/60, (c_A_0 - c_A_1) / c_A_0, 'b', ...
    t/60, (c_A_0 - c_A_2) / c_A_0, 'r', ...
    t/60, (c_A_0 - c_A_3) / c_A_0, 'k');
xlabel('Zeit [min]');
ylabel('Konzentration [mol/l]');
legend('$c_A$ 1. Ordnung', '$c_A$ 2. Ordnung $c_A$ = $c_B$',...
    '$c_A$ 2. Ordnung $c_A$ != $c_B$', 'interpreter', 'latex');
%set(gca, 'YLim', [0 2]);


