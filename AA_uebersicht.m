% Übersicht Codeblöcke

%% Geschwindigkeitsgesetze
% Erste Ordnung; A -> Produkte
r = k_1 * c_A; % k_1 = [1/s]; r =[mol/(l*s)]
r_A = -r;
c_A = c_A_0 * exp(-k_1 .* t); % symbolische Intergration
log(c_A) = log(c_A_0) - k_1 .* t; % Linearisierung

% Zweite Ordnung; 2A -> Produkte
r = k_2 * (c_A)^2; % k_2 =[(l/mol)*1/s]
r_A = -2 * r;
c_A = 1/(2 * k_2 * t + (1/c_A_0)); % symbolische Integration
1/c_A = 2 * k_2 * t + 1/c_A_0; % Linearisierung

% Zweite Ordnung Sonderfall c_A = c_B; A + B -> Produkte
% Bedingung: Die Anfangskonzentrationen von A und B müssen gleich sein
r = k_2 * c_AB; % k_2 =[(l/mol)*1/s]
r_A = -r;
c_A = 1/(k_2 * t * + 1/c_AB_0); % symbolische Integration
1/c_AB = k_2 * t + 1/c_AB_0; % Linearisierung

% Zweite Ordnung mit c_A != c_B; A + B -> Produkte
r = k_2AB * c_A * c_B; % k_2AB = [(l/mol)*1/s]
% wenn a = b = 1:
c_B = c_B_0 - (c_A_0 - c_A);
r_A = -r;
c_A = (c_A_0 * (c_A_0 - c_B_0) * exp((c_A_0 - c_B_0) * k_2AB * t)) / ...
    (c_A_0 * exp((c_A_0 - c_B_0) * k_2AB * t)) - c_B_0); % symbolische Integration
log((c_B_0 * c_A)/(c_A_0 * (c_B_0 - c_A_0 + c_A)) = k_2AB * (c_A_0 -c_B_0) * t; % Linearisierung

% Grafisch
figure; grid on
title('XXXXXXXXXXX');
plot(t, c_A, 'x'); % Zeit richtig umformen (s)
hold on;
plot(t, c_B, 'o');
xlabel('Zeit [min]');
ylabel('Konzentration [mol/l]');
legend('$c_A$', '$c_B$', 'interpreter', 'latex');
% set(gca, 'YLim', [0 2]); % Optional

%% Geschwindigkeitsgesetze als Funktion des Umsatzes einer Komponente
%Umsatz generell
X_A = (n_A_0 - n_A) / n_A_0;
X_A = (c_A_0 - c_A) / c_A_0; % solange das Volumen der Reaktionsmasse konstant bleibt

% Konzentrationen der Komponenten als Funktion von X_A bei einer Reaktion nach a*A + b*B -> c*C + d*D
% Bedingung: Komponente A ist stöchiometrisch limitierend
c_A = c_A_0 * (1 - X_A);
c_B = c_B_0 - b/a * c_A_0 * X_A;
c_C = c_C_0 + c/a * c_A_0 * X_A;
c_D = c_D_0 + d/a * c_A_0 * X_A;

% Grafisch; Geschwindigkeitsgesetze mit Umsatz kombiniert
figure, grid on, hold on
title('XXXXXXXXX');
xlabel('Zeit [s]');
% set(gca, 'XLim', [0 900]); % Limitierung der Achse

yyaxis left;
ylabel('$c_i$ [mol/l]', 'interpreter', 'latex');
plot(t, c_A); % Zeit richtig umformen (s)
plot(t, c_B); % Zeit richtig umformen (s)
plot(t, c_D); % Zeit richtig umformen (s)
% set(gca, 'YLim', [0 1]); % Limitierung der Achse

yyaxis right
ylabel('Umsatz $X_A$', 'interpreter', 'latex');
plot(t, X_A); % Zeit richtig umformen (s)
% set(gca, 'YLim', [0 1]); % Limitierung der Achse

legend('$c_A$','$c_B$', '$c_D$', '$X_A$', 'interpreter', 'latex');










