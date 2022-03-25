%% Matlab Konfiguration

clearvars; close all; clc;
format COMPACT

%% Physikalische Konstanten

R = 8.314;       % [J/(mol*K)] ideale Gaskonstante

%% Komponenten

% Molmassen der einzelnen Komponenten
mw_i = [200, 100, 160, 70]; % [g/mol]

%% Versuchsbedingungen gemäss Kapitel 6.5.2 (Aufgabe 6.2)

% A + B -> C + 2D 

% Anfangskonzentrationen [mol/l]
c_A_0 = 2; % [mol/l]
c_B_0 = 2; % [mol/l]
c_C_0 = 0; % [mol/l]
c_D_0 = 0; % [mol/l]

%% Reaktionsparameter

k_2 = 0.005;     % [l/(mol*s)] Geschwindigkeitskonstante bei 25 °C, 2. Ordnung
nu_i = [-1, -1, 1, 2, 0]; % [-] Vektor mit stöchiometrische Koeffizienten der Komponenten
n = 1;       % [-] Teilreaktionsordnung der Komponente "A"
m = 1;       % [-] Teilreaktionsordnung der Komponente "B"

a = 1;
b = 1;
c = 1;
d = 2;

%% Anforderungen an die Reaktion und den Reaktor

% Anforderungen an die Anlage:

Cap_CR = 100; % [t/a] Produktionskapazität
Verf = 8000;    % [h/a] Verfügbarkeit der Anlage

% Gewünschter Umsatz bezüglich der Komponente A am Ausgang des Reaktors
X_A_f = 0.8; % [-]

%% Feed-Konzentrationen

c_A_in = c_A_0 / 2.5; % [mol/l]
c_B_in = c_B_0 / 2.5 * 1.5; % [mol/l]

%% Aufgabe 7.17, Lösungsvariante 1: Berechnung der mittleren Verweilzeit via symbolische Integration

% Erforderliche mittlere Verweil-Zeit nach Gl. (7.63):
c_A_f = c_A_in * (1 - X_A_f);    % [mol/l] erforderliche Endkonzentration
% tau_1 = ((1 / c_A_f) - (1 / c_A_in)) / k_2; %[s] mittlere Verweilzeit
tau_1 = log(c_B_in * c_A_f / (c_A_in * (c_B_in - c_A_in + c_A_f)))/...
    (k_2 * (c_A_in - c_B_in)); %[s] mittlere Verweilzeit
disp(['Lösungsvariante 1: tau = ', num2str(tau_1/60, '%.4g'), ' min']);

% Minimales Reaktorvolumen "V_R" (7.35) auflösen nach V_R -> (7.48)
mw_Prod = mw_i(4); % Molmasse des Produktes, Komponente "C"
a = 1; % [-] Stöchiometrischer Koeffizient der limitierenden Komponente (A)
d = 2; % [-] Stöchiometrischer Koeffizient des Produktes (D)
V_R = Cap_CR * 1000 * tau_1 / 3600 / (c_A_in * mw_Prod * Verf * X_A_f) * 1000 * a/d; %[l]
disp(['V_R = ', num2str(V_R/1000, '%.2g'), ' m3']);

% Volumen-Zeit-Ausbeute gemäss Gl. (7.14):
VT_Yield = Cap_CR * 1000 / V_R; % [kg/(l*a)]
disp(['VT_Yield = ', num2str(VT_Yield, '%.2g'), ' kg/(l*a)']);

% Volumenstrom in den Reaktor gemäss (7.32)
Vf_in = V_R/tau_1;   % [l/s]
disp(['Vf_in = ', num2str(Vf_in, '%.2g'), ' l/s']);

%% Aufgabe 7.17, Lösungsvariante 2: Berechnung der mittleren Verweilzeit via einfache grafische Integration

n_X = 100;
X_A = linspace(0, X_A_f, n_X);
c_A = c_A_in * (1 - X_A);       % [mol/l] Konzentration "A" als Funktion von X_A
c_B = c_B_in - c_A_in * X_A;  % [mol/l] Konzentration "B" als Funktion von X_A
r = k_2 * c_A .* c_B;         % [mol/(l*s)] Reaktionsgeschwindigkeit ()
r_A = -a * r;                   % [mol/(l*s)] Umsatzgeschwindigkeit bezüglich "A" (7.1)

inv_r_A = 1./r_A;           % [(l*s)/mol] Inverse Umsatzgewschindigkeit

%% Berechnen der Fläche unter der Kurve -1/r_A vs. X_A anhand der Trapez-Methode

% Flaeche_1 = sum(X_A_f/(n_X-1)*-mean([inv_r_A(1:end-1); inv_r_A(2:end)]));
Flaeche_1 = trapz(X_A, -inv_r_A)

% Berechnung der mittleren Verweilzeit gemäss (7.61)
tau_2a= c_A_in * Flaeche_1; %[s]
disp(['Lösungsvariante 2 mit Gl. (7.61): tau = ', num2str(tau_2a/60, '%.4g'), ' min']);

% Minimales Reaktorvolumen "V_R" (7.35) auflösen nach V_R -> (7.48)
mw_Prod = mw_i(4); % Molmasse des Produktes, Komponente "C"
a = 1; % [-] Stöchiometrischer Koeffizient der limitierenden Komponente (A)
d = 2; % [-] Stöchiometrischer Koeffizient des Produktes (D)
V_R_2a = Cap_CR * 1000 * tau_2a / 3600 / (c_A_in * mw_Prod * Verf * X_A_f) * 1000 * a/d; %[l]
disp(['V_R_2a = ', num2str(V_R_2a/1000, '%.2g'), ' m3']);

% Volumen-Zeit-Ausbeute gemäss Gl. (7.14):
VT_Yield = Cap_CR * 1000 / V_R_2a; % [kg/(l*a)]
disp(['VT_Yield_2a = ', num2str(VT_Yield, '%.2g'), ' kg/(l*a)']);

% Volumenstrom in den Reaktor gemäss (7.32)
Vf_in = V_R_2a/tau_2a;   % [l/s]
disp(['Vf_in_2a = ', num2str(Vf_in, '%.2g'), ' l/s']);

%% Berechnen der Fläche unter der Kurve -1/r_A vs. c_A anhand der Trapez-Methode

% Flaeche_2 = sum((c_A_in-c_A_f)/(n_X-1)*-mean(-[inv_r_A(1:end-1); inv_r_A(2:end)]))
Flaeche_2 = trapz(c_A, -inv_r_A)

% Berechnung der mittleren Verweilzeit gemäss (7.62)
tau_2b = -Flaeche_2; %[s]
disp(['Lösungsvariante 2 mit Gl. (7.62): tau  = ', num2str(tau_2b/60, '%.4g'), ' min']);

% Minimales Reaktorvolumen "V_R" (7.35) auflösen nach V_R -> (7.48)
mw_Prod = mw_i(4); % Molmasse des Produktes, Komponente "C"
a = 1; % [-] Stöchiometrischer Koeffizient der limitierenden Komponente (A)
d = 2; % [-] Stöchiometrischer Koeffizient des Produktes (D)
V_R_2b = Cap_CR * 1000 * tau_2b / 3600 / (c_A_in * mw_Prod * Verf * X_A_f) * 1000 * a/d; %[l]
disp('')
disp(['V_R_2b = ', num2str(V_R_2b/1000, '%.2g'), ' m3']);

% Volumen-Zeit-Ausbeute gemäss Gl. (7.14):
VT_Yield_2b = Cap_CR * 1000 / V_R_2b; % [kg/(l*a)]
disp(['VT_Yield = ', num2str(VT_Yield, '%.2g'), ' kg/(l*a)']);

% Volumenstrom in den Reaktor gemäss (7.32)
Vf_in_2b = V_R_2b/tau_2b;   % [l/s]
disp(['Vf_in = ', num2str(Vf_in, '%.2g'), ' l/s']);

%% Graphen

figure;
hold on; grid on;
plot(X_A, -inv_r_A, 'o-');
for i = 1:n_X-1
    X = [X_A(i), X_A(i+1), X_A(i+1), X_A(i)];
    Y = [0, 0, -inv_r_A(i+1), -inv_r_A(i)];
    patch('XData', X,'YData',Y, 'FaceColor','blue','FaceAlpha',.1);
end
xlabel('$X_A$ [mol/mol]', 'interpreter', 'latex');
ylabel('$-\frac{1}{r_A}$ [(l$\cdot$s)/mol]','interpreter', 'latex');
legend('$-\frac{1}{r_A(X_A)}$', ...
    ['Flaeche = ' num2str(Flaeche_1, '%.2g'), ' [(l$\cdot$s)/mol]', newline , ...
    '$\tau$ = $c_{A,in}\cdot$Flaeche = ', num2str(tau_2a/60, '%.4g'), ' min'], 'interpreter', 'latex');
title('Gl. (7.61) - analog zu Gl. (7.21)', 'interpreter', 'latex');
set(gca, 'XLim', [0 1]);

figure; hold on; grid on;
plot(c_A, -inv_r_A, 'o-');
for i = 1:n_X-1
    X = [c_A(i), c_A(i+1), c_A(i+1), c_A(i)];
    Y = [0, 0, -inv_r_A(i+1), -inv_r_A(i)];
    patch('XData', X,'YData',Y, 'FaceColor','blue','FaceAlpha',.1);
end
xlabel('$c_A$ [mol/l]', 'interpreter', 'latex');
ylabel('$-\frac{1}{r_A}$ [(l$\cdot$s)/mol]','interpreter', 'latex');
legend('$-\frac{1}{r_A(c_A)}$', ...
    ['Flaeche = ' num2str(Flaeche_2, '%.2g'), ' [s]', newline , ...
    '$\tau$ = -Flaeche = ', num2str(tau_2b/60, '%.4g'), ' min'], 'interpreter', 'latex');
title('Gl. (7.62) - analog zu Gl. (7.20)', 'interpreter', 'latex');
XLim = get(gca, 'XLim');
set(gca, 'XLim', [0 XLim(2)]);
