%% Matlab Konfiguration

clearvars; close all; clc;
format COMPACT

%% Physikalische Konstanten

R = 8.314;       % [J/(mol*K)] ideale Gaskonstante

%% Versuchsbedingungen 

% Anfangskonzentrationen [mol/l]
c_A_0 = 2.8; % [mol/l] 
c_B_0 = 1.6; % [mol/l] 
c_C_0 = 0; % [mol/l] 

%% Reaktionsparameter

k_1 = 12.5 / 60;     % [L^2/(mol^2*s)] Geschwindigkeitskonstante 
k_2 = 1.5 / 60;     % [1/s] Geschwindigkeitskonstante 

% A + 2B -> C
a = 1;
b = 2;
c = 1;

nu_i = [-1, -2, 1, 0, 0]; % [-] Vektor mit stöchiometrische Koeffizienten der Komponenten
n = 1;       % [-] Teilreaktionsordnung der Komponente "A"
m = 1;       % [-] Teilreaktionsordnung der Komponente "B"
o = 1;

%% Anforderungen an die Reaktion und den Reaktor

V_R = 6; % [l]

% Gewünschter Umsatz bezüglich der Komponente B am Ausgang des Reaktors
X_B_f = 0.75; % [-]

%% Feed-Konzentrationen

c_A_in = c_A_0 / 2; % [mol/l]
c_B_in = c_B_0 / 2; % [mol/l]

%% Aufgabe 7.16, Lösungsvariante 2: Berechnung der mittleren Verweilzeit via einfache grafische Integration

n_X = 100;
X_B = linspace(0, X_B_f, n_X);
c_B = c_B_in * (1 - X_B);  % [mol/l] Konzentration "B" als Funktion von X_B
c_A = c_A_in - (a / b) * c_B_in * X_B; % [mol/l] Konzentration "A" als Funktion von X_B
c_C = c_C_0 + (c / b) * c_B_in * X_B;    % [mol/l] Komponente C

r_B = -2 * (k_1 .* c_A .* c_B.^2 - k_2 .* c_C); % [mol/(l*s)] Umsatzgeschwindigkeit 

inv_r_A = 1 ./ r_B; % [(l*s)/mol] Inverse Umsatzgewschindigkeit

% Berechnen der Fläche unter der Kurve -1/r_B vs. X_B anhand der Trapez-Methode
% Flaeche_1 = sum(X_A_f/(n_X-1)*-mean([inv_r_A(1:end-1); inv_r_A(2:end)]));
Flaeche_1 = trapz(X_B, -inv_r_A)

% Berechnung der mittleren Verweilzeit gemäss (7.61)
tau_2a = c_B_in * Flaeche_1; %[s]
disp(['Lösungsvariante 2 mit Gl. (7.61): tau = ', num2str(tau_2a/60, '%.4g'), ' min']);

% Berechnen der Fläche unter der Kurve -1/r_A vs. c_A anhand der Trapez-Methode
% Flaeche_2 = sum((c_A_in-c_A_f)/(n_X-1)*-mean(-[inv_r_A(1:end-1); inv_r_A(2:end)]))
Flaeche_2 = trapz(c_B, -inv_r_A)

% Berechnung der mittleren Verweilzeit gemäss (7.62)
tau_2b = -Flaeche_2; %[s]
disp(['Lösungsvariante 2 mit Gl. (7.62): tau  = ', num2str(tau_2b/60, '%.4g'), ' min']);

% Volumenströme der beiden Feeds () V_R = V_RM
Vf_in = V_R / tau_2b / 2; % [l/s]
disp(['Vf_in = ', num2str(Vf_in, '%.2g'), ' l/s']);

figure;
hold on; grid on;
plot(X_B, -inv_r_A, 'o-');
for i = 1:n_X-1
    X = [X_B(i), X_B(i+1), X_B(i+1), X_B(i)];
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
plot(c_B, -inv_r_A, 'o-');
for i = 1:n_X-1
    X = [c_B(i), c_B(i+1), c_B(i+1), c_B(i)];
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

%% Aufgabe 7.15, Lösungsvariante 3: Berechnung der mittleren Verweilzeit mittels ODE-Solver. Basis: Gl. (7.62)

% Stöchiometrische Koeffizienten
% Hinreaktion
nu_i_1 = [-1, -2, 1]; % [-] Stöchiometrische Koeffizienten der Komponenten

% Rückreaktion
nu_i_2 = [1, 2, -1]; % [-] Stöchiometrische Koeffizienten der Komponenten

% Reaktor-Geometrie (Rohrreaktor)
d_R = 1; % [m]    Durchmesser des Rohrreaktors

% Modell-Parameter zur Beschreibung des kontinuierlichen Feeds
% Dichte des Feeds in den Reaktor
rho_in = 1000;          % [kg/m3]
% Feed-Konzentrationen
% c_Lsm_i = (rho_in-c_A_in*mw_i(1)-c_B_in*mw_i(2))/mw_i(5);
c_i_in = [c_A_in, c_B_in, 0]; % [mol/l] Feed-Konzentrationen

% Volumenstrom in den Reaktor gemäss (7.32)
Vf_in = 1/60;   % [l/s] Zahl hat keinen Einfluss auf die Berechnung von tau
disp(['Vf_in = ', num2str(Vf_in, '%.2g'), ' l/s']);

% Stoffströme in den Reaktor
F_i_in = c_i_in.*Vf_in; %[mol/s]

% Definieren der strukturierten Varialbe "MP" mit den Modell-Parametern für des Reaktor-Modell
% MP.mw_i = mw_i;         % [g/mol] Vektor mit Molmassen der einzelnen Komponenten
MP.k_1 = k_1;           % [...] Geschwindigkeitskonstante
MP.k_2 = k_2;
MP.nu_i_1 = nu_i_1;         % [-] Vektor mit stöchiometrischen Koeffizienten der Komponenten
MP.nu_i_2 = nu_i_2;
MP.n = n;               % [-] Teilreaktionsordnung der Komponente "A"
MP.m = m;               % [-] Teilreaktionsordnung der Komponente "B"
MP.o = o;
MP.Vf_in = Vf_in;       % [l/s] Volumenstrom in den Reaktor in den Reaktor
MP.rho_in = rho_in;     % [kg/m3] Dichte des Feeds in den Reaktor
MP.d_R  = d_R;          % [m] Durchmesser des Reaktorrohres

% Anfangswerte der ortsabhängigen Variablen im Zeilenvektor "y"
y0 = [F_i_in];

% Bedingungen für die Simulation
n_l = 1000;       % [-] Anzahl der Punkte für die Resultat-Matrix
tau_max = 2 * 3600;                    % [s] Maximal simulierte Verweilzeit
V_R_max = tau_max*Vf_in;             % [l] Maximal simuliertes Reaktorvolumen
l_max = V_R_max / 1000 / ((d_R / 2)^2 * pi); % [m] Maximal simulierte Rohrlänge
lspan = linspace(0, l_max, n_l); % [s] Vektor mit Ortspunkte für die Resultat-Matrix

% Genauigkeit der numerischen Integration
options = odeset('AbsTol', 1e-10);  % Absolutes Abbruchkriterium
options = odeset(options, 'RelTol', 1e-6); % Relatives Abbruchkriterium

% Aufruf des ode-solvers
[Sim.l, Sim.y] = ode23s(@RMS_PFR_AnBm_it_7_16, lspan, y0, options, MP);

% Resultate der numerischen Integration
Sim.F_l_i = Sim.y;                  % [mol] Stoffströme
Sim.c_l_i = Sim.F_l_i ./ Vf_in;     % [mol/l] Konzentrations-Ort-Kurven

% Umsatz der söchiometrisch limiteirenden Komponente A gemäss Gl.(7.49)
i = 1;
Sim.X_A_l = (F_i_in(i) - Sim.F_l_i(:, i)) ./ F_i_in(i); % [mol/mol]
Sim.tau = (Sim.l*(d_R/2)^2*pi)*1000/Vf_in; % [s]

% Grafische Darstellung der Resultate
figure
grid on; grid minor; hold on
plot(Sim.l, Sim.c_l_i(:, 1), 'k-', 'linewidth', 2);
plot(Sim.l, Sim.c_l_i(:, 2), 'k', 'linewidth', 0.5);
plot(Sim.l, Sim.c_l_i(:, 3), 'k:', 'linewidth', 1);
plot(Sim.l, Sim.c_l_i(:, 4), 'k-.', 'linewidth', 1);
xlabel('Rohrlaenge $l$ [m]', 'interpreter', 'latex');
ylabel('Konzentrationen $c_{i}$ [mol/l]', 'interpreter', 'latex');
title(['Reaktion im PFR mit $d_R$ = ', num2str(d_R, '%.2g'), ...
    ' m und $\dot{V}_{in}$ = ', num2str(Vf_in*60, '%.2g'), ' l/min'],'interpreter', 'latex');
yyaxis right
plot(Sim.l, Sim.X_A_l, 'linewidth', 2);
ylabel('$X_A$ [-]', 'interpreter', 'latex');
legend('c_A', 'c_B', 'c_C', 'c_D', 'X_A');
set(gca, 'YLim', [0 1]);

figure
grid on; grid minor; hold on
plot(Sim.tau/3600, Sim.c_l_i(:, 1), 'k-', 'linewidth', 2);
plot(Sim.tau/3600, Sim.c_l_i(:, 2), 'k', 'linewidth', 0.5);
plot(Sim.tau/3600, Sim.c_l_i(:, 3), 'k:', 'linewidth', 1);
plot(Sim.tau/3600, Sim.c_l_i(:, 4), 'k-.', 'linewidth', 1);
xlabel('Verweilzeit $\tau$ [h]', 'interpreter', 'latex');
ylabel('Konzentrationen $c_{i}$ [mol/l]', 'interpreter', 'latex');
title(['Reaktion im PFR bei 25 °C']);
yyaxis right
plot(Sim.tau/3600, Sim.X_A_l, 'linewidth', 2);
ylabel('$X_A$ [-]', 'interpreter', 'latex');
legend('c_A', 'c_B', 'c_C', 'c_D', 'X_A');
set(gca, 'YLim', [0 1]);
set(gca ,'XLim', [0 tau_max]/3600);

% Besimmung der erforderlichen Verweilzeit für den gewünschten Umsatz X_A_f
% Index bei welchem Sim.X_A_l > X_A_f erfüllt ist
tmp = find(Sim.X_A_l >= X_B_f);
i_X_Af = tmp(1);

tau_3 = Sim.tau(i_X_Af); % [s]
disp(['Lösungsvariante 3 mit ode-Solver: tau = ', num2str(tau_3/3600, '%.2g'), ' h']);
disp(['Lösungsvariante 3 mit ode-Solver: tau = ', num2str(tau_3/60, '%.4g'), ' min']);
disp(['Lösungsvariante 3 mit ode-Solver: tau = ', num2str(tau_3, '%.4g'), ' sec']);


