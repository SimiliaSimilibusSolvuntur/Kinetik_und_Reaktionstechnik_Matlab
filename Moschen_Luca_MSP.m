%% Matlab Konfiguration
clearvars; close all; clc;      
format COMPACT

%% a.)
% 1 SA + 1 AcOAc -> 1 ASS + 1 AcOH
mw_i = [138.12, 102.09, 180.16, 60.05, 41.05]; % [g/mol] 

V_m = 20 + 0.2 + 3; % [ml]

% Dichte
rho_RM = 780;   % [kg/m3] 

%% b.)
n_SA = 3.5/mw_i(1) % mol = 0.0253mol
n_A = n_SA;
n_AcOc = 3*1.08/mw_i(2) % mol = 0.0317mol
n_B = n_AcOc;
% stöch limitierend ist SA

%% c.)
% Die Graphik zeigt AcOAc dieses wird nicht komplett aufgebracht also
% wahrscheinlich zweiter ordnung da zwei Reatktanden und c_A != c_B
%% d.)

%% Physikalische Konstanten
R = 8.314;       % [J/(mol*K)] ideale Gaskonstante

%% Versuchsbedingungen gemäss Kapitel 6.5.2 (Aufgabe 6.2)
% Volumen der Reaktionsmasse [l]
V_RM_0 = V_m*1000; % [l]

% Anfangskonzentrationen [mol/l]
c_A_0 = n_A/V_RM_0; % [mol/l] SA
c_B_0 = n_B; % [mol/l] rein AcOAc
c_C_0 = 0; % [mol/l] NaOAc
c_D_0 = 0; % [mol/l] EtOH

%% Reaktions-Parameter
%% Durch ausprobieren k_2AB herausfinden
k_2AB = 1;     % [l/(mol*s)] Geschwindigkeitskonstante bei 25 °C, 2. Ordnung
nu_i = [-1, -1, 1, 1, 0]; % [-] Vektor mit stöchiometrische Koeffizienten der Komponenten 
n = 1;       % [-] Teilreaktionsordnung der Komponente "A"
m = 1;       % [-] Teilreaktionsordnung der Komponente "B"

%% Anfangsbedingungen zum Zeitpunkt t = 0 im Reaktor (=vorgelegt)
V_RM_0 = 1;                         % [l] Reaktionsvolumen 
m_RM_0 = V_RM_0*rho_RM/1000;     % [kg] Reaktionsmasse 
% n_B_0 = c_B_0*V_RM_0;               % [mol] Stoffmenge Komponente B 
n_Lsm_0 = 20*0.78/mw_i(5); % [mol] Stoffmenge Lsm
n_i_0 = [n_A, 0, 0, 0, n_Lsm_0];  % [mol] Stoffmenge Komponente i
 
%% Reaktionsparameter zur Beschreibung der Dosierung
t_dos_start = 0; % [s] Startzeit der Dosierung
t_dos = 10*60;   % [s] Dauer der Dosierung
n_dos_B = n_B; % [mol] Stoffmenge Komponente B dosiert
n_dos_i = [0, n_dos_B, 0, 0, 0]; % [mol] Stoffmenge Komponente i dosiert
F_dos_i = n_dos_i / t_dos;  %  [mol/s] Stoffstrom während der Dosierung 

%% Definieren der strukturierten Varialbe "MP" mit den Modell-Parametern für des Reaktor-Modell
MP.mw_i = mw_i;         % [g/mol] Vektor mit Molmassen der einzelnen Komponenten
MP.k = k_2AB;           % [...] Geschwindigkeitskonstante 
MP.nu_i = nu_i;         % [-] Vektor mit stöchiometrischen Koeffizienten der Komponenten
MP.n = n;               % [-] Teilreaktionsordnung der Komponente "A"
MP.m = m;               % [-] Teilreaktionsordnung der Komponente "B"
MP.rho_RM = rho_RM;     % [kg/m3] Dichte der Reaktionsmasse 
MP.t_dos_start = t_dos_start; % [s] Startzeit der Dosierung
MP.t_dos = t_dos;       % [s] Dauer der Dosierung
MP.F_dos_i = F_dos_i;   % [mol/s] Vektor mit Stoffströmen während der Dosierung 

%% Anfangswerte der zeitabhängigen Grössen im Zeilenvektor "y" 
y0 = [n_i_0, m_RM_0];

%% Bedingungen für die Simulation
t_start = 0;      % [s] Start-Zeit der Simulation.
t_end = 0.5 * 3600; % [s] Stop-Zeit der Simulation
n_t = 1000;       % [-] Anzahl der Punkte für die Resultat-Matrix 
tspan = linspace(t_start, t_end, n_t); % [s] Vektor mit Zeitpunkten für die Resultat-Matrix

% Genauigkeit der numerischen Integration
options = odeset('AbsTol', 1e-10);  % Absolutes Abbruchkriterium
options = odeset(options, 'RelTol', 1e-6); % Relatives Abbruchkriterium

%% Aufrufen des ODE-Solvers zur numerischen Integration des Reaktor-Modells
[Sim.t, Sim.y] = ode23s(@RM_SBR_AnBm_it, tspan, y0, options, MP);

%% Resultate der numerischen Integration 
Sim.n_t_i = Sim.y(:, 1:end-1);               % [mol]
Sim.m_RM_t = Sim.y(:, end);                  % [kg]
% Reaktionsvolumen als Funktion der Zeit
Sim.V_RM_t = Sim.m_RM_t / rho_RM * 1000;     % [l]
% Konzentrations-Zeit-Kurven
Sim.c_t_i = Sim.n_t_i ./ Sim.V_RM_t;         % [mol/l]
 
% Umsatz der stöchiometrisch limitierenden Komponente A 
% Achtung: Dieser lässt sich für den Semi-Batch-Reaktor nicht mehr direkt 
% nach Gl.(4.22) berechnen.
i = 1;  % Index der stöchiometrisch limitierenden Komponente, hier "A"
n_i_0_t = n_i_0(i) + n_dos_i(i);  % Gesamte Zugabe an limitierender Komponente 
i_P = 3;          % Index der Produkt-Komponente, in diesem Fall "C"
Sim.X_A_t = 1-(n_i_0_t - abs(nu_i(i)/nu_i(i_P))*Sim.n_t_i(:, i_P)) ./ n_i_0_t; % [mol/mol]

%% Grafische Darstellung der Resultate - Abbildung 7.5
figure
plot(Sim.t/3600, Sim.m_RM_t)
ylabel('m_{RM} [kg]');
grid on; grid minor
xlabel('Zeit [h]');
title(['Abb. 7.5: Reaktion bei 25 °C im Semibatch']);
legend('m_{RM} Reaktionsmasse');
 
figure
plot(Sim.t/60, Sim.c_t_i(:, 1), 'k-', 'linewidth', 2);
hold on
plot(Sim.t/60, Sim.c_t_i(:, 2), 'k', 'linewidth', 0.5);
plot(Sim.t/60, Sim.c_t_i(:, 3), 'k:', 'linewidth', 1);
plot(Sim.t/60, Sim.c_t_i(:, 4), 'k-.', 'linewidth', 1);
grid on
grid minor
xlabel('Zeit [h]');
ylabel('Konzentrationen [mol/l]');
title(['Konz-Zeit Verläufe 60 °C im Semibatch']);
yyaxis right
plot(Sim.t/60, Sim.X_A_t, 'linewidth', 2);
ylabel('X_A [-]');
legend('c_A', 'c_B', 'c_C', 'c_D', 'X_A', 'location', 'east');

%% Anforderungen
% Gewünschter Umsatz bezüglich der Komponente A
X_A_f = 0.99; % [-] 
% Anforderungen an die Anlage:
Cap_BR = 125; % [t/a] Produktionskapazität
Verf = 7500; % [h/a] Verfügbarkeit der Anlage

%% Besimmung der erforderlichen Batch-Zeit für den gewünschten Umsatz X_A_f
% Durch grafischen Ablesen aus Abbildung 7.5
% Index bei welchem Sim.X_A_t > X_A_f erfüllt ist
tmp = find(Sim.X_A_t >= X_A_f);
i_X_Af = tmp(1);

t_Batch = Sim.t(i_X_Af); % [s]
disp(['Batch-Zeit = ', num2str(t_Batch/3600, '%.2g'), ' h']);
disp(['Batch-Zeit = ', num2str(t_Batch/60, '%.3g'), ' min']);
disp(['Batch-Zeit = ', num2str(t_Batch, '%.4g'), ' sec']);

%% Erforderliche Menge Produkt pro Charge gemäss Gl. (7.15):
m_Prod_XAf = Cap_BR *1000 / Verf * t_Batch/3600; % [kg]
disp(['m_Prod = ', num2str(m_Prod_XAf, '%.2g'), ' kg']);

%% Minimal erforderliches Reaktorvolumen gemäss Gl.(7.31):
i_P = 3;  % Index der Produkt-Komponente C
mw_Prod = mw_i(i_P); % [g/mol] Molmasse des Produktes
V_R = m_Prod_XAf * 1000 /mw_Prod /Sim.c_t_i(i_X_Af, i_P); % [l]
disp(['V_R = ', num2str(V_R/1000, '%.2g'), ' m3']);

%% Volumen-Zeit-Ausbeute gemäss Gl.(7.14):
VT_Yield = Cap_BR * 1000 /V_R; % [kg/(l*a)]
disp(['VT_Yield = ', num2str(VT_Yield, '%.2g'), ' kg/(l*a)']);


%% e.)
% Abbildung 1 im Aufgabenstellung entspricht dem dosierten Reaktanden AcOAc
% k schein in Ordnung zu sein

%% f.)
% 1 scheint am besten zu passen, obwohl auf der Aufgaben stellung die
% Konzentrationen um Faktor 10 höher liegen

%% g.) siehe ab Zeile 128

%% h.) 

%% Reaktionsparameter

k_2AB = 0.1;     % [l/(mol*s)] Geschwindigkeitskonstante bei 25 °C, 2. Ordnung
nu_i = [-1, -1, 1, 1, 0]; % [-] Vektor mit stöchiometrische Koeffizienten der Komponenten
n = 1;       % [-] Teilreaktionsordnung der Komponente "A"
m = 1;       % [-] Teilreaktionsordnung der Komponente "B"

%% Anforderungen an die Reaktion und den Reaktor

% Anforderungen an die Anlage:

Cap_CR = 125; % [t/a] Produktionskapazität
Verf = 7500;    % [h/a] Verfügbarkeit der Anlage

% Gewünschter Umsatz bezüglich der Komponente A am Ausgang des Reaktors
X_A_f = 0.99; % [-]

%% Feed-Konzentrationen

c_A_in = n_A/0.0202; % [mol/l] SA
c_B_in = n_B; % [mol] AcOAc rein

%% Aufgabe 7.15, Lösungsvariante 1: Berechnung der mittleren Verweilzeit via symbolische Integration

% Erforderliche mittlere Verweil-Zeit nach Gl. (7.63):
c_A_f = c_A_in*(1-X_A_f);    % [mol/l] erforderliche Endkonzentration
tau_1 = log(c_B_in*c_A_f/(c_A_in*(c_B_in-c_A_in+c_A_f)))/...
    (k_2AB*(c_A_in-c_B_in)); %[s] mittlere Verweilzeit
disp(['PFR Lösungsvariante 1: tau = ', num2str(tau_1/60, '%.4g'), ' min']);

% Minimales Reaktorvolumen "V_R" (7.35) auflösen nach V_R -> (7.48)
mw_Prod = mw_i(3); % Molmasse des Produktes, Komponente "C"
a = 1; % [-] Stöchiometrischer Koeffizient der limitierenden Komponente (A)
c = 1; % [-] Stöchiometrischer Koeffizient des Produktes (C)
V_R = Cap_CR * 1000 * tau_1 / 3600 / (c_A_in*mw_Prod*Verf* X_A_f) * 1000 * a/c; %[l]
disp(['PFR V_R = ', num2str(V_R/1000, '%.2g'), ' m3']);

% Volumen-Zeit-Ausbeute gemäss Gl. (7.14):
VT_Yield = Cap_CR * 1000 /V_R; % [kg/(l*a)]
disp(['PFR VT_Yield = ', num2str(VT_Yield, '%.2g'), ' kg/(l*a)']);

% Volumenstrom in den Reaktor gemäss (7.32)
Vf_in = V_R/tau_1;   % [l/s]
disp(['PFR Vf_in = ', num2str(Vf_in, '%.2g'), ' l/s']);

%% Aufgabe 7.15, Lösungsvariante 2: Berechnung der mittleren Verweilzeit via einfache grafische Integration

n_X = 100;
X_A = linspace(0, X_A_f, n_X);
c_A = c_A_in*(1-X_A);       % [mol/l] Konzentration "A" als Funktion von X_A
c_B = c_B_in - c_A_in*X_A;  % [mol/l] Konzentration "B" als Funktion von X_A
r = k_2AB*c_A.*c_B;         % [mol/(l*s)] Reaktionsgeschwindigkeit (2.13)
r_A = -r;                   % [mol/(l*s)] Umsatzgeschwindigkeit bezüglich "A" (7.1)

inv_r_A = 1./r_A;           % [(l*s)/mol] Inverse Umsatzgewschindigkeit

% Berechnen der Fläche unter der Kurve -1/r_A vs. X_A anhand der Trapez-Methode
% Flaeche_1 = sum(X_A_f/(n_X-1)*-mean([inv_r_A(1:end-1); inv_r_A(2:end)]));
Flaeche_1 = trapz(X_A, -inv_r_A)

% Berechnung der mittleren Verweilzeit gemäss (7.61)
tau_2a= c_A_in * Flaeche_1; %[s]
disp(['PFR Lösungsvariante 2 mit Gl. (7.61): tau = ', num2str(tau_2a/60, '%.4g'), ' min']);

% Berechnen der Fläche unter der Kurve -1/r_A vs. c_A anhand der Trapez-Methode
% Flaeche_2 = sum((c_A_in-c_A_f)/(n_X-1)*-mean(-[inv_r_A(1:end-1); inv_r_A(2:end)]))
Flaeche_2 = trapz(c_A, -inv_r_A)

% Berechnung der mittleren Verweilzeit gemäss (7.62)
tau_2b = -Flaeche_2; %[s]
disp(['PFR Lösungsvariante 2 mit Gl. (7.62): tau  = ', num2str(tau_2b/60, '%.4g'), ' min']);

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
title('PFR - analog zu Gl. (7.21)', 'interpreter', 'latex');
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
title('PFR Gl. (7.62) - analog zu Gl. (7.20)', 'interpreter', 'latex');
XLim = get(gca, 'XLim');
set(gca, 'XLim', [0 XLim(2)]);




