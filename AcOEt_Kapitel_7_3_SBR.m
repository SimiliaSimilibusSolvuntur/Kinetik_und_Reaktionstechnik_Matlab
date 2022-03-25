%% AcOEt_Kapitel 7.3 - Semi-Batch Reaktor
% Modul Chemische Kinetik und Reaktionstechnik.
% Teil: Reaktionstechnik
%
% Fallstudie Ethylacetat. 
%
% _A. Zogg, 11.04.2021_

%% Matlab Konfiguration
clearvars; close all; clc;      
format COMPACT

%% Physikalische Konstanten
R = 8.314;       % [J/(mol*K)] ideale Gaskonstante

%% Komponenten
% 1 = A: Ethylacetat
% 2 = B: NaOH
% 3 = C: Natriumacetat
% 4 = D: Ethanol
% 5 = F: L�sungsmittel = Wasser

%% Versuchsbedingungen gem�ss Kapitel 6.5.2 (Aufgabe 6.2)
% Volumen der Reaktionsmasse [l]
V_RM_0 = (25 +25)/1000; % [l]

% Anfangskonzentrationen [mol/l]
c_A_0 = 197/1000*0.894 / 88.11 /V_RM_0; % [mol/l] SA
c_B_0 = 0; % [mol/l] AcOAc rein
c_C_0 = 0; % [mol/l] 
c_D_0 = 0; % [mol/l] 

%% Stoffdaten der reinen Komponenten
% Molmassen der einzelnen Komponenten
mw_i = [88.11, 40, 82.03, 46.07, 18.02]; % [g/mol] 

%% Stoffdaten der Reaktionsmasse
% Dichte der Reaktionsmasse - wird als konstant angenommen
rho_RM = 1000;   % [kg/m3] 
 
%% Reaktions-Parameter
%% Durch ausprobieren k_2AB herausfinden
k_2AB = 0.1;     % [l/(mol*s)] Geschwindigkeitskonstante bei 25 �C, 2. Ordnung
nu_i = [-1, -1, 1, 1, 0]; % [-] Vektor mit st�chiometrische Koeffizienten der Komponenten 
n = 1;       % [-] Teilreaktionsordnung der Komponente "A"
m = 1;       % [-] Teilreaktionsordnung der Komponente "B"

%% Anfangsbedingungen zum Zeitpunkt t = 0 im Reaktor (=vorgelegt)
V_RM_0 = 1;                         % [l] Reaktionsvolumen 
m_RM_0 = V_RM_0*rho_RM/1000;     % [kg] Reaktionsmasse 
n_B_0 = c_B_0*V_RM_0;               % [mol] Stoffmenge Komponente B 
n_Lsm_0 = (m_RM_0*1000-n_B_0*mw_i(2))/mw_i(5); % [mol] Stoffmenge Lsm
n_i_0 = [0, n_B_0, 0, 0, n_Lsm_0];  % [mol] Stoffmenge Komponente i
 
%% Reaktionsparameter zur Beschreibung der Dosierung
t_dos_start = 0; % [s] Startzeit der Dosierung
t_dos = 10*60;   % [s] Dauer der Dosierung
n_dos_A =  c_A_0 / c_B_0 * n_B_0; % [mol] Stoffmenge Komponente A dosiert
n_dos_i =    [n_dos_A, 0, 0, 0, 0]; % [mol] Stoffmenge Komponente i dosiert
F_dos_i = n_dos_i / t_dos;  %  [mol/s] Stoffstrom w�hrend der Dosierung 

%% Definieren der strukturierten Varialbe "MP" mit den Modell-Parametern f�r des Reaktor-Modell
MP.mw_i = mw_i;         % [g/mol] Vektor mit Molmassen der einzelnen Komponenten
MP.k = k_2AB;           % [...] Geschwindigkeitskonstante 
MP.nu_i = nu_i;         % [-] Vektor mit st�chiometrischen Koeffizienten der Komponenten
MP.n = n;               % [-] Teilreaktionsordnung der Komponente "A"
MP.m = m;               % [-] Teilreaktionsordnung der Komponente "B"
MP.rho_RM = rho_RM;     % [kg/m3] Dichte der Reaktionsmasse 
MP.t_dos_start = t_dos_start; % [s] Startzeit der Dosierung
MP.t_dos = t_dos;       % [s] Dauer der Dosierung
MP.F_dos_i = F_dos_i;   % [mol/s] Vektor mit Stoffstr�men w�hrend der Dosierung 

%% Anfangswerte der zeitabh�ngigen Gr�ssen im Zeilenvektor "y" 
y0 = [n_i_0, m_RM_0];

%% Bedingungen f�r die Simulation
t_start = 0;      % [s] Start-Zeit der Simulation.
t_end = 2 * 3600; % [s] Stop-Zeit der Simulation
n_t = 1000;       % [-] Anzahl der Punkte f�r die Resultat-Matrix 
tspan = linspace(t_start, t_end, n_t); % [s] Vektor mit Zeitpunkten f�r die Resultat-Matrix

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
 
% Umsatz der st�chiometrisch limitierenden Komponente A 
% Achtung: Dieser l�sst sich f�r den Semi-Batch-Reaktor nicht mehr direkt 
% nach Gl.(4.22) berechnen.
i = 1;  % Index der st�chiometrisch limitierenden Komponente, hier "A"
n_i_0_t = n_i_0(i) + n_dos_i(i);  % Gesamte Zugabe an limitierender Komponente 
i_P = 3;          % Index der Produkt-Komponente, in diesem Fall "C"
Sim.X_A_t = 1-(n_i_0_t - abs(nu_i(i)/nu_i(i_P))*Sim.n_t_i(:, i_P)) ./ n_i_0_t; % [mol/mol]

%% Grafische Darstellung der Resultate - Abbildung 7.5
figure
plot(Sim.t/3600, Sim.m_RM_t)
ylabel('m_{RM} [kg]');
grid on; grid minor
xlabel('Zeit [h]');
title(['Abb. 7.5: Reaktion bei 25 �C im Semibatch']);
legend('m_{RM} Reaktionsmasse');
 
figure
plot(Sim.t/3600, Sim.c_t_i(:, 1), 'k-', 'linewidth', 2);
hold on
plot(Sim.t/3600, Sim.c_t_i(:, 2), 'k', 'linewidth', 0.5);
plot(Sim.t/3600, Sim.c_t_i(:, 3), 'k:', 'linewidth', 1);
plot(Sim.t/3600, Sim.c_t_i(:, 4), 'k-.', 'linewidth', 1);
grid on
grid minor
xlabel('Zeit [h]');
ylabel('Konzentrationen [mol/l]');
title(['Abb. 7.5: Reaktion bei 25 �C im Semibatch']);
yyaxis right
plot(Sim.t/3600, Sim.X_A_t, 'linewidth', 2);
ylabel('X_A [-]');
legend('c_A', 'c_B', 'c_C', 'c_D', 'X_A', 'location', 'east');

%% Aufgabe 7.4: Anforderungen
% Gew�nschter Umsatz bez�glich der Komponente A
X_A_f = 0.9; % [-] 
% Anforderungen an die Anlage:
Cap_BR = 100; % [t/a] Produktionskapazit�t
Verf = 8000; % [h/a] Verf�gbarkeit der Anlage

%% Aufgabe 7.4: Besimmung der erforderlichen Batch-Zeit f�r den gew�nschten Umsatz X_A_f
% Durch grafischen Ablesen aus Abbildung 7.5
% Index bei welchem Sim.X_A_t > X_A_f erf�llt ist
tmp = find(Sim.X_A_t >= X_A_f);
i_X_Af = tmp(1);

t_Batch = Sim.t(i_X_Af); % [s]
disp(['Batch-Zeit = ', num2str(t_Batch/3600, '%.2g'), ' h']);
disp(['Batch-Zeit = ', num2str(t_Batch/60, '%.3g'), ' min']);
disp(['Batch-Zeit = ', num2str(t_Batch, '%.4g'), ' sec']);

%% Aufgabe 7.4: Erforderliche Menge Produkt pro Charge gem�ss Gl. (7.15):
m_Prod_XAf = Cap_BR *1000 / Verf * t_Batch/3600; % [kg]
disp(['m_Prod = ', num2str(m_Prod_XAf, '%.2g'), ' kg']);


%% Aufgabe 7.4: Minimal erforderliches Reaktorvolumen gem�ss Gl.(7.31):
i_P = 3;  % Index der Produkt-Komponente C
mw_Prod = mw_i(i_P); % [g/mol] Molmasse des Produktes
V_R = m_Prod_XAf * 1000 /mw_Prod /Sim.c_t_i(i_X_Af, i_P); % [l]
disp(['V_R = ', num2str(V_R/1000, '%.2g'), ' m3']);

%% Aufgabe 7.4: Volumen-Zeit-Ausbeute gem�ss Gl.(7.14):
VT_Yield = Cap_BR * 1000 /V_R; % [kg/(l*a)]
disp(['VT_Yield = ', num2str(VT_Yield, '%.2g'), ' kg/(l*a)']);

