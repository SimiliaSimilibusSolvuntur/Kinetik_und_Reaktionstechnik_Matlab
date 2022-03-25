clearvars; close all; clc;

%%
%{
Reaktion 4: 
c_A_0 = 1.25
c_B_0 = 1.5
c_C_0 = 0
c_D_0 = 0

1 A + 2 B -> 1 C + 2 D

%}

%% Reaktion 4

%% Physikalische Konstanten

R = 8.314;       % [J/(mol*K)] ideale Gaskonstante

%% Komponenten
%{ 
1 = A: Ethylacetat 2 = B: NaOH 3 = C: Natriumacetat 4 = D: 
Ethanol 5 = F: Lösungsmittel = Wasser
%}

%% Versuchsbedingungen gemäss Kapitel 6.5.2 (Aufgabe 6.2)

% Volumen der Reaktionsmasse [l]
V_RM_0 = (25 + 25)/1000; % [l]

% Anfangskonzentrationen [mol/l]
c_A_0 = 1.25; % [mol/l] AcOEt
c_B_0 = 1.5; % [mol/l] NaOH
c_C_0 = 0; % [mol/l] NaOAc
c_D_0 = 0; % [mol/l] EtOH

%% Stoffdaten der reinen Komponenten

% Molmassen der einzelnen Komponenten

mw_i = [200, 50, 160, 70, 130]; % [g/mol]

%% Stoffdaten der Reaktionsmasse

% Dichte der Reaktionsmasse - wird als konstant angenommen
rho_RM = 1000;   % [kg/m3]

%% Reaktions-Parameter

% t_1 = 500;
% c_A_ln = -0.66;
% k_1 = (c_A_ln - log(c_A_0)) / (-t_1) %[1/s]
% t_1 = 500;
% AB_01 = -0.21;
% k_2_AB2 = AB_01 / ((c_A_0 - c_B_0) * t_1) %[(L / mol) * (1 / s)]
k_2_AB2 = 0.005; %[(L/mol)*(1/s)]
nu_i = [-1, -2, 1, 2, 0]; % [-] Vektor mit stöchiometrische Koeffizienten der Komponenten
n = 1;       % [-] Teilreaktionsordnung der Komponente "A"
m = 1;       % [-] Teilreaktionsordnung der Komponente "B"

%% Anfangsbedingungen zum Zeitpunkt t = 0 im Reaktor (=vorgelegt)

V_RM_0 = 1;                         % [l] Reaktionsvolumen
m_RM_0 = V_RM_0*rho_RM/1000;     % [kg] Reaktionsmasse
n_B_0 = c_B_0 * V_RM_0;               % [mol] Stoffmenge Komponente B
n_A_0 = c_A_0 * V_RM_0;               % [mol] Stoffmenge Komponente A
n_Lsm_0 = (m_RM_0*1000-n_B_0*mw_i(2))/mw_i(5); % [mol] Stoffmenge Lsm
n_i_0 = [n_A_0, n_B_0, 0, 0, n_Lsm_0];  % [mol] Stoffmenge Komponente i

%% Reaktionsparameter zur Beschreibung der Dosierung

t_dos_start = 0; % [s] Startzeit der Dosierung
t_dos = 0;   % [s] Dauer der Dosierung
n_dos_A =  0; % [mol] Stoffmenge Komponente A dosiert
n_dos_i =    [n_dos_A, 0, 0, 0, 0]; % [mol] Stoffmenge Komponente i dosiert
F_dos_i = n_dos_i;  %  [mol/s] Stoffstrom während der Dosierung

%% Definieren der strukturierten Varialbe "MP" mit den Modell-Parametern für des Reaktor-Modell

MP.mw_i = mw_i;         % [g/mol] Vektor mit Molmassen der einzelnen Komponenten
MP.k = k_2_AB2;           % [...] Geschwindigkeitskonstante
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
t_end = 600; % [s] Stop-Zeit der Simulation
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
Sim.X_A_t = 1 - (n_i_0_t - abs(nu_i(i) / nu_i(i_P)) * Sim.n_t_i(:, i_P)) ./ n_i_0_t; % [mol/mol]

%% Grafische Darstellung der Resultate - Abbildung 7.5

figure
plot(Sim.t/3600, Sim.m_RM_t)
ylabel('m_{RM} [kg]');
grid on; grid minor
xlabel('Zeit [h]');
title(['Abb. 7.5: Reaktion bei 25 °C im Semibatch']);
legend('m_{RM} Reaktionsmasse');

figure
plot(Sim.t, Sim.c_t_i(:, 1), 'k-', 'linewidth', 2);
hold on
plot(Sim.t, Sim.c_t_i(:, 2), 'k', 'linewidth', 0.5);
plot(Sim.t, Sim.c_t_i(:, 3), 'k:', 'linewidth', 1);
plot(Sim.t, Sim.c_t_i(:, 4), 'k-.', 'linewidth', 1);
grid on
grid minor
xlabel('Zeit [s]');
ylabel('Konzentrationen [mol/l]');
title(['Abb. 7.5: Reaktion bei 25 °C im Semibatch']);
yyaxis right
plot(Sim.t, Sim.X_A_t, 'linewidth', 2);
ylabel('X_A [-]');
legend('c_A', 'c_B', 'c_C', 'c_D', 'X_A');
