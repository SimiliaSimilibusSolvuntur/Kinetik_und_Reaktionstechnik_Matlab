%% AcOEt_Kapitel 7.4 - CSTR Reaktor
% Modul Chemische Kinetik und Reaktionstechnik.
% Teil: Reaktionstechnik
%
% Dimensionieurng des CSTR-Reaktors gem�ss Kapitel 7.4. Fallstudie Ethylacetat. 
%
% Der Reaktor ist am Anfang mit L�sungsmittel gef�llt.
%
% Bedienung zur Dimensionierung (ohne Performance-Gleichung):
%
% # Variieren Sie die Verweilzeit tau, bis Sie im station�ren
% Zustand den gew�nschten Umsatz von X_A = z.B. 90% erreichen. 
%       -> tau = ...; 
% # Durch manuelle Iteration finden Sie 108 Minuten als die gew�nschte Verweilzeit. 
% Der station�re Zustand stellt sich dann allerdings erst nach ca. 10 Stunden ein.
% Verl�ngern Sie deshalb die Simulationszeit auf 10 Stunden: 
%       -> t_end = 10*3600; 
%
% _A. Zogg, 04.05.2021_

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
% Volumen der Reaktionsmasse im Laborexperiment [l]
V_RM_0 = (25 +25)/1000; % [l]

% Anfangskonzentrationen [mol/l]
c_A_0 = 197/1000*0.894 / 88.11 /V_RM_0; % [mol/l] AcOEt
c_B_0 = 25/1000*0.1 / V_RM_0; % [mol/l] NaOH
c_C_0 = 0; % [mol/l] NaOAc
c_D_0 = 0; % [mol/l] EtOH

%% Anforderungen bez�glich der Produktion 
% Anforderungen an die Anlage:
Cap_CR = 100; % [t/a] Produktionskapazit�t
Verf = 8000;    % [h/a] Verf�gbarkeit der Anlage

% Gew�nschter Umsatz bez�glich der Komponente A am Ausgang des Reaktors
X_A_out = 0.9; % [-] 

%% Durch ausprobieren mittlere Verweilzeit
tau = 10*60;            % [s] {10 oder 108 Minuten}
disp(['Reaktion bei 25 �C im CSTR mit tau = ', num2str(tau/3600, '%.2g'), 'h']);

%% Stoffdaten der reinen Komponenten
% Molmassen der einzelnen Komponenten
mw_i = [88.11, 40, 82.03, 46.07, 18.02]; % [g/mol] 
 
%% Reaktions-Parameter
k_2AB = 0.1;     % [l/(mol*s)] Geschwindigkeitskonstante bei 25 �C, 2. Ordnung
nu_i = [-1, -1, 1, 1, 0]; % [-] Vektor mit st�chiometrische Koeffizienten der Komponenten 
n = 1;       % [-] Teilreaktionsordnung der Komponente "A"
m = 1;       % [-] Teilreaktionsordnung der Komponente "B"
 
%% Feed-Konzentrationen
rho_in = 1000;          % [kg/m3] Dichte des Feeds in den Reaktor
c_Lsm_i = (rho_in-c_A_0*mw_i(1)-c_B_0*mw_i(2))/mw_i(5);
c_i_in = [c_A_0, c_B_0, 0, 0, c_Lsm_i]; % [mol/l] Feed-Konzentrationen
 
%% Reaktor-Parameter
% Minimales Reaktorvolumen "V_R" dann �berlauf gem�ss (7.48) 
a = 1; % [-] St�chiometrischer Koeffizient der limitierenden Komponente (A)
c = 1; % [-] St�chiometrischer Koeffizient des Produktes (C)
mw_Prod = mw_i(3); % [g/mol] Molmasse der Komponente C (Produkt)
V_R = Cap_CR * 1000 * tau / 3600 / (c_i_in(1)*mw_Prod*Verf* X_A_out) * 1000 * a/c; %[l]
disp(['Minimal erforderliches Reaktorvolumen V_R = ', num2str(V_R, '%.2g'), ' l']);

%% Anfangsbedingungen zum Zeitpunkt t = 0 im Reaktor (=vorgelegt)
rho_RM_0 = 1000;                  % [kg/m3] Dichte der vorgelegten Reaktionsmasse
V_RM_0 = V_R;                     % [l] Reaktionsvolumen 
m_RM_0 = V_RM_0*rho_RM_0/1000;    % [kg] Reaktionsmasse 
n_Lsm_0 = (m_RM_0*1000)/mw_i(5);  % [mol] Stoffmenge Lsm
n_i_0 = [0, 0, 0, 0, n_Lsm_0];    % [mol] Stoffmenge Komponente i
 
%% Modell-Parameter zur Beschreibung des kontinuierlichen Feeds
Vf_in = V_R/tau;   % [l/s] Volumenstrom in den Reaktor gem�ss (5.32)
disp(['Vf_in = ', num2str(Vf_in, '%.2g'), ' l/s']);
mf_in = Vf_in * rho_in /1000; % [kg/s] Gesamtmassenstrom in den Reaktor (5.43)
F_i_in = c_i_in.*Vf_in; %[mol/s] Stoffstr�me in den Reaktor

%% Definieren der strukturierten Varialbe "MP" mit den Modell-Parametern f�r des Reaktor-Modell
MP.mw_i = mw_i;         % [g/mol] Vektor mit Molmassen der einzelnen Komponenten
MP.k = k_2AB;           % [...] Geschwindigkeitskonstante 
MP.nu_i = nu_i;         % [-] Vektor mit st�chiometrischen Koeffizienten der Komponenten
MP.n = n;               % [-] Teilreaktionsordnung der Komponente "A"
MP.m = m;               % [-] Teilreaktionsordnung der Komponente "B"
MP.F_i_in = F_i_in;     % [mol/s] Vektort mit Stoffstr�men der einzelnen Komponenten in den Reaktor
MP.rho_in = rho_in;     % [kg/m3] Dichte des Feeds in den Reaktor 
MP.V_R = V_R;           % [l] Volumen des Reaktors 

%% Anfangswerte der zeitabh�ngigen Gr�ssen im Zeilenvektor "y" 
y0 = [n_i_0, m_RM_0, V_RM_0];   
 
%% Bedingungen f�r die Simulation
t_start = 0;      % [s] Start-Zeit der Simulation.
t_end = 2 * 3600; % [s] Stop-Zeit der Simulation {2 oder 10 h}
n_t = 1000;       % [-] Anzahl der Punkte f�r die Resultat-Matrix 
tspan = linspace(t_start, t_end, n_t); % [s] Vektor mit Zeitpunkten f�r die Resultat-Matrix
 
% Genauigkeit der numerischen Integration
options = odeset('AbsTol', 1e-10);  % Absolutes Abbruchkriterium
options = odeset(options, 'RelTol', 1e-6); % Relatives Abbruchkriterium
 
%% Aufrufen des ODE-Solvers zur numerischen Integration des Reaktor-Modells
[Sim.t, Sim.y] = ode23s(@RM_CSTR_AnBm_it, tspan, y0, options, MP);


%% Resultate der numerischen Integration
Sim.n_t_i = Sim.y(:, 1:end-2);  % [mol]
Sim.m_RM_t = Sim.y(:, end-1);   % [kg]
Sim.V_RM_t = Sim.y(:, end);     % [l]

% Konzentrations-Zeit-Kurven
Sim.c_t_i = Sim.n_t_i ./ Sim.V_RM_t;         % [mol/l]
 
% Umsatz der s�chiometrisch limiteirenden Komponente A gem�ss Gl.(5.49)
i = 1;
Sim.X_A_t = (F_i_in(i) - Sim.c_t_i(:, i)*Vf_in) ./ F_i_in(i); % [mol/mol]

% Umsatz der st�chiometrisch limitierenden Komponente A 
i = 1;  % Index der st�chiometrisch limitierenden Komponente, hier "A"
i_P = 3;          % Index der Produkt-Komponente, in diesem Fall "C"
Sim.X_A_t = 1-(F_i_in(i) - abs(nu_i(i)/nu_i(i_P))*Sim.c_t_i(:, i_P)*Vf_in) ./ F_i_in(i); % [mol/mol]

% Alternative:
% Umsatz der s�chiometrisch limiteirenden Komponente A gem�ss Gl.(4.10)
% Sim.X_A_t = (c_i_in(1) - Sim.c_t_i(:, 1)) / c_i_in(1); % [mol/mol]

% Station�rer Umsatz
disp(['Station�rer Umsatz = ', num2str(Sim.X_A_t(end), '%.2g'), ' -']);

%% Grafische Darstellung der Resultate - Abbildung 7.6
figure
plot(Sim.t/3600, Sim.m_RM_t)
ylabel('m_{RM} [kg]');
grid on; grid minor
xlabel('Zeit [h]');
title(['Abb. 7.6: Reaktion im CSTR mit \tau = ', num2str(tau/3600, '%.2g'), ' h']);
ytickformat('%.0g')
legend(['m_{RM} Reaktionsmasse', newline, ...
    'Reaktorvolumen V_R = ', num2str(V_R/1000, '%0.2g'), ' m3'], ...
    'location', 'best');
 
figure
grid on; grid minor; hold on
plot(Sim.t/3600, Sim.c_t_i(:, 1), 'k-', 'linewidth', 2);
plot(Sim.t/3600, Sim.c_t_i(:, 2), 'k', 'linewidth', 0.5);
plot(Sim.t/3600, Sim.c_t_i(:, 3), 'k:', 'linewidth', 1);
plot(Sim.t/3600, Sim.c_t_i(:, 4), 'k-.', 'linewidth', 1);
xlabel('Zeit [h]');
ylabel('Konzentrationen [mol/l]');
title(['Abb. 7.6: Reaktion im CSTR mit \tau = ', num2str(tau/3600, '%.2g'), ' h']);
yyaxis right
plot(Sim.t/3600, Sim.X_A_t, 'linewidth', 2);
ylabel('X_A [-]');
legend('c_A', 'c_B', 'c_C', 'c_D', 'X_A', 'location', 'best');


%% Aufgabe 7.9: Resultierende Produktionskapazit�t 
i = 3; % Index der Produktkomponente C

mf_out_Prod = Sim.c_t_i(end, i) * Vf_in * mw_i(i) / 1000; % [kg]
Cap_ist = mf_out_Prod * 3600 * Verf/1000; % [t/a] (7.34)
disp(['Cap_ist = ', num2str(Cap_ist, '%.3g'), ' t/a']);

% Volumen-Zeit-Ausbeute gem�ss Gl. (7.14):
VT_Yield = Cap_ist * 1000 /V_R; % [kg/(l*a)] 
disp(['VT_Yield = ', num2str(VT_Yield, '%.2g'), ' kg/(l*a)']);

%% Kapitel 7.4.2: Anwendung der Performance Gleichungen zum station�ren CSTR-Reaktor 
disp('Auslegung gem�ss Performance-Gleichungen f�r den station�ren CSTR-Reaktor');
% Anforderungen an den Reaktionsumsatz im station�ren Zustand
X_A_f = X_A_out; % [-] Umsatz bez�glich Komponente A 

% Konzentrationen im Feed sind identisch mit den Anfangskonzentrationen aus
% Quelle [2]
c_A_in = c_A_0; % [mol/l] Komponente A
c_B_in = c_B_0; % [mol/l] Komponente B

% Konzentrationen in station�ren Zustand beim Umsatz X_A 
c_A = c_A_in *(1 - X_A_f);     % [mol/l] Komponente A
c_B = c_B_in - c_A_in*X_A_f;    % [mol/l] Komponente B

% Reaktionsgeschwindigkeit gem�ss Gl.(4.17)
r = k_2AB*c_A*c_B; % [mol/(s*l)] 

% Umsatzrate der s�chiometrisch limitierenden Komponente A gem�ss Gl.(7.1)
% beim gew�nschten station�ren Umsatz X_A_f:
r_A = -1 * r    % [mol/(s*l)] 

% Mittlere Verweilzeit gem�ss Gl. (7.53)
tau = c_A_in*X_A_f/-r_A;
disp(['tau = ', num2str(tau/60, '%.3g'), ' min']);

% Minimales Reaktorvolumen "V_R" (7.48)
a = 1; % [-] St�chiometrischer Koeffizient der limitierenden Komponente (A)
c = 1; % [-] St�chiometrischer Koeffizient des Produktes (C)
V_R = Cap_CR * 1000 * tau / 3600 / (c_A_in*mw_Prod*Verf* X_A_f) * 1000 * a/c; %[l]     
disp(['V_R = ', num2str(V_R/1000, '%.2g'), ' m3']);

% Volumen-Zeit-Ausbeute gem�ss Gl. (7.14):
VT_Yield = Cap_CR * 1000 /V_R; % [kg/(l*a)]
disp(['VT_Yield = ', num2str(VT_Yield, '%.2g'), ' kg/(l*a)']);

% Volumenstrom in den Reaktor gem�ss (7.32)
Vf_in = V_R/tau;   % [l/s] 
disp(['Vf_in = ', num2str(Vf_in, '%.2g'), ' l/s']);
