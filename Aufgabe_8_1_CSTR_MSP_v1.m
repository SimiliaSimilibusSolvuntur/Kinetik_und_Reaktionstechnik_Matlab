clearvars; clc;
close all;

% Dichte = 1 kg/l
% B wird als Reinstoff zudosiert
% E Lösungsmittel
% T = 60°C

%  1 A + 1 B -> 1 C + 2 D
a = 1;
b = 1;
c = 1;
d = 2;

m_i_0 = [200, 0, 0, 0, 650]; % [g] i = 5 lsm
mw_i = [200, 100, 250, 25, 72.11]; % [] i = 5 lsm

% Anfangskonzentrationen [mol/l]
c_A_0 = 1; % [mol]
c_B_0 = 1.5; % [mol] Dosierung rein
c_C_0 = 0; % [mol]
c_D_0 = 0; % [mol]

% Dichte der Reaktionsmasse - wird als konstant angenommen
rho_RM = 1000;   % [kg/m^3]

% Volumen der Reaktionsmasse [L]
V_RM_0 = (0.85)/1; % [L]

%% Physikalische Konstanten
R = 8.314;       % [J/(mol*K)] ideale Gaskonstante

%% Reaktions-Parameter
k_2AB = 0.0035     % [l/(mol*s)] durch ausprobieren, 2. Ordnung
nu_i = [-1, -1, 1, 2, 0]; % [-] Vektor mit stöchiometrische Koeffizienten der Komponenten 
n = 1;       % [-] Teilreaktionsordnung der Komponente "A"
m = 1;       % [-] Teilreaktionsordnung der Komponente "B"

%% Anforderungen bezüglich der Produktion 
% Anforderungen an die Anlage:
Cap_CR = 100; % [t/a] Produktionskapazität
Verf = 8000;    % [h/a] Verfügbarkeit der Anlage

% Gewünschter Umsatz bezüglich der Komponente A am Ausgang des Reaktors
X_A_out = 0.9; % [-] 

% Geschätzte mittlere Verweilzeit
tau = 1.2*60*60;            % [s] {10 oder 108 Minuten}
disp(['Reaktion bei 60 °C im CSTR mit tau = ', num2str(tau/3600, '%.2g'), 'h']);

%% Feed-Konzentrationen
rho_in = 1000;          % [kg/m3] Dichte des Feeds in den Reaktor
% c_Lsm_i = 650/mw_i(5);
c_Lsm_i = (rho_in-c_A_0*mw_i(1)-c_B_0*mw_i(2))/mw_i(5);
c_i_in = [c_A_0, c_B_0, 0, 0, c_Lsm_i] % [mol/l] Feed-Konzentrationen

%% Reaktor-Parameter
% Minimales Reaktorvolumen "V_R" dann Überlauf gemäss (7.48) 
% a = 1; % [-] Stöchiometrischer Koeffizient der limitierenden Komponente (A)
% c = 1; % [-] Stöchiometrischer Koeffizient des Produktes (C)
% mw_Prod = mw_i(3); % [g/mol] Molmasse der Komponente C (Produkt)
% V_R_1 = Cap_CR * 1000 * tau / 3600 / (c_i_in(1)*mw_Prod*Verf* X_A_out) * 1000 * a/c %[l]
% disp(['Minimal erforderliches Reaktorvolumen C V_R_1 = ', num2str(V_R_1, '%.2g'), ' l']);

mw_Prod = mw_i(4); % [g/mol] Molmasse der Komponente D (Produkt)
V_R_2 = Cap_CR * 1000 * tau / 3600 / (c_i_in(1)*mw_Prod*Verf* X_A_out) * 1000 * a/d; %[l]
disp(['Minimal erforderliches Reaktorvolumen D V_R_2 = ', num2str(V_R_2, '%.2g'), ' l']);

V_R = V_R_2; %[l]

%% Anfangsbedingungen zum Zeitpunkt t = 0 im Reaktor (=vorgelegt)
rho_RM_0 = 1000;                  % [kg/m3] Dichte der vorgelegten Reaktionsmasse
V_RM_0 = V_R;                     % [l] Reaktionsvolumen 
m_RM_0 = V_RM_0*rho_RM_0/1000;    % [kg] Reaktionsmasse 
n_Lsm_0 = (m_RM_0*1000)/mw_i(5);  % [mol] Stoffmenge Lsm
n_i_0 = [0, 0, 0, 0, n_Lsm_0];    % [mol] Stoffmenge Komponente i
 
%% Modell-Parameter zur Beschreibung des kontinuierlichen Feeds
Vf_in = V_R/tau;   % [l/s] Volumenstrom in den Reaktor gemäss (5.32)
disp(['Vf_in = ', num2str(Vf_in, '%.2g'), ' l/s']);
mf_in = Vf_in * rho_in / 1000; % [kg/s] Gesamtmassenstrom in den Reaktor (5.43)
F_i_in = c_i_in.*Vf_in; %[mol/s] Stoffströme in den Reaktor

%% Definieren der strukturierten Varialbe "MP" mit den Modell-Parametern für des Reaktor-Modell
MP.mw_i = mw_i;         % [g/mol] Vektor mit Molmassen der einzelnen Komponenten
MP.k = k_2AB;           % [...] Geschwindigkeitskonstante 
MP.nu_i = nu_i;         % [-] Vektor mit stöchiometrischen Koeffizienten der Komponenten
MP.n = n;               % [-] Teilreaktionsordnung der Komponente "A"
MP.m = m;               % [-] Teilreaktionsordnung der Komponente "B"
MP.F_i_in = F_i_in;     % [mol/s] Vektort mit Stoffströmen der einzelnen Komponenten in den Reaktor
MP.rho_in = rho_in;     % [kg/m3] Dichte des Feeds in den Reaktor 
MP.V_R = V_R;           % [l] Volumen des Reaktors 

%% Anfangswerte der zeitabhängigen Grössen im Zeilenvektor "y" 
y0 = [n_i_0, m_RM_0, V_RM_0];   
 
%% Bedingungen für die Simulation
t_start = 0;      % [s] Start-Zeit der Simulation.
t_end = 10*3600; % [s] Stop-Zeit der Simulation {2 oder 10 h}
n_t = 1000;       % [-] Anzahl der Punkte für die Resultat-Matrix 
tspan = linspace(t_start, t_end, n_t); % [s] Vektor mit Zeitpunkten für die Resultat-Matrix

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

% Umsatz der söchiometrisch limiteirenden Komponente A gemäss Gl.(5.49)
i = 1;
Sim.X_A_t = (F_i_in(i) - Sim.c_t_i(:, i)*Vf_in) ./ F_i_in(i); % [mol/mol]

% Alternative:
% Umsatz der söchiometrisch limiteirenden Komponente A gemäss Gl.(4.10)
% Sim.X_A_t = (c_i_in(1) - Sim.c_t_i(:, 1)) / c_i_in(1); % [mol/mol]

% Stationärer Umsatz
disp(['Stationärer Umsatz = ', num2str(Sim.X_A_t(end), '%.2g'), ' -']);

%% Grafische Darstellung der Resultate - Abbildung 7.6
figure
plot(Sim.t/3600, Sim.m_RM_t)
ylabel('m_{RM} [kg]');
grid on; grid minor
xlabel('Zeit [h]');
title(['Abb. 7.6: Reaktion im CSTR mit \tau = ', num2str(tau/3600, '%.2g'), ' h']);
ytickformat('%.0g')
legend(['m_{RM} Reaktionsmasse', newline, ...
    'Reaktorvolumen V_R = ', num2str(V_R/1000, '%0.2g'), ' m3']);

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
legend('c_A', 'c_B', 'c_C', 'c_D', 'X_A');


%% Aufgabe 7.9: Resultierende Produktionskapazität 
i = 4; % Index der Produktkomponente D

mf_out_Prod = Sim.c_t_i(end, i) * Vf_in * mw_i(i) / 1000; % [kg]
Cap_ist = mf_out_Prod * 3600 * Verf/1000; % [t/a] (7.34)
disp(['Cap_ist = ', num2str(Cap_ist, '%.3g'), ' t/a']);

% Volumen-Zeit-Ausbeute gemäss Gl. (7.14):
VT_Yield = Cap_ist * 1000 / V_R; % [kg/(l*a)] 
disp(['VT_Yield = ', num2str(VT_Yield, '%.2g'), ' kg/(l*a)']);











