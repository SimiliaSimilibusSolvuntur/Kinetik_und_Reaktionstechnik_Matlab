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
k_2AB = 0.003     % [l/(mol*s)] durch ausprobieren, 2. Ordnung
nu_i = [-1, -1, 1, 2, 0]; % [-] Vektor mit stöchiometrische Koeffizienten der Komponenten 
n = 1;       % [-] Teilreaktionsordnung der Komponente "A"
m = 1;       % [-] Teilreaktionsordnung der Komponente "B"

%% Anforderungen bezüglich der Produktion 
% Anforderungen an die Anlage:
Cap_CR = 100; % [t/a] Produktionskapazität
verf = 8000;    % [h/a] Verfügbarkeit der Anlage

% Gewünschter Umsatz bezüglich der Komponente A am Ausgang des Reaktors
X_A_f = 0.9; % [-] 

%%
c_A_in = c_A_0;
c_B_in = c_B_0;

% Konzentrationen in stationären Zustand beim Umsatz X_A 
c_A = c_A_in * (1 - X_A_f);     % [mol/l] Komponente A gemäss (4.25)
c_B = c_B_in - b/a * c_A_in * X_A_f;    % [mol/l] Komponente B gemäss (4.26)

% Berechnung der Umsatzrate r_A gemäss (4.17) und (4.19)
r = k_2AB * c_A * c_B;
r_A = -r;

% Mittlere Verweilzeit gemäss Gl. (7.53)
tau = c_A_in * X_A_f / -r_A;
disp(['tau = ', num2str(tau/60, '%.3g'), ' min']);

% Minimales Reaktorvolumen "V_R" = V_RM (7.48) mittels Produkt "D"
V_R = Cap_CR * 1000 * tau / 3600 / (c_A_in * mw_i(4) * verf * X_A_f) * 1000 * a/d; %[l]
disp(['V_R = ', num2str(V_R/1000, '%.2g'), ' m3']);

% Volumen-Zeit-Ausbeute gemäss Gl. (7.14):
VT_Yield = Cap_CR * 1000 / V_R; % [kg/(l*a)]
disp(['VT_Yield = ', num2str(VT_Yield, '%.2g'), ' kg/(l*a)']);








