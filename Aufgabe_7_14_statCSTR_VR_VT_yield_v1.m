clearvars; close all; clc;
format COMPACT

% lösung 23 min

mw_i = [200, 100, 160, 70]; % [g/mol]
dichte = 1000; % [kg/m^3]

c_A_in = 2 / 2.5; % [mol/l] Komponente A
c_B_in = 2 / 2.5 * 1.5; % [mol/l] Komponente B

% A + B -> C + 2D

cap_CR = 100; % [t/a]
verf = 8000; % [h/a]
k_2 = 0.005; %[(l/mol)*(1/s)]

disp('Auslegung gemäss Performance-Gleichungen für den stationären CSTR-Reaktor');
% Anforderungen an den Reaktionsumsatz im stationären Zustand
X_A_f = 0.8; % [-] Umsatz bezüglich Komponente A

% Konzentrationen in stationären Zustand beim Umsatz X_A
c_A = c_A_in * (1 - X_A_f);     % [mol/l] Komponente A
c_B = c_B_in - c_A_in * X_A_f;    % [mol/l] Komponente B

% Reaktionsgeschwindigkeit gemäss Gl.(4.17)
r = k_2 * c_A * c_B; % [mol/(s*l)]

% Umsatzrate der söchiometrisch limitierenden Komponente A gemäss Gl.(7.1)
% beim gewünschten stationären Umsatz X_A_f:
r_A = -1 * r    % [mol/(s*l)]

% Mittlere Verweilzeit gemäss Gl. (7.53)
tau = c_A_in * X_A_f / -r_A;
disp(['tau = ', num2str(tau/60, '%.3g'), ' min']);

% Minimales Reaktorvolumen "V_R" (7.48)
a = 1; % [-] Stöchiometrischer Koeffizient der limitierenden Komponente (A)
d = 2; % [-] Stöchiometrischer Koeffizient des Produktes (D)
V_R = cap_CR * 1000 * tau / 3600 / (c_A_in * mw_i(4) * verf * X_A_f) * 1000 * a/d; %[l]
disp(['V_R = ', num2str(V_R/1000, '%.2g'), ' m3']);

% Volumen-Zeit-Ausbeute gemäss Gl. (7.14):
VT_Yield = cap_CR * 1000 / V_R; % [kg/(l*a)]
disp(['VT_Yield = ', num2str(VT_Yield, '%.2g'), ' kg/(l*a)']);



















