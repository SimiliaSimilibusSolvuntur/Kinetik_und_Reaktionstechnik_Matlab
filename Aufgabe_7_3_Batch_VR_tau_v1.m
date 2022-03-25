clearvars; close all; clc;
format compact

% Dichte: konstant bei 1000 % [kg/m3]
% Produkt: D
% Kapazität: cap_BR = 100 % =cap [t/a]
% Anlagenverfügbarkeit Zeit pro Jahr: verf = 8000 % [h/a]
mw_i = [200 100 160 70]; % [g/mol] A ... D
a = 1; % [-] Stöchiometrischer Koeffizient der limitierenden Komponente (A)
d = 2; % [-] Stöchiometrischer Koeffizient der Komponente D (Produkt)

%% a.)

% Reaktion 3.2
% Reaktion 3: 2. Ordnung
% A + B -> C + 2D
% c_A != c_B

c_A_0 = 1.25; %[mol/L]
c_B_0 = 1.5; %[mol/L]
c_C_0 = 0; %[mol/L]
c_D_0 = 0; %[mol/L]
AB_01 = -0.5; %[L/mol]
t = 400; %[s]

k_2AB = AB_01 / ((c_A_0 - c_B_0) * t) %[(L/mol)*(1/s)]

X_A_f = 0.8; % [-] finaler Umsatz von A

c_A_f = c_A_0 * (1 - X_A_f); % [mol/l] erforderliche Endkonzentration

t_Batch = log(c_B_0 * c_A_f / (c_A_0 * (c_B_0 - c_A_0 + c_A_f)))/...
(k_2AB * (c_A_0 - c_B_0)) %[s] Batch-Zeit
disp(['Batch-Zeit = ', num2str(t_Batch/3600, '%.2g'), ' h']);
disp(['Batch-Zeit = ', num2str(t_Batch/60, '%.3g'), ' min']);
disp(['Batch-Zeit = ', num2str(t_Batch, '%.4g'), ' sec']);

verf = 8000; % [h/a]

cap_BR = 100; %  [t/a]

% Erforderliche Menge Produkt pro Charge gemäss Gl. (7.15):
m_Prod_XAf = cap_BR * 1000 / verf * t_Batch / 3600; % [kg]
disp(['m_Prod = ', num2str(m_Prod_XAf, '%.2g'), ' kg']);

% Minimal erforderliches Reaktorvolumen Gl. (7.23):
mw_Prod = mw_i(4); % [g/mol] Molmasse der Komponente D (Produkt)
V_R = m_Prod_XAf * 1000 * a / (mw_Prod * X_A_f * c_A_0 * d); % [l]
disp(['V_R = ', num2str(V_R / 1000, '%.2g'), ' m3']);

% Volumen-Zeit-Ausbeute gemäss Gl. (7.14):
VT_Yield = cap_BR * 1000 / V_R; % [kg/(l*a)]
disp(['VT_Yield = ', num2str(VT_Yield, '%.2g'), ' kg/(l*a)']);




















