%% AcOEt_Kapitel 7.2 - Batch Reaktor
% Modul Chemische Kinetik und Reaktionstechnik.
% Teil: Reaktionstechnik
%
% Fallstudie Ethylacetat. 
%
% _A. Zogg, 11.02.2021_

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
% 5 = F: Lösungsmittel = Wasser

%% Anfangsbedingungen, gemäss Kapitel 6.5.2 (Aufgabe 6.2)
% Volumen der Reaktionsmasse [l]
V_RM_0 = (25 +25)/1000; % [l]

% Anfangskonzentrationen [mol/l]
c_A_0 = 197/1000*0.894 / 88.11 /V_RM_0; % [mol/l] AcOEt
c_B_0 = 25/1000*0.1 / V_RM_0; % [mol/l] NaOH
c_C_0 = 0; % [mol/l] NaOAc
c_D_0 = 0; % [mol/l] EtOH

%% Erforderliche Batch-Zeit 
X_A_f = 0.9; % [-] finaler Umsatz bezüglich Komponente A

% Erforderliche Batch-Zeit gemäss Gl. (7.24):
k_2AB = 0.1; % [l/(mol*s)] Geschwindigkeitskonstante bei 25 °C gemäss Kapitel 6.6
c_A_f = c_A_0*(1-X_A_f); % [mol/l] erforderliche Endkonzentration
t_Batch = log(c_B_0*c_A_f/(c_A_0*(c_B_0-c_A_0+c_A_f)))/...
(k_2AB*(c_A_0-c_B_0)); %[s] Batch-Zeit
disp(['Batch-Zeit = ', num2str(t_Batch/3600, '%.2g'), ' h']);
disp(['Batch-Zeit = ', num2str(t_Batch/60, '%.3g'), ' min']);
disp(['Batch-Zeit = ', num2str(t_Batch, '%.4g'), ' sec']);

% Anforderungen an die Anlage:
Cap_BR = 100; % [t/a] Produktionskapazität
Verf = 8000; % [h/a] Verfügbarkeit der Anlage

% Erforderliche Menge Produkt pro Charge gemäss Gl. (7.15):
m_Prod_XAf = Cap_BR *1000 / Verf * t_Batch/3600; % [kg]
disp(['m_Prod = ', num2str(m_Prod_XAf, '%.2g'), ' kg']);

% Minimal erforderliches Reaktorvolumen Gl. (7.23):
mw_Prod = 82.03; % [g/mol] Molmasse der Komponente C (Produkt)
a = 1; % [-] Stöchiometrischer Koeffizient der limitierenden Komponente (A)
c = 1; % [-] Stöchiometrischer Koeffizient der Komponente C (Produkt)
V_R = m_Prod_XAf * 1000 /mw_Prod /X_A_f/c_A_0*a/c; % [l]
disp(['V_R = ', num2str(V_R/1000, '%.2g'), ' m3']);

% Volumen-Zeit-Ausbeute gemäss Gl. (7.14):
VT_Yield = Cap_BR * 1000 /V_R; % [kg/(l*a)]
disp(['VT_Yield = ', num2str(VT_Yield, '%.2g'), ' kg/(l*a)']);

