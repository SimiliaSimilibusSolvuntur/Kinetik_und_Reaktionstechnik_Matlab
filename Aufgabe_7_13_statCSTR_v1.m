clearvars; close all; clc;
format COMPACT

% Lösung 3 min

disp('Auslegung gemäss Performance-Gleichungen für den stationären CSTR-Reaktor');
% Anforderungen an den Reaktionsumsatz im stationären Zustand
X_i_f = 0.75; % [-] Umsatz bezüglich Komponente A

a = 1; b = 2; c = 1;

V_R = 6; % [l]

% Konzentrationen im Feed sind identisch mit den Anfangskonzentrationen aus
% Quelle [2]
c_A_in = 2.8 / 2; % [mol/l] Komponente A
c_B_in = 1.6 / 2; % [mol/l] Komponente B
c_R_theo = 0; % [mol/l]

k_1 = 12.5; % [L^2/(mol^2*min)]
k_2 = 1.5; % [1/min]

% Konzentrationen in stationären Zustand beim Umsatz X_B
c_B_f = c_B_in * (1 - X_i_f);     % [mol/l] Komponente B
c_A_f = c_A_in - (a / b) * c_B_in * X_i_f;    % [mol/l] Komponente A
c_C_f = 0 + (c / b) * c_B_in * X_i_f;    % [mol/l] Komponente B

%% 
% Umsatzrate von B (limitierende Komponente)
r_B = -2 * (k_1 * c_A_f * c_B_f^2 - k_2 * c_C_f); % [mol/(l*min)]
r_A = -(k_1 * c_A_f * c_B_f^2 - k_2 * c_C_f); % [mol/(l*min)]

% Verweilzeit (7.54)
tau_B = (c_B_in - c_B_f)/-r_B % [min]
tau_A = (c_A_in - c_A_f)/-r_A % [min]

% Volumenströme der beiden Feeds (7.32) V_R = V_RM
Vf = V_R / tau_B / 2 % [l/min]


%% 
% % Reaktionsgeschwindigkeit gemäss Gl.(4.17)
% r = k_1 * c_A_f * (c_B_f)^2 - k_2 * c_R; % [mol/(s*l)]
% 
% % Umsatzrate der söchiometrisch limitierenden Komponente B gemäss Gl.(7.1)
% % beim gewünschten stationären Umsatz X_A_f:
% r_A = -1 * r    % [mol/(s*l)]
% 
% % Mittlere Verweilzeit gemäss Gl. (7.53)
% tau = c_B_in * X_B_f / -r_A;
% disp(['tau = ', num2str(tau, '%.3g'), ' s']);
% 
% % Minimales Reaktorvolumen "V_R" (7.48)
% b = 2; % [-] Stöchiometrischer Koeffizient der limitierenden Komponente (A)
% c = 1; % [-] Stöchiometrischer Koeffizient des Produktes (C)
% V_R = 6; %[l]
% disp(['V_R = ', num2str(V_R / 1000, '%.2g'), ' m3']);
% 
% % Volumenstrom in den Reaktor gemäss (7.32)
% Vf_in = V_R / tau;   % [l/s]
% disp(['Vf_in = ', num2str(Vf_in, '%.2g'), ' l/s']);
% 
% Vf_in_A = V_R / tau;   % [l/s]
% disp(['Vf_in = ', num2str(Vf_in_A / 2, '%.2g'), ' l/s']);
% 
% Vf_in_B = V_R / tau;   % [l/s]
% disp(['Vf_in = ', num2str(Vf_in_B / 2, '%.2g'), ' l/s']);



