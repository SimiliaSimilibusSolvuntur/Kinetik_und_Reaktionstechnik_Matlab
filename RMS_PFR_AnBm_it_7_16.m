function dy_dt = RMS_PFR_AnBm_it_7_16(l, y, MP)

% Reaktor-Model für einen isothermen, idealen Plug-Flow-Reaktor mit
% konstanter Dichte - analog zu Gleichung (7.62)
%
% Es können folgende Reaktionstypen abgebildet werden:
%   nu_1*A + nu_2*B -> nu_4*C + nu_4*D
%   nu_1*A -> nu_2*B + nu_4*C + nu_4*D
%
% Lösungsmittel = Komponente E (Nr. 5). Es gilt nu_5 = 0.
%
% Funktionsargumente
% l           [m] Rohlänge "l" zwingend das erste Funktionsargument
% y           Vektor mit den ortsabhängigen Variablen
%             Element 1 ... 5: [mol/s] Stoffstrom der Komponente A ... D + Lsm
% MP:         Strukturierte Varialbe "MP" mit den Modell-Parametern
%
% Felder der strukturierten Variablen "MP":
% MP.mw_i        [g/mol] Vektor mit Molmassen der einzelnen Komponenten
% MP.k           [...] Geschwindigkeitskonstante
%                [l/s] Geschwindigkeitskonstante 1. Ordnung
%                [l/(mol*s)] Geschwindigkeitskonstante 2. Ordnung
%                [l^2/(mol^2*s)] Geschwindigkeitskonstante 3. Ordnung
% MP.nu_i        [-] Vektor mit stöchiometrischen Koeffizienten der Komponenten
% MP.n           [-] Teilreaktionsordnung der Komponente "A"
% MP.m           [-] Teilreaktionsordnung der Komponente "B"
% MP.Vf_in       [l/s] Volumenstrom in den Reaktor in den Reaktor
% MP.rho_in      [kg/m3] Dichte des Feeds in den Reaktor
% MP.d_R         [m] Durchmesser des Reaktorrohres

%% Momentan-Werte für die ortsabhängigen Variablen

% Wir setzen alle Momentan-Werte im Spaltenvektor "y" auf >= 0

y(y<0) =0;

% Auslesen des Spaltenvektors "y" mit allen Momentan-Werten
F_i = y';  % [mol/s] Stoffstrom der Komponente A ... D + Lsm

%% Berechnung Konzentrationen

c_i = F_i / MP.Vf_in; % [mol/l] Konzentrationen der einzelnen Komponenten

%% Stoffströme und Stoffstrombilanz

% Reaktionsgeschwindigkeit bzw. Geschwindigkeitsgesetz: Durch geeignete Wahl 
% der Modell-Parameter können alle Geschwindigkeitsgesetze der Gleichungen (4.2), (4.8) und (4.17) abgebildet werden.

r_1 = MP.k_1 * c_i(1)^MP.n * c_i(2)^MP.m; % [mol/(s*l)]
r_2 = MP.k_2 * c_i(3)^MP.o; % [mol/(s*l)]

% Umsatzrate der Komponenten
r_i = MP.nu_i_1 * r_1 + MP.nu_i_2 * r_2; % [mol/(s*l)] Entsrpicht den Gleichungen (7.1) ... (7.4)

% 7.16 MP.nu_i_1 = [-1 -2 1] MP.nu_i_2 = [1 2 -1]

% Stoffstrombilanz gemäss Gl. (7.57)
% Wobei dV_RM = dl * (d_R/2)^2 * pi()
dF_i_dl =  r_i * 1000 * (MP.d_R / 2)^2 * pi; % [mol/(s*m)]

%% Übergabe des Funktionswertes

dy_dt = [dF_i_dl]';

end
