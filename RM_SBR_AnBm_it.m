function dy_dt = RM_SBR_AnBm_it(t, y, MP)

% Reaktor-Model für einen idealen, isothermen Semi-Batch-Reaktor
% Kann auch zur Simulation eines Batch-Reaktors verwendet werden
%
% Es können folgende Reaktionstypen abgebildet werden:
%   nu_1*A + nu_2*B -> nu_4*C + nu_4*D
%   nu_1*A -> nu_2*B + nu_4*C + nu_4*D
%
% Lösungsmittel = Komponente E (Nr. 5). Es gilt nu_5 = 0.
%
% Funktionsargumente
% t           [s] Zeitpunkt "t" zwingend das erste Funktionsargument
% y           Vektor mit den Zeitabhängigen Variablen
%             Element 1 ... 5: [mol] Molmenge der Komponenten A ... D + Lsm
%             Element 6: [kg] Reaktionsmasse
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
% MP.rho_RM      [kg/m3] Dichte der Reaktionsmasse
% MP.t_dos_start [s] Startzeit der Dosierung
% MP.t_dos       [s] Dauer der Dosierung
% MP.F_dos_i     [mol/s] Vektor mit Stoffströmen während der Dosierung

%% Momentan-Werte für die zeitabhängigen Variablen im Vektor "y"

% Wir setzen alle Momentan-Werte im Spaltenvektor "y" auf >= 0
y(y<0) = 0;

% Auslesen des Spaltenvektors "y" mit allen Momentan-Werten
n_i = y(1:end-1)';  % [mol] Molmenge der Komponenten A ... D + Lsm
m_RM = y(end);      % [kg] Reaktionsmasse

%% Berechnung des Reaktionsvolumens und der Konzentrationen

V_RM = m_RM / MP.rho_RM * 1000; % [l] Reaktionsvolumen gemäss Gleichung (5.30)
c_i = n_i / V_RM; % [mol/l] Konzentrationen der einzelnen Komponenten

%% Stoffströme und Stoffstrombilanz

% Reaktionsgeschwindigkeit bzw. Geschwindigkeitsgesetz: Durch geeignete Wahl
% der Modell-Parameter können alle Geschwindigkeitsgesetze der
% Gleichungen (4.2), (4.8) und (4.17) abgebildet werden.
r = MP.k * c_i(1)^MP.n * c_i(2)^MP.m; % [mol/(s*l)]

% Umsatzrate der Komponenten
r_i = MP.nu_i * r;    % [mol/(s*l)] Entsrpicht den Gleichungen (7.1) ... (7.4)

% Berechnung von F_i_in für den Semi-Batch-Reaktor
if (t >= MP.t_dos_start) && (t <= MP.t_dos_start + MP.t_dos)
    F_i_in = MP.F_dos_i; % [mol/s]
else
    F_i_in = 0; % [mol/s]
end

% Stoffstrombilanz gemäss Gl. (7.26) bzw. (7.27)
dn_i_dt =  F_i_in + r_i * V_RM;  % [mol/s]

%% Massenströme und Massenstrombilanz

% Berechnung von mf_i_in gemäss Gl. (7.29)
mf_i_in = F_i_in .* MP.mw_i/1000; % [kg/s]

% Massenflussbilanz gemäss Gl. (7.29)
dm_RM_dt = sum(mf_i_in);  % [kg/s]

%% Übergabe des Funktionswertes = Ableitungen aller Variablen im Vektor "y" nach der Zeit

dy_dt = [dn_i_dt, dm_RM_dt]';


end




















