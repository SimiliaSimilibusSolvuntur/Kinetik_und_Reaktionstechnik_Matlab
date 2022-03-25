function dy_dt = RM_CSTR_AnBm_it(t, y, MP)

% Reaktor-Model für einen isothermen, idealen CSTR-Reaktor
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
% MP.F_i_in      [mol/s] Vektort mit Stoffströmen der einzelnen Komponenten in den Reaktor
% MP.rho_in      [kg/m3] Dichte des Feeds in den Reaktor
% MP.V_R         [l] Volumen des Reaktors

%% Momentan-Werte für die zeitabhängigen Variablen

% Wir setzen alle Momentan-Werte im Spaltenvektor "y" auf >= 0

y(y<0) =0;

% Auslesen des Spaltenvektors "y" mit allen Momentan-Werten
n_i = y(1:end-2)';  % [mol] Molmenge der Komponenten A ... D + Lsm
m_RM = y(end-1);    % [kg] Reaktionsmasse
V_RM = y(end);      % [l] Reaktionsvolumen

%% Berechnung der Dichte der Reaktionsmasse und der Konzentrationen

rho_RM = m_RM / V_RM * 1000; % [kg/m3] Dichte der Reaktionsmasse Gl. (7.39)
if isnan(rho_RM) | isinf(rho_RM)
    rho_RM = 1000;
end
c_i = n_i / V_RM;   % [mol/l] Konzentrationen der einzelnen Komponenten
c_i(isnan(c_i))=0;
c_i(isinf(c_i))=0;
m_i = n_i .* MP.mw_i / 1000;   % [kg] Menge der Komponente i
w_i = m_i ./ sum(m_i);      % [kg/kg] Gewichtsprozent der Komponente i
w_i(isnan(w_i))=0;
w_i(isinf(w_i))=0;

%% Gesamtmassenstrom in den Reaktor

% Berechnung von mf_i_in gemäss Gl. (7.41)

mf_i_in = MP.F_i_in .* MP.mw_i/1000; % [kg/s]
mf_in = sum(mf_i_in); % [kg/s] Gesamtmassenstrom in den Reaktor

%% Volumenströme und Volumenstrombilanz

Vf_in = mf_in / MP.rho_in * 1000; % [l/s] Gl. (7.43)

% Festlegen des Betriebszustandes.
if V_RM > MP.V_R
    % Überlaufen der Reaktionsmasse bei konstantem Volumen
    Vf_RM_out = Vf_in; % [l/s] Gl. (7.47)
    dV_RM_dt = 0; % [l/s] Gl. (7.46)
else
    % Füllen des Reaktors, kein Überlauf.
    Vf_RM_out = 0; % [l/s] Gl. (7.44)
    dV_RM_dt = Vf_in; % [l/s] Gl. (7.43)
end

%% Stoffströme und Stoffstrombilanz

% Reaktionsgeschwindigkeit bzw. Geschwindigkeitsgesetz: Durch geeignete Wahl
% der Modell-Parameter können alle Geschwindigkeitsgesetze der
% Gleichungen (4.2), (4.8) und (4.17) abgebildet werden.
r = MP.k*c_i(1)^MP.n*c_i(2)^MP.m; % [mol/(s*l)]

% Umsatzrate der Komponenten
r_i = MP.nu_i * r;    % [mol/(s*l)] Entsrpicht den Gleichungen (7.1) ... (7.4)

% Stoffstrom aus dem Reaktor gemäss Gleichung (7.38)
F_i_out = w_i .* Vf_RM_out .* rho_RM ./ MP.mw_i; % [mol/s]

% Stoffstrombilanz gemäss Gl. (7.33) bzw. (7.37)
dn_i_dt =  MP.F_i_in - F_i_out + r_i*V_RM;  % [mol/s]

%% Gesamtmassenstrom aus dem Reaktor und Massenstrombilanz

% Berechnung von mf_in_i gemäss Gl. (7.41)

mf_i_out = F_i_out .* MP.mw_i/1000; % [kg/s]
mf_out = sum(mf_i_out);          % [kg/s] Gesamtmassenstrom aus dem Reaktor

% Massenflussbilanz gemäss Gl. (7.41)
dm_RM_dt = mf_in - mf_out;  % [kg/s]

%% Übergabe des Funktionswertes

dy_dt = [dn_i_dt, dm_RM_dt, dV_RM_dt]';

end
