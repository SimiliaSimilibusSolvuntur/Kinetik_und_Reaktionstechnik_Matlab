function dy_dt = RM_CSTR_2AB_it(t, y, mw_i, k_2AB, F_i_in, rho_in, V_R)
% Reaktor-Model für einen isothermen, idealen CSTR-Reaktor

% Es kann nur eine Reaktion abgebildet werden:
%   a*A + b*B -> c*C + d*D

% Anpassen der Stöchiometrie: Zeile xy.
% Anpassen des Geschwindigkeitsgesetzes: Zeile xy.
% Lösungsmittel = Komponente E (Nr. 5)

% Funktionsargumente 
% t           [s] Zeitpunkt "t" zwingend das erste Funktionsargument
% y           Vektor mit den Zeitabhängigen Variablen
%             Element 1 ... 5: [mol] Molmenge der Komponenten A ... D + Lsm
%             Element 6: [kg] Reaktionsmasse
%             Element 7: [l] Reaktionsvolumen
% mw_i        [g/mol] Molmassen der einzelnen Komponenten
% k_2AB       [l/(mol*s)] Geschwindigkeitskonstante 
% F_i_in      [mol/s] Stoffströme der einzelnen Komponenten in den Reaktor
% rho_in      [kg/m3] Dichte des Feeds in den Reaktor 

% A. Zogg 09.04.2020

%% Momentan-Werte für die zeitabhängigen Variablen
% Wir setzen alle Momentan-Werte im Spaltenvektor "y" auf >= 0
y(y<0) =0;

% Auslesen des Spaltenvektors "y" mit allen Momentan-Werten
n_i = y(1:end-2)';  % [mol] Molmenge der Komponenten A ... D + Lsm
m_RM = y(end-1);   % [kg] Reaktionsmasse
V_RM = y(end);     % [l] Reaktionsvolumen

%% Modell-Parameter
nu_i = [-1, -1, 1, 2, 0]; % [-] Stöchiometrische Koeffizienten der Komponenten 
  
%% Berechnung der Dichte der Reaktionsmasse und der Konzentrationen
rho_RM = m_RM / V_RM * 1000; % [kg/m3] Dichte der Reaktionsmasse Gl. (5.39)
if isnan(rho_RM) | isinf(rho_RM)
    rho_RM = 1000;
end
c_i = n_i / V_RM;   % [mol/l] Konzentrationen der einzelnen Komponenten
c_i(isnan(c_i))=0;
c_i(isinf(c_i))=0;
m_i = n_i .* mw_i / 1000;   % [kg] Menge der Komponente i
w_i = m_i ./ sum(m_i);      % [kg/kg] Gewichtsprozent der Komponente i
w_i(isnan(w_i))=0;
w_i(isinf(w_i))=0;

%% Gesamtmassenstrom in den Reaktor
% Berechnung von mf_i_in gemäss Gl. (5.41)
mf_i_in = F_i_in .* mw_i/1000; % [kg/s]
mf_in = sum(mf_i_in); % [kg/s] Gesamtmassenstrom in den Reaktor

%% Volumenströme und Volumenstrombilanz
Vf_in = mf_in / rho_in * 1000; % [l/s] Gl. (5.43)

% Festlegen des Betriebszustandes. 
if V_RM > V_R 
    % Überlaufen der Reaktionsmasse bei konstantem Volumen
    Vf_RM_out = Vf_in; % [l/s] Gl. (5.47)
    dV_RM_dt = 0; % [l/s] Gl. (5.46)
else
    % Füllen des Reaktors, kein Überlauf.
    Vf_RM_out = 0; % [l/s] Gl. (5.44) 
    dV_RM_dt = Vf_in; % [l/s] Gl. (5.43) 
end


%% Stoffströme und Stoffstrombilanz
% Reaktionsgeschwindigkeit 
r = k_2AB*c_i(1)*c_i(2); % [mol/(s*l)] Gleichung (2.13)
% Umsatzrate der Komponente i 
r_i = nu_i * r;    % [mol/(s*l)] Gleichung (5.1) ... (5.4)

% Stoffstrom aus dem Reaktor gemäss Gleichung (5.38)
F_i_out = w_i .* Vf_RM_out .* rho_RM ./ mw_i; % [mol/s]

% Stoffstrombilanz gemäss Gl. (5.33) bzw. (5.37)
dn_i_dt =  F_i_in - F_i_out + r_i*V_RM;  % [mol/s]

%%  Gesamtmassenstrom aus dem Reaktor und Massenstrombilanz
% Berechnung von mf_in_i gemäss Gl. (5.41)
mf_i_out = F_i_out .* mw_i/1000; % [kg/s]
mf_out = sum(mf_i_out);          % [kg/s] Gesamtmassenstrom aus dem Reaktor

% Massenflussbilanz gemäss Gl. (5.41)
dm_RM_dt = mf_in - mf_out;  % [kg/s]

%% Übergabe des Funktionswertes
dy_dt = [dn_i_dt, dm_RM_dt, dV_RM_dt]';

end

