function dy_dt = RM_SBR_2AB_it(t, y, mw_i, k_2AB, rho_RM, t_dos_start, t_dos, F_dos_i)
% Reaktor-Model f�r einen isothermen Semi-Batch-Reaktor
% Kann auch zur Simulation eines Batch-Reaktors verwendet werden
% Gem�ss S. 47 Fallstudie Ethylacetat, Vorlesung Reaktionstechnik

% Es kann nur eine Reaktion abgebildet werden:
%   a*A + b*B -> c*C + d*D

% St�chiometrie: Zeile 36.
% Geschwindigkeitsgesetz: Zeile 44.
% L�sungsmittel = Komponente E (Nr. 5)

% Funktionsargumente 
% t           [s] Zeitpunkt "t" zwingend das erste Funktionsargument
% y           Vektor mit den Zeitabh�ngigen Variablen
%             Element 1 ... 5: [mol] Molmenge der Komponenten A ... D + Lsm
%             Element 6: [kg] Reaktionsmasse
% mw_i        [g/mol] Molmassen der einzelnen Komponenten
% k_2AB       [l/(mol*s)] Geschwindigkeitskonstante 
% rho_RM      [kg/m3] Dichte der Reaktionsmasse 
% t_dos_start [s] Startzeit der Dosierung
% t_dos       [s] Dauer der Dosierung
% F_dos_i     [mol/s] Stoffstr�me w�hrend der Dosierung 

% A. Zogg 01.04.2020

%% Momentan-Werte f�r die zeitabh�ngigen Variablen
% Wir setzen alle Momentan-Werte im Spaltenvektor "y" auf >= 0
y(y<0) =0;

% Auslesen des Spaltenvektors "y" mit allen Momentan-Werten
n_i = y(1:end-1)';  % [mol] Molmenge der Komponenten A ... D + Lsm
m_RM = y(end);     % [kg] Reaktionsmasse

%% Modell-Parameter
nu_i = [-1, -1, 1, 2, 0]; % [-] St�chiometrische Koeffizienten der Komponenten 
  
%% Berechnung des Reaktionsvolumen und der Konzentrationen
V_RM = m_RM /rho_RM * 1000; % [l] Reaktionsvolumen gem�ss Gleichung (5.30)
c_i = n_i / V_RM; % [mol/l] Konzentrationen der einzelnen Komponenten

%% Stoffstr�me und Stoffstrombilanz
% Reaktionsgeschwindigkeit
r = k_2AB*c_i(1)*c_i(2); % [mol/(s*l)] Gleichung (2.13)
% Umsatzrate der Komponente A 
r_i = nu_i * r;    % [mol/(s*l)] Gleichung (5.1) ... (5.4)

% Berechnung von F_i_in f�r den Semi-Batch-Reaktor
if (t >= t_dos_start) && (t <= t_dos_start+t_dos)
    F_i_in = F_dos_i; % [mol/s]
else
    F_i_in = 0; % [mol/s]
end

% Stoffstrombilanz gem�ss Gl. (5.26) bzw. (5.27)
dn_i_dt =  F_i_in + r_i*V_RM;  % [mol/s]

%%  Massenstr�me und Massenstrombilanz
% Berechnung von mf_i_in gem�ss Gl. (5.29)
mf_i_in = F_i_in .* mw_i/1000; % [kg/s]

% Massenflussbilanz gem�ss Gl. (5.29)
dm_RM_dt = sum(mf_i_in);  % [kg/s]

%% �bergabe des Funktionswertes
dy_dt = [dn_i_dt, dm_RM_dt]';

end



%%


