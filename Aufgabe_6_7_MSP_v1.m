clearvars; close all; clc
% A Essigsäureanhydrid 0.25 mol 23.6 ml
% B 1-Butanol 0.375 mol 35.4 ml

% A + B -> C

%% a.)
% Es wird die Zudosierung nicht beachtet 
% sehr kurze Dosierzeit

%% b.) 

load('Absorptionsverlauf.mat');

% erstellen eines Graphen
figure, hold on, grid on, grid minor;
title('Aufgabe 6.7 Experiment 1 bei 40°C');
plot(t, A, 'o');
xlabel('Zeit [s]');
ylabel('Absorption [-]');
legend('1701 $cm^{-1}$', 'interpreter', 'latex');

%% c.) Essigsäureanhydrid 


%% d.) 

c_A_0 = 0.25/(23.6+35.4)*1000 % [mol/l] Anfangskonzentration Essigsäureanhydrid
c_B_0 = 0.375/(23.6+35.4)*1000 % [mol/l] Anfangskonzentration 1-Butanol
c_C_0 = 0; % [mol] Anfangskonzentration

A_C_max = 0.4606 

a = 1;
b = 1;
c = 1;

X_A =  A / A_C_max;

figure; hold on; grid on; grid minor;
title('Umsatz für stöch. limitierend');
plot(t, X_A, 'o'); 
xlabel('Zeit [s]');
ylabel('X_A');
legend('$X_A$', ...
    'interpreter', 'latex');


%% e.) Konzentrations-Zeit-Kurve

% c_A = 1 - c_A_0 * (X_A);
c_A = c_A_0 * (1 - X_A)

figure; hold on; grid on; grid minor;
title('Konz-Zeit für stöch. limitierend');
plot(t, c_A, 'or'); 
xlabel('Zeit [s]');
ylabel('c_A');
legend('$c_A$', ...
    'interpreter', 'latex');

%% f.)

% Ausschneiden der Messdaten
c_A_cut = c_A(1:70);
t_cut = t(1:70);

% Linearisierung 2. Ordnung c_A != c_B
lin_A_2AB = log(c_B_0*c_A_cut./(c_A_0*(c_B_0-c_A_0 + c_A_cut))); % [-]

% Lineare Kurvenanpassung -> Gerade y = p(1)*x + p(2)
p = polyfit(t_cut, lin_A_2AB, 1);

% Geschwindigkeitskonstante 2. Ordnung:
k_2AB = p(1)/(c_A_0-c_B_0) % [L/(mol*s)]

% Bestimmtheitsmass für die lineare Regression
R_square = corrcoef(t_cut, lin_A_2AB).^2;
R_square(1, 2);

% Plot figure 3
figure; hold on; grid on; grid minor;
plot(t_cut, lin_A_2AB, ':o', 'linewidth', 1);
plot(t_cut, polyval(p, t_cut), '-', 'linewidth', 1);
xlabel('Zeit t [s]');
ylabel('[-]', 'interpreter', 'latex');

title('Linearisiert');
legend('$ln(\frac{c_{B,0} \cdot c_A}{c_{A,0} \cdot (c_{B,0}-c_{A,0}+c_A)})$', ...
    'Regressionsgerade', 'interpreter', 'latex');
text(1000, -0.4, ['R^2 = ', num2str(R_square(1, 2)), newline, ...
    'k_{2AB} = ', num2str(k_2AB), ' [l/(mol*s)]']);


%% g.)

R = 8.314510; % [J/(mol*K)] R = N_A * k_B
E_A = 100; % [kJ/mol]
T_1 = 273 + 40; % [K]
T_2 = 273 + 20; % [K]

Arr = k_2AB;
T = 273 + 20;
% Arrhenius
k = Arr * exp(-E_A/(R*T))

ln_k_2 = log(k_2AB)-(E_A/R)*(1/T_1-1/T_2);
k_2 = exp(ln_k_2)












