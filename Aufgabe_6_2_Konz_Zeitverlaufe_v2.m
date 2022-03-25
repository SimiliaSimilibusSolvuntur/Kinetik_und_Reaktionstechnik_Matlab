clearvars; close all; clc

% Zellkonstante der Leitfähigkeitsmesszelle K = l/A [1/cm]
K_cell = 1; % [1/cm]

% Einlesen der Messdaten aus drei Datenfiles
file_names = {'AcOEt_25.txt'};
delimiter = ';';  % \t für Tabulator
Kopfzeilen = 6 + [56, 23, 47]; % Kürzen der Daten, am Anfang wird ein konstanter Wert gemessen, da noch nichts passiert
% Datenformat für jede einzelne Spalte
%   double: %f
%	datum:  %{....}D
%   zeit:   %{....}T
%   Zeilenabschluss: %[^\n\r]
formatSpec = '%f %{hh:mm:ss.SSS}T %f %f %f';
for i = 1:length(file_names)
    fileID = fopen(char(file_names(i)),'r');
    dataArray = textscan(fileID, formatSpec, 'Headerlines', Kopfzeilen(i), ...
        'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
    fclose(fileID);
    time = dataArray{:, 2};
    % falsch: time_s = hours(time)*3600+minutes(time)*60+seconds(time);
    time_s = seconds(time); %[s]
    Exp(i).t = time_s - time_s(1);         % [s]
    Exp(i).kappa = dataArray{:, 4}*K_cell; % Spezifische Leitfähigkeit [mS/cm]
    Exp(i).T = dataArray{:, 5};            % Temperatur T [°C]
end
clearvars dataArray;

c_A_0 = 0.04; % [mol/L]
c_B_0 = 0.05; % [mol/L]
c_C_0 = 0; % [mol/L]
c_D_0 = 0; % [mol/L]
% t = linspace(0, 2500, 100);  %[s] Zeit t als Zeilenvektor

i = 1;
Exp(i).X_A = (Exp(i).kappa(1)-Exp(i).kappa) / (Exp(i).kappa(1)-Exp(i).kappa(end));
X_A = Exp(i).X_A;
c_A = c_A_0 * (1 - X_A);
c_B = c_B_0 - c_A_0 * X_A;
c_C = c_C_0 + c_A_0 * X_A;
c_D = c_D_0 + c_A_0 * X_A;

% Graph plot
figure; hold on; grid on; grid minor;
title('Aufgabe 6.2 bei 25°');
xlabel('Zeit [s]');
set(gca, 'XLim', [0 900]);

yyaxis right;
plot(Exp(i).t, X_A);
ylabel('Umsetzung X_A [-]');
set(gca, 'YLim', [0 1]);

yyaxis left;
plot(Exp(i).t, c_A, Exp(i).t, c_B, Exp(i).t, c_C, Exp(i).t, c_D);
ylabel('Konzentration [mol/l]');
set(gca, 'YLim', [0 0.1]);

legend('$c_A$', '$c_B$', '$c_C$', '$c_D$',...
    '$X_A$', 'interpreter', 'latex');




