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
    % time_s = hours(time)*3600+minutes(time)*60+seconds(time);
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

i = 1;
Exp(i).X_A = (Exp(i).kappa(1)-Exp(i).kappa) / (Exp(i).kappa(1)-Exp(i).kappa(end));
X_A = Exp(i).X_A;
c_A = c_A_0 * (1 - X_A);
c_B = c_B_0 - c_A_0 * X_A;
c_C = c_C_0 + c_A_0 * X_A;
c_D = c_D_0 + c_A_0 * X_A;

%% Linearisierung 1. Ordnung
lin_A_1 = log(c_A); % [-]

% Lineare Kurvenanpassung -> Gerade y = p(1)*x + p(2)
p = polyfit(Exp(i).t, lin_A_1, 1);

% Geschwindigkeitskonstante 1. Ordnung:
k_1 = -p(1) % [1/s]

% Bestimmtheitsmass für die lineare Regression
R_square = corrcoef(Exp(i).t, lin_A_1).^2;
R_square(1, 2);

% Plot figure 1
figure; hold on; grid on; grid minor;
plot(Exp(i).t, lin_A_1, ':o', 'linewidth', 1);
hold on
plot(Exp(i).t, polyval(p, Exp(i).t), '-', 'linewidth', 1);
xlabel('Zeit t [s]');
ylabel('$ln(c_A)$ [-]', 'interpreter', 'latex');

title('Linearisierte Messdaten aus [2]: 1. Ordnung');
legend('$ln(c_A)$', 'Regressionsgerade', 'interpreter', 'latex');
text(200, -5.5, ['R^2 = ', num2str(R_square(1, 2)), newline, ...
    'k_1 = ', num2str(k_1), ' [1/s]']);


%% Linearisierung 2. Ordnung c_A_0 = c_B_0
lin_A_2 = 1 ./ c_A; % [-]

% Lineare Kurvenanpassung -> Gerade y = p(1)*x + p(2)
p = polyfit(Exp(i).t, lin_A_2, 1);

% Geschwindigkeitskonstante 2. Ordnung:
k_2 = p(1) % [L/(mol*s)]

% Bestimmtheitsmass für die lineare Regression
R_square = corrcoef(Exp(i).t, lin_A_2).^2;
R_square(1, 2);

% Plot figure 2
figure; hold on; grid on; grid minor;
plot(Exp(i).t, lin_A_2, ':o', 'linewidth', 1);
hold on
plot(Exp(i).t, polyval(p, Exp(i).t), '-', 'linewidth', 1);
xlabel('Zeit t [s]');
ylabel('$1/c_A$ [-]', 'interpreter', 'latex');

title('Linearisierte Messdaten aus [2]: 2. Ordnung c_{A,0} = c_{B,0}');
legend('$1/c_A$', 'Regressionsgerade', 'interpreter', 'latex');
text(200, 2 * 10^4, ['R^2 = ', num2str(R_square(1, 2)), newline, ...
    'k_2 = ', num2str(k_2), ' [l/(mol*s)]']);


%% Linearisierung 2. Ordnung c_A != c_B
lin_A_2AB = log(c_B_0*c_A./(c_A_0*(c_B_0-c_A_0 + c_A))); % [-]

% Lineare Kurvenanpassung -> Gerade y = p(1)*x + p(2)
p = polyfit(Exp(i).t, lin_A_2AB, 1);

% Geschwindigkeitskonstante 2. Ordnung:
k_2AB = p(1)/(c_A_0-c_B_0) % [L/(mol*s)]

% Bestimmtheitsmass für die lineare Regression
R_square = corrcoef(Exp(i).t, lin_A_2AB).^2;
R_square(1, 2);

% Plot figure 3
figure; hold on; grid on; grid minor;
plot(Exp(i).t, lin_A_2AB, ':o', 'linewidth', 1);
hold on
plot(Exp(i).t, polyval(p, Exp(i).t), '-', 'linewidth', 1);
xlabel('Zeit t [s]');
ylabel('[-]', 'interpreter', 'latex');

title('Linearisierte Messdaten aus [2]: 2. Ordnung c_{A,0} != c_{B,0}');
legend('$ln(\frac{c_{B,0} \cdot c_A}{c_{A,0} \cdot (c_{B,0}-c_{A,0}+c_A)})$', ...
    'Regressionsgerade', 'interpreter', 'latex');
text(200, -1.5, ['R^2 = ', num2str(R_square(1, 2)), newline, ...
    'k_{2AB} = ', num2str(k_2AB), ' [l/(mol*s)]']);











