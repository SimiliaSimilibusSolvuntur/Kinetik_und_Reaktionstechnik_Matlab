clearvars; clc; 
close all;
%% Physikalische Konstanten
R = 8.314;       % [J/(mol*K)] ideale Gaskonstante

%%
%  1 A -> 1 B
a = 1;
b = 1;

c_A_in_list = [2, 1.2, 2, 1, 0.48, 1, 0.48, 0.48]; % [mol/l]
c_A_out_list = [1, 0.8, 0.65, 0.56, 0.42, 0.37, 0.28, 0.2]; % [mol/l]
t_list = [50, 16, 60, 22, 4.8, 72, 40, 112]; % [s]

% berchnung des Umsatzes X_A gemäss (4.23)
X_A_list = (c_A_in_list-c_A_out_list)./c_A_in_list;

% r_A ergibt sich nach Umstellen der Gleichung (7.53)
r_A2 = -(c_A_in_list .* X_A_list) ./ t_list;

% r_A ergibt sich nach Umstellen der Gleichung (7.54)
r_A = -(c_A_in_list - c_A_out_list) ./ t_list;

inv_r_A = 1./r_A;           % [(l*s)/mol] Inverse Umsatzgewschindigkeit

figure; hold on; grid on;
plot(c_A_out_list, -inv_r_A, 'o-');
for i = 1:8-1
    X = [c_A_out_list(i), c_A_out_list(i+1), c_A_out_list(i+1), c_A_out_list(i)];
    Y = [0, 0, -inv_r_A(i+1), -inv_r_A(i)];
    patch('XData', X,'YData',Y, 'FaceColor','blue','FaceAlpha',.1);
end
xlabel('$c_A$ [mol/l]', 'interpreter', 'latex');
ylabel('$-\frac{1}{r_A}$ [(l$\cdot$s)/mol]','interpreter', 'latex');
title('Gl. (7.62) - analog zu Gl. (7.20)', 'interpreter', 'latex');
XLim = get(gca, 'XLim');
set(gca, 'XLim', [0 XLim(2)]);

% Eintrittskonzentration [mol/l]
c_A_in = 0.8; % [mol/l]
c_B_in = 0; % [mol/l]

% Volumenstrom
V_in = 1; % [l/s]
F_A_in = V_in * c_A_in; % [mol/s]

% Gewünschter Umsatz bezüglich der Komponente A am Ausgang des Reaktors
X_A_f = 0.7; % [-] 

c_A_f = c_A_in * (1 - X_A_f);

%% Berechnung Fläche
c_A_in_list_cut = c_A_in_list(2:7); % [mol/l]
c_A_out_list_cut = c_A_out_list(2:7); % [mol/l]
t_list_cut = t_list(2:7); % [s]

% r_A ergibt sich nach Umstellen der Gleichung (7.54)
r_A_3 = -(c_A_in_list_cut - c_A_out_list_cut) ./ t_list_cut;

inv_r_A_3 = 1./r_A_3;           % [(l*s)/mol] Inverse Umsatzgewschindigkeit

% Einfügen von Punkt (0.024, 300) wenn X_A = 0.7 erreicht ist
inv_r_A_3 = [inv_r_A_3, -300];
c_A_out_list_2 = [0.8, 0.65, 0.56, 0.42, 0.37, 0.28, c_A_f]; % [mol/l]

% Berechnen der Fläche unter der Kurve -1/r_A vs. c_A anhand der Trapez-Methode
% Flaeche_2 = sum((c_A_in-c_A_f)/(n_X-1)*-mean(-[inv_r_A(1:end-1); inv_r_A(2:end)]))
Flaeche_3 = trapz(c_A_out_list_2, -inv_r_A_3);

% Berechnung der mittleren Verweilzeit gemäss (7.62)
tau_3 = -Flaeche_3; %[s]
disp(['tau_PFR  = ', num2str(tau_3, '%.4g'), ' s']);

% V_RM gemäss (7.32)
V_RM_PFR = V_in * tau_3;
disp(['V_RM_PFR  = ', num2str(V_RM_PFR, '%.4g'), ' l']);

%% CSTR 
tau_2 = (c_A_in - c_A_f) * 300;
disp(['tau_CSTR  = ', num2str(tau_2/60, '%.4g'), ' min']);

% V_RM gemäss (7.32)
V_RM_CSTR = V_in * tau_2;
disp(['V_RM_CSTR  = ', num2str(V_RM_CSTR, '%.4g'), ' l']);


