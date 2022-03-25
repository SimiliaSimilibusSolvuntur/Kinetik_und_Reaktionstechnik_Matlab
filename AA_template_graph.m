%% erstellen eines Graphen
figure; hold on; grid on; grid minor;
title('Titel Graph');
plot(t, A, 'o'); % oder farbe 'r' = rot
xlabel('Titel Achse');
ylabel('Titel Achse');
legend('Legende mit Latex-Interpreter, z B: $cm^{-1}$ oder $c_A$', ...
    'interpreter', 'latex');

% plot(t_cut, lin_A_2AB, ':o', 'linewidth', 1);

%% erstellen eines Graphen
figure; hold on; grid on;
title('Titel Graph');
plot(t/60, A, 'b', ...
    t/60, B, 'r', ...
    t/60, C, 'k');
xlabel('Zeit [min]');
ylabel('Konzentration [mol/l]');
set(gca, 'YLim', [0 2]); % Limitierung der Achse

legend('Legende mit Latex-Interpreter, z B: $c_A$', ...
    'interpreter', 'latex');

%% erstellen eines Graphen mit zwei Y-Achsen
figure; hold on; grid on; grid minor;
title('Titel Graph');

xlabel('Titel Achse');
set(gca, 'XLim', [0 900]); % Limitierung der Achse

yyaxis right;
plot(t, A, 'o');
ylabel('Titel Achse');
set(gca, 'YLim', [0 1]); % Limitierung der Achse

yyaxis left;
plot(t, A, t, B, t, C, t, D);
ylabel('Titel Achse');
set(gca, 'YLim', [0 0.1]); % Limitierung der Achse

legend('Legende mit Latex-Interpreter, z B: $c_A$', ...
    'interpreter', 'latex');

