clear; clc; close all;

Fs = 48;

%% Wywołanie filtrów
Hd2 = filtr2();
Hd8 = filtr8();

[b2, a2] = tf(Hd2);
[b8, a8] = tf(Hd8);
%% Charakterystyka amplitudowa
[H2, f] = freqz(b2, a2, 2048, Fs);
[H8, ~] = freqz(b8, a8, 2048, Fs);

mag2 = 20*log10(abs(H2));
mag8 = 20*log10(abs(H8));

figure;
plot(f, mag2, 'b', 'DisplayName', 'Filtr II rzędu'); hold on;
plot(f, mag8, 'r', 'DisplayName','Filtr VIII rzędu');

grid on;
xlabel('Częstotliwość [Hz]');
ylabel('Amplituda [dB]');
title('Charakterystyka amplitudowa (porównanie)');
legend show;


%% Linie odniesienia
yline(-3, '--k', '-3 dB', 'LineWidth', 1.2, 'DisplayName', 'linia -3dB');
yline(-65, '--b', '-65 dB', 'LineWidth', 1.2, 'DisplayName', 'linia -65dB');

%% Punkt -3 dB
idx2_3 = find(mag2 <= -3, 1);
idx8_3 = find(mag8 <= -3, 1);


%% Punkt -65 dB
idx2_65 = find(mag2 <= -65, 1);
idx8_65 = find(mag8 <= -65, 1);


if ~isempty(idx2_3)
    plot(f(idx2_3), mag2(idx2_3), 'bo', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName','Częstotliwość odcięcia Filtra II rzędu');
    text(f(idx2_3), mag2(idx2_3)+10, sprintf('  %.2f Hz', f(idx2_3)));
end

if ~isempty(idx8_3)
    plot(f(idx8_3), mag8(idx8_3), 'ro', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName','Częstotliwość odcięcia Filtra VIII rzędu');
    text(f(idx8_3), mag8(idx8_3)+10, sprintf('  %.2f Hz', f(idx8_3)));
end

%% Wypisanie wyników
fprintf('\n=== WYNIKI ===\n');

fprintf('Filtr rzędu 2:\n');
fprintf('  fc (-3dB): %.2f Hz\n', f(idx2_3));

fprintf('\nFiltr rzędu 8:\n');
fprintf('  fc (-3dB): %.2f Hz\n', f(idx8_3));

if ~isempty(idx2_65)
    plot(f(idx2_65), mag2(idx2_65), 'ys', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Częstotliwość graniczna pasma zaporowego dla filtra II rzędu');
    text(f(idx2_65), mag2(idx2_65)+10, sprintf('  %.2f Hz', f(idx2_65)));
end

if ~isempty(idx8_65)
    plot(f(idx8_65), mag8(idx8_65), 'gs', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Częstotliwość graniczna pasma zaporowego dla filtra VIII rzędu');
    text(f(idx8_65), mag8(idx8_65)+10, sprintf('  %.2f Hz', f(idx8_65)));
end

hold off
%% Nachylenie zbocza
f1 = 6; f2 = 12;

[~, i1_2] = min(abs(f - f1));
[~, i2_2] = min(abs(f - f2));

[~, i1_8] = min(abs(f - f1));
[~, i2_8] = min(abs(f - f2));

slope2 = (mag2(i2_2) - mag2(i1_2)) / (log10(f2) - log10(f1));
slope8 = (mag8(i2_8) - mag8(i1_8)) / (log10(f2) - log10(f1));

fprintf('\nNachylenie:\n');
fprintf('  Rząd 2: %.2f dB/dekadę\n', slope2);
fprintf('  Rząd 8: %.2f dB/dekadę\n', slope8);

%% Opóźnienie grupowe
[gd2, fg] = grpdelay(b2, a2, 2048, Fs);
[gd8, ~] = grpdelay(b8, a8, 2048, Fs);

figure;
plot(fg, gd2, 'b', 'LineWidth', 1.5); hold on;
plot(fg, gd8, 'r', 'LineWidth', 1.5);
grid on;
xlabel('Częstotliwość [Hz]');
ylabel('Opóźnienie [próbki]');
title('Opóźnienie grupowe');
legend('Rząd 2', 'Rząd 8');

%% Maksymalne opóźnienie
fprintf('\nOpóźnienie grupowe (max):\n');
fprintf('  Rząd 2: %.2f próbek\n', max(gd2));
fprintf('  Rząd 8: %.2f próbek\n', max(gd8));

%%  Odpowiedź skokowa
figure;
stepz(b2, a2); hold on;
stepz(b8, a8);
title('Odpowiedź skokowa');
legend('Rząd 2', 'Rząd 8');
grid on;

%% Bieguny i zera
figure;
subplot(1,2,1);
zplane(b2, a2); title('Rząd 2');

subplot(1,2,2);
zplane(b8, a8); title('Rząd 8');

%% Charakterystyka fazowa (w radianach)
% Obliczanie fazy
phi2 = angle(H2);
phi8 = angle(H8);

% Odwijanie fazy (usuwanie skoków o 2*pi)
phi2_unwrapped = unwrap(phi2);
phi8_unwrapped = unwrap(phi8);

% Ograniczamy zakres do 99% Nyquista, żeby uciąć błąd na końcu
figure;
idx_limit = round(length(f) * 0.99); 
plot(f(1:idx_limit), phi2_unwrapped(1:idx_limit), 'b', 'DisplayName', 'Filtr II rząd '); hold on;
plot(f(1:idx_limit), phi8_unwrapped(1:idx_limit), 'r', 'DisplayName', 'Filtr VIII rząd');

% Zaznaczenie fazy w punktach odcięcia (-3dB)
if ~isempty(idx2_3)
    plot(f(idx2_3), phi2_unwrapped(idx2_3), 'o', 'MarkerFaceColor', 'b', 'DisplayName', 'Punkt -3dB dla filtra II rzędu');
    text(f(idx2_3), phi2_unwrapped(idx2_3)-0.5, sprintf(' %.2f rad', phi2_unwrapped(idx2_3)));
end

if ~isempty(idx8_3)
    plot(f(idx8_3), phi8_unwrapped(idx8_3), 'o', 'MarkerFaceColor', 'r',  'DisplayName', 'Punkt -3dB dla filtra VIII rzędu');
    text(f(idx8_3), phi8_unwrapped(idx8_3)-0.5, sprintf(' %.2f rad', phi8_unwrapped(idx8_3)));
end

grid on;
xlabel('Frequency [Hz]');
ylabel('Phaze [rad]');
title('Phaze Response');
legend show;

% Wypisanie wyników do konsoli
fprintf('\nPrzesunięcie fazowe przy fc (-3dB):\n');
fprintf('  Rząd 2: %.2f rad\n', phi2_unwrapped(idx2_3));
fprintf('  Rząd 8: %.2f rad\n', phi8_unwrapped(idx8_3));

%% Charakterystyka opóźnienia fazowego (Phase Delay) - z wycięciem Nyquista
% Obliczanie opóźnienia fazowego w próbkach
w_vec = 2*pi*f/Fs; 
pd2 = -phi2_unwrapped ./ w_vec;
pd8 = -phi8_unwrapped ./ w_vec;

% Poprawka dla DC (częstotliwość 0), gdzie faza=0
pd2(1) = pd2(2); 
pd8(1) = pd8(2);

% --- WYCINANIE ARTEFAKTÓW PRZY NYQUISTCIE ---
% Wycinamy ostatnie 1-2% zakresu, gdzie obliczenia numeryczne "wariują"
limit_idx = round(length(f) * 0.98); 
f_plot = f(1:limit_idx);
pd2_plot = pd2(1:limit_idx);
pd8_plot = pd8(1:limit_idx);

figure;
plot(f_plot, pd2_plot, 'b',  'DisplayName', 'Filtr II rząd'); hold on;
plot(f_plot, pd8_plot, 'r',  'DisplayName', 'Filtr VIII rząd');

% Zaznaczenie punktów dla fc (tylko jeśli mieszczą się w zakresie plotu)
if ~isempty(idx2_3) && idx2_3 <= limit_idx
    plot(f(idx2_3), pd2(idx2_3), 'bo', 'MarkerFaceColor', 'b', 'DisplayName','Opóźnienie dla częstotliwości odcięcia filtra II rzędu');
    text(f(idx2_3), pd2(idx2_3)+0.5, sprintf(' %.2f', pd2(idx2_3)));
end
if ~isempty(idx8_3) && idx8_3 <= limit_idx
    plot(f(idx8_3), pd8(idx8_3), 'ro', 'MarkerFaceColor', 'r', 'DisplayName','Opóźnienie dla częstotliwości odcięcia filtra II rzędu');
    text(f(idx8_3), pd8(idx8_3)+0.5, sprintf(' %.2f', pd8(idx8_3)));
end
grid on;
title("Phaze Delay");
ylabel("Delay [sample]");
xlabel("Frequency [Hz]");
fprintf('\nOpóźnienie fazowe przy fc (-3dB):\n');
fprintf('  Rząd 2: %.2f próbek\n', pd2(idx2_3));
fprintf('  Rząd 8: %.2f próbek\n', pd8(idx8_3));

%% Charakterystyka opóźnienia grupowego (Group Delay)
[gd2, fg] = grpdelay(b2, a2, 2048, Fs);
[gd8, ~] = grpdelay(b8, a8, 2048, Fs);

% Wycinanie końcówki (Nyquista), aby uniknąć artefaktów
limit_idx_gd = round(length(fg) * 0.98);
fg_plot = fg(1:limit_idx_gd);
gd2_plot = gd2(1:limit_idx_gd);
gd8_plot = gd8(1:limit_idx_gd);

figure;
plot(fg_plot, gd2_plot, 'b', 'DisplayName', 'Filtr II rząd'); hold on;
plot(fg_plot, gd8_plot, 'r', 'DisplayName', 'Filtr VIII rząd');

% Zaznaczenie maksymalnego opóźnienia dla rzędu 8
[max_gd8, max_idx8] = max(gd8_plot);
plot(fg_plot(max_idx8), max_gd8, 'ro', 'MarkerFaceColor', 'r', 'DisplayName','maksymaplne opóźnienie grupowe dla filtrza VIII rzędu');
text(fg_plot(max_idx8), max_gd8 + 0.5, sprintf(' Max: %.2f pr.', max_gd8));

% Zaznaczenie maksymalnego opóźnienia dla rzędu 2
[max_gd2, max_idx2] = max(gd2_plot);
plot(fg_plot(max_idx2), max_gd2, 'bo', 'MarkerFaceColor', 'b','DisplayName','maksymaplne opóźnienie grupowe dla filtrza II rzędu');
text(fg_plot(max_idx2), max_gd2 + 0.5, sprintf(' Max: %.2f pr.', max_gd2));

grid on;
xlabel('Frequency [Hz]');
ylabel('Delay [sample]');
title('Group Delay');
legend show;

fprintf('\nMaksymalne opóźnienie grupowe:\n');
fprintf('  Rząd 2: %.2f próbek\n', max(gd2_plot));
fprintf('  Rząd 8: %.2f próbek\n', max_gd8);
