clc;
close all;


%=== 3. Sformułowanie problemu – „Dyskretny regulator Kesslera w
%serwomechanizmie DC ,czyli Perdita Durango znudzona na technikach komputerowych.
% numer grupy 4.

%% === zdefiniwanie modelu silnika 
CePhi = 2.62;
Tm = 0.18;

% Przyjęcie 𝛿1 = 0,0𝑋[𝑠] 𝛿2 = 0,0𝑌[𝑠] 𝛿3 = 0,0𝑍[𝑠] na podstawie
% readd reakcji
% na pojawiający się prąd na wejściu przetwornika a po czasie T4 dopiero
% prąd pojawi się na wyjściu przetwornika Węc 
% delta1 = 0.7μs + 2 μs = 2.7μs = 2.7*10^(-6)s
delta1 = 2.7*10^(-6);

% Człon przetwornika C/A
% Dobrane urządzenie DAC0808 8-Bit D/A Converter datasheet https://www.ti.com/lit/ds/symlink/dac0808.pdf
% do obliczenia stałej czasowej przyjęto t ( settling time ) + t ( propagation time delay )=
% = 150ns +100ns = 250 *10^-9s 
delta2 = 250*10^(-9);


%Człon filtra 
% Jako człon filtra (δ3​) przyjęto aktywny filtr dolnoprzepustowy zbudowany
% na wzmacniaczu operacyjnym TL082. datasheet https://www.ti.com/lit/ds/symlink/tl082-n.pdf
% Stała czasowa filtra wynika z doboru
% elementów R=50 kΩ oraz C=100 nF (δ3​=R⋅C=0,005 [s]). Filtr ten ma na celu
% eliminację tętnień napięcia prądnicy tachometrycznej wynikających z
% komutacji szczotkowej silnika oraz zapobieganie zjawisku aliasingu w
% cyfrowym torze pomiarowym 
delta3 = 0.005;

delta_sum = delta1+ delta2+delta3; % zsumowanie wszystkich członów układu
s = tf('s'); %zdefiniowanie zmiennej Laplace'a 
G_0 = (1/CePhi)/ ( (Tm*s+1) * (delta_sum*s +1) ); %zdfiniowanie transmitancji układu

gen_wykresy_uni({G_0}, {0}, 'Wykres modelu wyjściwego', {'G_01'}, 1.5);
gen_tabela_parametry({G_0}, {'G_0'},'tab 1. parametry silnika');
%% === Wyznaczenie regulatorów keslera w domecie ciągłej ==
% === Metoda modelu układu zamkniętego 1.1a
[UR1_ukl, Reg1_ukl] = gen_kessler('ukl', 'modul');

% === Metoda modelu układu zamkniętego 1.1a
[UR2_ukl, Reg2_ukl] = gen_kessler('ukl', 'sym');



gen_wykresy_uni({UR1_ukl, UR2_ukl}, {Reg1_ukl, Reg2_ukl},'Wykres metody układu zamkniętego', {'UR1 moduł optimum', 'UR2 symetrycznye optimum'}, 1);

% Metoda uproszczona przeliczeniowa na podstawie tabel
% a) kryterium optimum modułu 
[UR1_tab, Reg1_tab] = gen_kessler('tab', 'modul');

% b) kryterium symetryczne
[UR2_tab, Reg2_tab] = gen_kessler('tab', 'sym');

gen_wykresy_uni({UR1_tab, UR2_tab}, {Reg1_tab, Reg2_tab}, 'Wykres metody uproszczonej tablicowej', {'UR1 tab optimum_modułu','UR2 tab symetryczny'}, 1);

tab_ciagle = gen_tabela_parametry({UR1_ukl, UR2_ukl, UR1_tab, UR2_tab}, {'UR1_ukl', 'UR2_ukl', 'UR1_tab', 'UR2_tab'}, 'Tab 2. Paramtery dla regulatorów ciągłych');

%% === Przeprowadzenie dyskretyzacji Regulatorów ===

% Przyjęcie czasu próbkowania na podstawie najbardziej dominującej stałej
% czasowej którą jest stała mechaniczna Tm = 0.18s czyli Tp = 0.25Tm =
% 0.045

% === Metoda modelu układu zamkniętego 1.1
% a) b) Metoda optimum modułu (filtr drugiego rzędu ) 
[UR1_ukl_D, Reg1_ukl_D] = gen_kessler('ukl', 'modul', 0.001);



% === Metoda modelu układu zamkniętego 1.1 
% b) metoda symetryczna (filtr trzeciego rzędu)
[UR2_ukl_D, Reg2_ukl_D] = gen_kessler('ukl', 'sym', 0.001);


gen_wykresy_uni({UR1_ukl_D, UR2_ukl_D}, {Reg1_ukl_D, Reg2_ukl_D}, 'Wykres metody układu zamkniętego dyskretny', {'UR1 optimum modułu', 'UR2 symetrycznego potimum'}, 1);


% Metoda uproszczona przeliczeniowa na podstawie tabel
% a) kryterium optimum modułu 
[UR1_tab_D, Reg1_tab_D] = gen_kessler('tab', 'modul', 0.001);

% b) kryterium symetryczne
[UR2_tab_D, Reg2_tab_D] = gen_kessler('tab', 'sym', 0.001);


% === Wygenerowanie wykresów porównawczych ===
gen_wykresy_uni({UR1_tab_D, UR2_tab_D}, {Reg1_tab_D, Reg2_tab_D}, 'Wykres metody uproszczonej tablicowej', {'UR1 optimum modułu', 'UR2 symetrycznego potimum'}, 1);

% === porównanie pomiędzy metodami obliczeniowymi
gen_wykresy_uni({UR1_tab_D, UR1_ukl_D}, {Reg1_tab_D, Reg1_ukl_D}, 'Porównianie medody układu zamkniętego i uproszczonej dla kryterium optimum modułu', {'UR1 ukł zamkniety', 'UR2 uproszczony tablicowa'}, 1);

gen_wykresy_uni({UR2_tab_D, UR2_ukl_D}, {Reg2_tab_D, Reg2_ukl_D}, 'Porównianie medody układu zamkniętego i uproszczonej dla kryterium symetrycznego optimum ', {'UR1 ukł zamkniety', 'UR2 uproszczony tablicowa'}, 1);


% == Wygenerowanie tabeli z parametrami
tab_dyskretne = gen_tabela_parametry({UR1_ukl_D, UR2_ukl_D, UR1_tab_D, UR2_tab_D}, {'UR1_ukl_D', 'UR2_ukl_D', 'UR1_tab_D', 'UR2_tab_D'}, 'Tab 3. Paramtery dla regulatorów dyskretnych');

% --- ANALIZA PORÓWNAWCZA (Ciągłe vs Dyskretne) ---
disp('Porównanie dla metody Optimum Modułu');
% 1. Porównanie dla metody Optimum Modułu (układowej)
gen_tabela_porownawcza(tab_ciagle, tab_dyskretne, 'UR1_ukl', 'UR1_ukl_D', 'tab 4. Wpływ dyskretyzacji: Modułowe Układowe');

% 2. Porównanie dla metody Symetrycznej (układowej)
gen_tabela_porownawcza(tab_ciagle, tab_dyskretne, 'UR2_ukl', 'UR2_ukl_D', 'tab 5. Wpływ dyskretyzacji: Symetryczne Układowe');

% 3. Porównanie metod między sobą (wewnątrz domeny dyskretnej)
% Tutaj sprawdzamy jak wypada Symetryczne względem Modułowego po dyskretyzacji
gen_tabela_porownawcza(tab_dyskretne, tab_dyskretne, 'UR1_ukl_D', 'UR2_ukl_D', 'tab 6. Porównanie metod dyskretnych: Modułowe vs Symetryczne');

disp('Porównanie dla metody uproszczonej tablicowej');
% 4. Wpływ dyskretyzacji konkretnie dla regulatorów tablicowych
gen_tabela_porownawcza(tab_ciagle, tab_dyskretne, 'UR1_tab', 'UR1_tab_D', 'tab 7. Wpływ dyskretyzacji: Optimum Modułu Tablicowe');
gen_tabela_porownawcza(tab_ciagle, tab_dyskretne, 'UR2_tab', 'UR2_tab_D', 'tab 8, Wpływ dyskretyzacji: Symetryczne Tablicowe');

% 5. Optimum Modułu: Czy metoda uproszczona (tab) różni się od dokładnej (ukl)?
gen_tabela_porownawcza(tab_dyskretne, tab_dyskretne, 'UR1_ukl_D', 'UR1_tab_D', 'tab 9. Dokładność metody: Modułowe Ukl vs Tablicowe');

% 6. Kryterium Symetryczne: Czy wzory uproszczone dają te same wyniki co G_muz?
gen_tabela_porownawcza(tab_dyskretne, tab_dyskretne, 'UR2_ukl_D', 'UR2_tab_D', 'tab 10. Dokładność metody: Symetryczne Ukl vs Tablicowe');



% === testowe ====
tf_DAC = (delta1*s+1); 
G_0 = G_0 *tf_DAC;
G_0 = c2d(G_0, 0.005,'zoh');

%% === Wprowadzenie ograniczenia prądowego bez modelu windup ===

% == Dla metody uproszczonej tablicowej ==
[~, ~,  Reg1_tab_D_sym] = symulacja_obiekt_nieliniowy(Reg1_tab_D, G_0, 0.001, 1.0, 'Metoda: tablicowa Optimum Modułu', false);
[~, ~,  Reg2_tab_D_sym] = symulacja_obiekt_nieliniowy(Reg2_tab_D, G_0, 0.001, 1.0, 'Metoda: tablicowa Symetryczny', false);

% == Dla metody układu zamkniętego ==
[~, ~,  Reg1_ukl_D_sym] = symulacja_obiekt_nieliniowy(Reg1_ukl_D, G_0, 0.001, 1.0, 'Metoda: Układu zamknietego Optimum Modułu', false);
[~, ~,  Reg2_ukl_D_sym] = symulacja_obiekt_nieliniowy(Reg2_ukl_D, G_0, 0.001, 1.0, 'Metoda: Układu zamknietego Symetryczny', false);

% == Wygenerowanie tabpic z parametrami dla symulacji regulatorów == 
tab_symulacje = gen_tabela_parametry({Reg1_tab_D_sym, Reg2_tab_D_sym, Reg1_ukl_D_sym, Reg2_ukl_D_sym}, {'Reg1_tab_D_sym', 'Reg2_tab_D_sym', 'Reg1_ukl_D_sym', 'Reg2_ukl_D_sym'}, 'tab. 11 Parametry Regulatorów w przeprowadzonej symulacji');


% == Wygenerowanie tabel porównawczych danych z symulacji == 
% == metoda uproszczona tablicowa 
gen_tabela_porownawcza(tab_dyskretne, tab_symulacje, 'UR1_tab_D', 'Reg1_tab_D_sym', 'tab 12. Bez ograniczenia vs Z ograniczeniem bez windup');
gen_tabela_porownawcza(tab_dyskretne, tab_symulacje, 'UR2_tab_D', 'Reg2_tab_D_sym', 'tab 13. Bez ograniczenia vs Z ograniczeniem bez windup');
% == metoda układu zamknietego
gen_tabela_porownawcza(tab_dyskretne, tab_symulacje, 'UR1_ukl_D', 'Reg1_ukl_D_sym', 'tab 14. Bez ograniczenia vs Z ograniczeniem bez windup');
gen_tabela_porownawcza(tab_dyskretne, tab_symulacje, 'UR2_ukl_D', 'Reg2_ukl_D_sym', 'tab 15. Bez ograniczenia vs Z ograniczeniem bez windup');

%% === Wprowadzenie ograniczenia prądowego z modelem windup ===

% == Dla metody uproszczonej tablicowej ==
[~, ~,  Reg1_tab_D_sym_wind] = symulacja_obiekt_nieliniowy(Reg1_tab_D, G_0, 0.001, 1.0, 'Metoda: tablicowa Optimum Modułu z windup', true);
[~, ~,  Reg2_tab_D_sym_wind] = symulacja_obiekt_nieliniowy(Reg2_tab_D, G_0, 0.001, 1.0, 'Metoda: tablicowa Symetryczny z windup', true);

% == Dla metody układu zamkniętego ==
[~, ~,  Reg1_ukl_D_sym_wind] = symulacja_obiekt_nieliniowy(Reg1_ukl_D, G_0, 0.001, 1.0, 'Metoda: Układu zamknietego Optimum Modułu z windup', true);
[~, ~,  Reg2_ukl_D_sym_wind] = symulacja_obiekt_nieliniowy(Reg2_ukl_D, G_0, 0.001, 1.0, 'Metoda: Układu zamknietego Symetryczny z windup', true);

% == Wygenerowanie tabpic z parametrami dla symulacji regulatorów == 
tab_symulacje_wind = gen_tabela_parametry({Reg1_tab_D_sym_wind, Reg2_tab_D_sym_wind, Reg1_ukl_D_sym_wind, Reg2_ukl_D_sym_wind}, {'Reg1_tab_D_sym_wind', 'Reg2_tab_D_sym_wind', 'Reg1_ukl_D_sym_wind', 'Reg2_ukl_D_sym_wind'}, 'tab. 16 Parametry Regulatorów w przeprowadzonej symulacji z modułem windup');


% == Wygenerowanie tabel porównawczych danych z symulacji == 
% == metoda uproszczona tablicowa 
gen_tabela_porownawcza(tab_dyskretne, tab_symulacje_wind, 'UR1_tab_D', 'Reg1_tab_D_sym_wind', 'tab 17. Bez ograniczenia vs Z ograniczeniem z wind up');
gen_tabela_porownawcza(tab_dyskretne, tab_symulacje_wind, 'UR2_tab_D', 'Reg2_tab_D_sym_wind', 'tab 18. Bez ograniczenia vs Z ograniczeniem z wind up');
% == metoda układu zamknietego
gen_tabela_porownawcza(tab_dyskretne, tab_symulacje_wind, 'UR1_ukl_D', 'Reg1_ukl_D_sym_wind', 'tab 19. Bez ograniczenia vs Z ograniczeniem z wind up');
gen_tabela_porownawcza(tab_dyskretne, tab_symulacje_wind, 'UR2_ukl_D', 'Reg2_ukl_D_sym_wind', 'tab 20. Bez ograniczenia vs Z ograniczeniem z wind up');

% == porównanie regulatora symulowanego z windup i bez

gen_tabela_porownawcza(tab_symulacje, tab_symulacje_wind, 'Reg1_tab_D_sym', 'Reg1_tab_D_sym_wind', 'tab 21. Regulator symulowany Bez Windup vs  Z wind up');
gen_tabela_porownawcza(tab_symulacje, tab_symulacje_wind, 'Reg2_tab_D_sym', 'Reg2_tab_D_sym_wind', 'tab 22. Regulator symulowany Bez Windup vs  Z wind up');
% == metoda układu zamknietego
gen_tabela_porownawcza(tab_symulacje, tab_symulacje_wind, 'Reg1_ukl_D_sym', 'Reg1_ukl_D_sym_wind', 'tab 23. Regulator symulowany Bez Windup vs  Z wind up');
gen_tabela_porownawcza(tab_symulacje, tab_symulacje_wind, 'Reg2_ukl_D_sym', 'Reg2_ukl_D_sym_wind', 'tab 24. Regulator symulowany Bez Windup vs  Z wind up');

function [Gz_out, Gr_out] = gen_kessler(metoda, kryterium, Tp)
    % Pobieranie parametrów bazowych
    % === parametry silnika ===
    CePhi = 2.62;
    Tm = 0.18;
    delta1 = 2.7*10^(-6);
    delta2 = 250*10^(-9);
    delta3 = 0.005;
    delta_sum = delta1 + delta2 + delta3;
    Ks = 1/CePhi;                                                                                               % Wyznaczenie wzmocnienia modelu silnika
    s = tf('s');                                                                                                % Zdefiniowanie zmiennej Laplace'a
    G_model = Ks / ((Tm*s + 1) * (delta_sum*s + 1));                                                            % Wyznaczenie modelu matematycznego silnika

   
    if strcmp(metoda, 'ukl')                                                                                    % Metoda układu zamkniętego
        if strcmp(kryterium, 'modul')                                                                           % Metoda układu zamkniętego optimum modułu
            G_muz = 1 / (2*delta_sum^2*s^2 + 2*delta_sum*s + 1);                                                % Wyznaczenie G_muz dla metody optimum modułu
        else                                                                                                    % Metoda układu zamkniętego 
            G_muz = (4*delta_sum*s + 1) / (8*delta_sum^3*s^3 + 8*delta_sum^2*s^2 + 4*delta_sum*s + 1);          % wyznaczenie G_muz dla metody układu zamkniętego symetrycznej
        end
        Gr_s = minreal(G_muz / (G_model * (1 - G_muz)));                                                        % Wyznaczenie transmitancji regulatora ciągłego metodą układu zamkniętego
        
    elseif strcmp(metoda, 'tab')                                                                                                        
        if strcmp(kryterium, 'modul')                                                                           % Metoda uproszczona tablicowa
            Kr = 1 / (2 * Ks * delta_sum);                                                                      % Wyznaczenie Kr dla metody tablicowej optimum modułu
            Tr = Tm;                                                                                            % Zdefiniowanie dominującej stałej       
            Gr_s = Kr * ((Tr*s + 1) / s);                                                                    % Wyznaczenie regulatora wyznaczonego metodą optimum modułu
        else                                                                                                    % Metoda uproszczona uproszczona tablicowa symetryczna
            Kr = Tm / (8 * Ks * delta_sum^2);                                                                   % Wyznaczenie Kr dla metody symetrycznej 
            Tr = 4 * delta_sum;                                                                                 % Wyznaczenie czasu zastępczego dla metody tablicowej symetrycznej
            Gr_s = (Kr/s) * (1+Tr*s);                                                              % Wyznaczenie transmitancji regulatora dla metody ablicowej symetrycznej
        end
    end

    
    % Do dyskretyzacji bierzemy transmitancję toru głównego bez C/A
    % ponieważ sama metodo dyskretyzacji ZOH wprowadzi transmitancję C/A do
    % układu
    tf_DAC = (delta1*s+1);                                                                                            % czas propagacji pzrec DAC
    G_real_cont = G_model *tf_DAC;                                                                             % Wyrzucenie opóźnienia DAC z toru głównego

    % --- 4. GENEROWANIE WYNIKU ---
    if nargin < 3 || isempty(Tp)                                                                                % sprawdzenie czy został podany czas próbkowania jeśli nie to generuj ukł ciągły
        Gz_out = feedback(Gr_s * G_model, 1);                                                                   % Zwrócenie transmitancji układu regulacji ciągłego
        Gr_out = Gr_s;                                                                                          % Zwrócenie transmitancji regulatora ciągłego
    else
        Gr_z = c2d(Gr_s, Tp, 'ZOH');                                                                            % Przeprowadzenie dyskretyzacji regulatora
        G_obj_d = c2d(G_real_cont, Tp, 'zoh');                                                                  % Przeprowadzenie dyskretyzacji toru gównego
        Gz_out = feedback(Gr_z * G_obj_d, 1);                                                                   % Połączenie regulatora z torem głównym oraz zamknięcie go zprzężeniem zwrotnym oraz zwrócenie dyskretnej wersji ukłądu reulacji
        Gr_out = Gr_z;                                                                                          % Zwrócenie transmitancji terulatora dyskretnego
    end
end



function gen_wykresy_uni(Gz_list, Gr_list, naglowek, wlasna_legenda, t_final)
    % Zabezpieczenie przed starym formatem wywołania (brak Gr_list)
    if ~iscell(Gr_list)
        if nargin > 3, t_final = wlasna_legenda; end
        if nargin > 2, wlasna_legenda = naglowek; end
        naglowek = Gr_list;
        Gr_list = cell(size(Gz_list));
        for k = 1:length(Gr_list), Gr_list{k} = tf(1); end
    end

    if nargin < 5 || isempty(t_final), t_final = 2; end
    
    n = length(Gz_list);
    figure('Color', 'white', 'Name', naglowek, 'Units', 'normalized', 'Position', [0.05 0.1 0.9 0.8]);
    sgtitle(naglowek, 'FontSize', 14, 'FontWeight', 'bold');
    colors = lines(n);
    w_range = logspace(-1, 3, 500); 

    for i = 1:n
        c = colors(i, :);
        Gz = Gz_list{i};
        Gr = Gr_list{i};
        
        % Obliczenia sygnałów
        [y, t] = step(Gz, t_final);
        y = squeeze(y);
        Gu = minreal(Gr * (1 - Gz));
        [u, tu] = step(Gu, t_final);
        [yi, ti] = impulse(Gz, t_final);
        
        % Wybór stylu rysowania
        func = @plot; if isdt(Gz), func = @stairs; end
        
        % --- Rząd 1: Dziedzina czasu ---
        
        % 1. Wyjście y(t)
        subplot(2, 3, 1); hold on;
        func(t, y, 'Color', c, 'LineWidth', 1.5);
        title('Wyjście y(t)'); grid on; ylabel('Amplituda');

        % 2. Uchyb e(t) = 1 - y(t)
        subplot(2, 3, 4); hold on;
        func(t, 1 - y, 'Color', c, 'LineWidth', 1.5);
        title('Uchyb e(t)'); grid on;

        % 3. Sterowanie u(t)
        subplot(2, 3, 3); hold on;
        func(tu, squeeze(u), 'Color', c, 'LineWidth', 1.5);
        title('Sterowanie u(t)'); grid on;

        % --- Rząd 2: Częstotliwość i Impuls ---

        % 4. Magnituda (Bode)
        subplot(2, 3, 5); hold on;
        [mag, ~] = bode(Gz, w_range);
        semilogx(w_range, 20*log10(squeeze(mag)), 'Color', c, 'LineWidth', 1.5);
        title('Bode: Magnituda [dB]'); grid on; xlabel('rad/s');

        % 5. Faza (Bode)
        subplot(2, 3, 6); hold on;
        [~, ph] = bode(Gz, w_range);
        semilogx(w_range, squeeze(ph), 'Color', c, 'LineWidth', 1.5);
        title('Bode: Faza [deg]'); grid on; xlabel('rad/s');

        % 6. Odpowiedź impulsowa h(t)
        subplot(2, 3, 2); hold on;
        func(ti, squeeze(yi), 'Color', c, 'LineWidth', 1.5);
        title('Impuls h(t)'); grid on;
    end
    
    subplot(2, 3, 1);
    legend(wlasna_legenda, 'Location', 'best', 'FontSize', 8);
end

function [T_wyniki] = gen_tabela_parametry(data_array, wlasna_legenda, tytul)
    % data_array    - komórka z modelami TF/ZPK LUB strukturami stepinfo
    % wlasna_legenda - nazwy dla kolumn tabeli
    % tytul          - nagłówek wyświetlany w konsoli

    n = length(data_array);
    
    % Przygotowanie kontenerów na dane
    rise_times = zeros(n, 1);
    settling_times = zeros(n, 1);
    overshoots = zeros(n, 1);
    ess_errors = zeros(n, 1);
    stability = cell(n, 1);

    disp(' ');
    disp('===========================================================');
    disp([upper(tytul)]);
    disp('===========================================================');

    for i = 1:n
        item = data_array{i};
        
        % Sprawdzenie czy wejście to model liniowy czy struktura stepinfo
        if isstruct(item)
            % PRZYPADEK: Struktura stepinfo
            info = item;
            rise_times(i) = info.RiseTime;
            settling_times(i) = info.SettlingTime;
            overshoots(i) = info.Overshoot;
            
            % Dane niedostępne w samej strukturze stepinfo (domyślne wartości)
            ess_errors(i) = NaN; 
            stability{i} = 'N/A';
            domena = 'Dane ze struktury';
        else
            % PRZYPADEK: Model (np. tf, zpk, ss)
            sys = item;
            info = stepinfo(sys);
            rise_times(i) = info.RiseTime;
            settling_times(i) = info.SettlingTime;
            overshoots(i) = info.Overshoot;
            
            % Obliczanie uchybu i stabilności
            k_stat = dcgain(sys);
            ess_errors(i) = abs(1 - k_stat);
            
            p = pole(sys);
            if isdt(sys)
                is_stable = all(abs(p) < 1.0);
                domena = ['Z (Ts=', num2str(sys.Ts), 's)'];
            else
                is_stable = all(real(p) < 0);
                domena = 'S (Ciągła)';
            end
            if is_stable; stability{i} = 'Tak'; else; stability{i} = 'NIE!'; end
        end
        
        fprintf('Pozycja %d [%s]: %s\n', i, wlasna_legenda{i}, domena);
    end

    % Budowanie tabeli wynikowej
    Metryka = {
        'Czas narastania (t_r) [s]'; 
        'Czas ustalania (t_us) [s]'; 
        'Przeregulowanie [%]'; 
        'Uchyb ustalony';
        'Stabilny'
    };

    dane_num = [rise_times'; settling_times'; overshoots'; ess_errors'];
    T_wyniki = table(Metryka, 'RowNames', {});
    
    for i = 1:n
        % Jeśli uchyb to NaN (dla struktur), wyświetlamy '-'
        ess_str = dane_num(4,i);
        if isnan(ess_str), ess_val = {'-'}; else, ess_val = {ess_str}; end
        
        col_data = {dane_num(1,i); dane_num(2,i); dane_num(3,i); ess_val{1}; stability{i}};
        T_wyniki.(wlasna_legenda{i}) = col_data;
    end

    disp(' ');
    disp(T_wyniki);
    disp('===========================================================');
end

function [T_delta] = gen_tabela_porownawcza(tab1, tab2, col_name1, col_name2, naglowek)
    % tab1, tab2 - tabele wynikowe z funkcji gen_tabela_parametry
    % col_name1  - nazwa kolumny z tab1 (odniesienie/baza)
    % col_name2  - nazwa kolumny z tab2 (badany model)
    % naglowek   - tytuł wyświetlany w konsoli

    % 1. Pobranie metryk (pomijamy ostatni wiersz 'Stabilny')
    Metryka = tab1.Metryka(1:end-1); 
    n_metryk = length(Metryka);
    eps_val = 1e-10;

    V1 = zeros(n_metryk, 1);
    V2 = zeros(n_metryk, 1);

    % 2. Bezpieczna ekstrakcja wartości numerycznych
    for i = 1:n_metryk
        val1 = tab1{i, col_name1};
        val2 = tab2{i, col_name2};
        
        % Jeśli wartość jest w komórce, wyciągnij ją
        if iscell(val1), val1 = val1{1}; end
        if iscell(val2), val2 = val2{1}; end
        
        % Konwersja na double (jeśli to '-' lub inny tekst, zamień na NaN)
        if ischar(val1) || isstring(val1), V1(i) = NaN; else, V1(i) = val1; end
        if ischar(val2) || isstring(val2), V2(i) = NaN; else, V2(i) = val2; end
    end

    % 3. Obliczenie różnicy procentowej
    % Używamy abs(V1), aby poprawnie liczyć zmiany dla wartości ujemnych
    Roznica_proc = ((V2 - V1) ./ (abs(V1) + eps_val)) * 100;

    % 4. Tworzenie tabeli wynikowej
    T_delta = table(Metryka, V1, V2, Roznica_proc, ...
        'VariableNames', {'Parametr', 'Baza_Referencyjna', 'Model_Badany', 'Zmiana_procentowa'});

    % 5. Wyświetlanie raportu w konsoli
    disp(' ');
    disp('================================================================');
    disp(['ANALIZA PORÓWNAWCZA: ', upper(naglowek)]);
    disp(['Baza: ', col_name1, ' | Test: ', col_name2]);
    disp('----------------------------------------------------------------');
    
    % Formatowanie wyświetlania: NaN zostaną pokazane jako blanki lub NaN
    disp(T_delta);
    
    disp('================================================================');
end

% function simulate_real_motor(Gr_z, Tp, t_final, tytul)
%     % Parametry fizyczne silnika 6V
%     R = 1.5;          % Rezystancja [Ohm]
%     ke = 0.02;        % Stała napięciowa [V*s/rad]
%     Tm = 0.18;        % Stała czasowa mechaniczna [s]
%     Ks = 1/2.62;      % Wzmocnienie statyczne
% 
%     % Limity fizyczne (L298 i zasilanie)
%     U_max = 6.0;      % Napięcie zasilania [V]
%     I_limit = 4.0;    % Ograniczenie prądu [A]
% 
%     t = 0:Tp:t_final;
%     n = length(t);
%     y = zeros(1, n);       % Prędkość
%     u_reg = zeros(1, n);   % Wyjście regulatora (idealne)
%     u_sat = zeros(1, n);   % Wyjście po nasyceniu (realne)
%     i_motor = zeros(1, n); % Prąd
%     e = zeros(1, n);       % Uchyb
% 
%     % Pobranie współczynników równania różnicowego regulatora
%     [num, den] = tfdata(Gr_z, 'v');
% 
%     for k = 2:n
%         % 1. Oblicz uchyb
%         e(k) = 1.0 - y(k-1);
% 
%         % 2. Równanie różnicowe regulatora (u[k] = f(e[k], e[k-1], u[k-1]))
%         % Uwzględniamy u_sat[k-1] zamiast u_reg[k-1] (uproszczony Anti-Windup)
%         u_reg(k) = (num(1)*e(k) + num(2)*e(k-1) - den(2)*u_sat(k-1)) / den(1);
% 
%         % 3. Nasycenie napięcia (Mostek H)
%         u_temp = u_reg(k);
%         if u_temp > U_max, u_temp = U_max; end
%         if u_temp < -U_max, u_temp = -U_max; end
% 
%         % 4. Ograniczenie prądu (Model uproszczony)
%         % i = (U - E) / R -> U = I*R + E
%         current_now = (u_temp - ke * y(k-1)) / R;
%         if abs(current_now) > I_limit
%             current_now = sign(current_now) * I_limit;
%             u_temp = current_now * R + ke * y(k-1); % Wyliczenie napięcia granicznego
%         end
% 
%         u_sat(k) = u_temp;
%         i_motor(k) = current_now;
% 
%         % 5. Model silnika (Euler: y[k] = y[k-1] + dy*Tp)
%         dy = (Ks * u_sat(k) - y(k-1)) / Tm;
%         y(k) = y(k-1) + dy * Tp;
%     end
% 
%     % --- WYKRESY ---
%     figure('Color', 'white', 'Name', tytul, 'Units', 'normalized', 'Position', [0.1 0.1 0.5 0.8]);
%     sgtitle(['Analiza nieliniowa: ', tytul], 'FontSize', 12, 'FontWeight', 'bold');
% 
%     % Wyjście (Prędkość)
%     subplot(3, 1, 1);
%     plot(t, y, 'b', 'LineWidth', 2); hold on;
%     line([0 t_final], [1 1], 'Color', 'r', 'LineStyle', '--'); % Wartość zadana
%     title('Prędkość obrotowa y(t)'); ylabel('\omega [rad/s]'); grid on;
% 
%     % Sterowanie (Napięcie)
%     subplot(3, 1, 2);
%     plot(t, u_reg, 'g--', 'LineWidth', 1); hold on;
%     plot(t, u_sat, 'r', 'LineWidth', 1.5);
%     title('Napięcie sterujące u(t)'); ylabel('U [V]'); grid on;
%     legend('Zapotrzebowanie reg.', 'Napięcie nasycone', 'Location', 'best');
% 
%     % Prąd
%     subplot(3, 1, 3);
%     plot(t, i_motor, 'k', 'LineWidth', 1.5);
%     title('Natężenie prądu i(t)'); ylabel('I [A]'); xlabel('Czas [s]'); grid on;
%     ylim([-I_limit*1.2, I_limit*1.2]);
% end

% function [t, y, info] = gen_silnik_nieliniowy(Gr_z, Tp, t_final, tytul)
%     % --- Parametry silnika ---
%     R = 1.5; ke = 0.02; Tm = 0.18; Ks = 1/2.62;
%     U_limit = 6.0; I_limit = 4.0;
% 
%     % --- Przygotowanie regulatora (State-Space dla stabilności) ---
%     [A, B, C, D] = ssdata(Gr_z);
%     x_reg = zeros(size(A,1), 1);
% 
%     t = 0:Tp:t_final;
%     n = length(t);
%     y = zeros(1, n); u_s = zeros(1, n); i_m = zeros(1, n);
% 
%     % --- Pętla symulacji ---
%     for k = 2:n
%         err = 1.0 - y(k-1);
%         u_ideal = C * x_reg + D * err;
% 
%         % Nieliniowość prądowo-napięciowa
%         i_pot = (u_ideal - ke * y(k-1)) / R;
%         if abs(i_pot) > I_limit
%             u_f = sign(i_pot) * I_limit * R + ke * y(k-1);
%             if abs(u_f) > U_limit, u_f = sign(u_f) * U_limit; end
%         else
%             u_f = max(min(u_ideal, U_limit), -U_limit);
%         end
% 
%         u_s(k) = u_f;
%         i_m(k) = (u_s(k) - ke * y(k-1)) / R;
% 
%         % Aktualizacja stanu regulatora i obiektu
%         x_reg = A * x_reg + B * err;
%         dy = (Ks * u_s(k) - y(k-1)) / Tm;
%         y(k) = y(k-1) + dy * Tp;
%     end
% 
%     % --- ANALIZA JAKOŚCI I TABELA ---
%     % stepinfo oblicza czasy: narastania, ustalania, przeregulowanie
%     info = stepinfo(y, t, 1.0);
% 
%     % Wyznaczenie maksymalnego prądu i napięcia z symulacji
%     max_I = max(abs(i_m));
%     max_U = max(abs(u_s));
% 
%     % Tworzenie tabeli
%     Parametr = {'Czas narastania [s]'; 'Czas ustalania (2%) [s]'; 'Przeregulowanie [%]'; 'Max Prąd [A]'; 'Max Napięcie [V]'};
%     Wartosc = [info.RiseTime; info.SettlingTime; info.Overshoot; max_I; max_U];
% 
%     TabelaWynikow = table(Wartosc, 'RowNames', Parametr);
% 
%     % Wyświetlenie nagłówka i tabeli w konsoli
%     fprintf('\n--- Wyniki symulacji dla: %s ---\n', tytul);
%     disp(TabelaWynikow);
% 
%     % --- Wykresy (bez zmian) ---
%     figure('Name', tytul, 'Color', 'w');
%     subplot(3,1,1); plot(t, y, 'b', 'LineWidth', 1.5); hold on;
%     line([0 t_final], [1 1], 'Color', 'r', 'LineStyle', '--');
%     title(['Wyjście (Prędkość): ', tytul]); grid on;
% 
%     subplot(3,1,2); plot(t, u_s, 'r', 'LineWidth', 1.2); 
%     title('Napięcie sterujące u(t) [V]'); grid on;
% 
%     subplot(3,1,3); plot(t, i_m, 'k', 'LineWidth', 1.2); 
%     title('Prąd silnika i(t) [A]'); grid on; xlabel('Czas [s]');
% end
% 


%próba teraz nie wiem która
function [t, y, info] = symulacja_obiekt_nieliniowy(Gr_z, G_obj_d, Tp, t_final, tytul, anti_windup)

    % --- Parametry fizyczne ---
    R = 1.5;
    ke = 0.02;
    U_limit = 6.0;
    I_limit = 4.0;

    % --- State-space regulatora ---
    [Ar, Br, Cr, Dr] = ssdata(ss(Gr_z));
    xr = zeros(size(Ar,1),1);

    % --- State-space obiektu ---
    [Ao, Bo, Co, Do] = ssdata(ss(G_obj_d));
    xo = zeros(size(Ao,1),1);

    % --- Czas ---
    t = 0:Tp:t_final;
    n = length(t);

    % --- Bufory ---
    y = zeros(1,n);
    u = zeros(1,n);
    i_m = zeros(1,n);

    % --- Symulacja ---
    for k = 2:n
        
        % Uchyb
        e = 1 - y(k-1);

        % --- Regulator ---
        u_ideal = Cr*xr + Dr*e;

        % ==========================================
        % NIELINIOWOŚĆ (prąd + napięcie)
        % ==========================================

        saturation = false;

        i_temp = (u_ideal - ke*y(k-1)) / R;

        % --- Ograniczenie prądu ---
        if abs(i_temp) > I_limit

            saturation = true;

            i_temp = sign(i_temp) * I_limit;
            u_sat = i_temp * R + ke*y(k-1);

        else
            u_sat = u_ideal;
        end

        % --- Ograniczenie napięcia ---
        if abs(u_sat) > U_limit

            saturation = true;

            u_sat = max(min(u_sat, U_limit), -U_limit);
        end

        % --- Sygnały ---
        u(k) = u_sat;
        i_m(k) = (u(k) - ke*y(k-1)) / R;

        % ==========================================
        % ANTI-WINDUP
        % ==========================================

        if anti_windup

            % aktualizacja regulatora tylko
            % gdy NIE ma saturacji
            if ~saturation
                xr = Ar*xr + Br*e;
            end

        else

            % klasyczny regulator bez anti-windup
            xr = Ar*xr + Br*e;

        end

        % --- Obiekt ---
        xo = Ao*xo + Bo*u(k);

        % --- Wyjście ---
        y(k) = Co*xo + Do*u(k);

    end

    % --- Analiza ---
    info = stepinfo(y, t, 1);

    % --- Wykresy ---
    figure('Name', tytul, 'Color', 'w');

    subplot(3,1,1)
    plot(t, y, 'b', 'LineWidth', 1.5); hold on;
    yline(1,'r--');
    title(['Wyjście y(t): ', tytul]);
    grid on;

    subplot(3,1,2)
    plot(t, u, 'r', 'LineWidth', 1.5);
    title('Sterowanie u(t)');
    grid on;

    subplot(3,1,3)
    plot(t, i_m, 'k', 'LineWidth', 1.5);
    title('Prąd i(t)');
    grid on;
    xlabel('Czas [s]');

end

%% ===  Wnioski Ilościowe === 
% ===========================================================
% TAB 1. PARAMETRY SILNIKA
% ===========================================================
% Pozycja 1 [G_0]: S (Ciągła)
% 
%                Metryka                  G_0    
%     _____________________________    __________
% 
%     {'Czas narastania (t_r) [s]'}    {[0.3956]}
%     {'Czas ustalania (t_us) [s]'}    {[0.7093]}
%     {'Przeregulowanie [%]'      }    {[     0]}
%     {'Uchyb ustalony'           }    {[0.6183]}
%     {'Stabilny'                 }    {'Tak'   }
% 
% ===========================================================
% 
% ===========================================================
% TAB 2. PARAMTERY DLA REGULATORÓW CIĄGŁYCH
% ===========================================================
% Pozycja 1 [UR1_ukl]: S (Ciągła)
% Pozycja 2 [UR2_ukl]: S (Ciągła)
% Pozycja 3 [UR1_tab]: S (Ciągła)
% Pozycja 4 [UR2_tab]: S (Ciągła)
% 
%                Metryka                UR1_ukl        UR2_ukl       UR1_tab        UR2_tab  
%     _____________________________    __________    ___________    __________    ___________
% 
%     {'Czas narastania (t_r) [s]'}    {[0.0152]}    {[ 0.0106]}    {[0.0152]}    {[ 0.0109]}
%     {'Czas ustalania (t_us) [s]'}    {[0.0422]}    {[ 0.0828]}    {[0.0422]}    {[ 0.0798]}
%     {'Przeregulowanie [%]'      }    {[4.3210]}    {[43.3919]}    {[4.3210]}    {[37.5717]}
%     {'Uchyb ustalony'           }    {[     0]}    {[      0]}    {[     0]}    {[      0]}
%     {'Stabilny'                 }    {'Tak'   }    {'Tak'    }    {'Tak'   }    {'Tak'    }
% 
% ===========================================================
% 
% ===========================================================
% TAB 3. PARAMTERY DLA REGULATORÓW DYSKRETNYCH
% ===========================================================
% Pozycja 1 [UR1_ukl_D]: Z (Ts=0.001s)
% Pozycja 2 [UR2_ukl_D]: Z (Ts=0.001s)
% Pozycja 3 [UR1_tab_D]: Z (Ts=0.001s)
% Pozycja 4 [UR2_tab_D]: Z (Ts=0.001s)
% 
%                Metryka                 UR1_ukl_D         UR2_ukl_D         UR1_tab_D         UR2_tab_D   
%     _____________________________    ______________    ______________    ______________    ______________
% 
%     {'Czas narastania (t_r) [s]'}    {[    0.0150]}    {[    0.0100]}    {[    0.0150]}    {[    0.0110]}
%     {'Czas ustalania (t_us) [s]'}    {[    0.0430]}    {[    0.0840]}    {[    0.0430]}    {[    0.0820]}
%     {'Przeregulowanie [%]'      }    {[    5.7977]}    {[   48.2611]}    {[    5.7977]}    {[   41.8535]}
%     {'Uchyb ustalony'           }    {[3.1131e-13]}    {[1.3462e-10]}    {[3.1131e-13]}    {[2.2105e-13]}
%     {'Stabilny'                 }    {'Tak'       }    {'Tak'       }    {'Tak'       }    {'Tak'       }
% 
% ===========================================================
% Porównanie dla metody Optimum Modułu
% 
% ================================================================
% ANALIZA PORÓWNAWCZA: TAB 4. WPŁYW DYSKRETYZACJI: MODUŁOWE UKŁADOWE
% Baza: UR1_ukl | Test: UR1_ukl_D
% ----------------------------------------------------------------
%               Parametr               Baza_Referencyjna    Model_Badany    Zmiana_procentowa
%     _____________________________    _________________    ____________    _________________
% 
%     {'Czas narastania (t_r) [s]'}        0.015204               0.015          -1.3412     
%     {'Czas ustalania (t_us) [s]'}        0.042188               0.043           1.9249     
%     {'Przeregulowanie [%]'      }           4.321              5.7977           34.176     
%     {'Uchyb ustalony'           }               0          3.1131e-13          0.31131     
% 
% ================================================================
% 
% ================================================================
% ANALIZA PORÓWNAWCZA: TAB 5. WPŁYW DYSKRETYZACJI: SYMETRYCZNE UKŁADOWE
% Baza: UR2_ukl | Test: UR2_ukl_D
% ----------------------------------------------------------------
%               Parametr               Baza_Referencyjna    Model_Badany    Zmiana_procentowa
%     _____________________________    _________________    ____________    _________________
% 
%     {'Czas narastania (t_r) [s]'}        0.010596                0.01          -5.6213     
%     {'Czas ustalania (t_us) [s]'}        0.082803               0.084           1.4452     
%     {'Przeregulowanie [%]'      }          43.392              48.261           11.221     
%     {'Uchyb ustalony'           }               0          1.3462e-10           134.62     
% 
% ================================================================
% 
% ================================================================
% ANALIZA PORÓWNAWCZA: TAB 6. PORÓWNANIE METOD DYSKRETNYCH: MODUŁOWE VS SYMETRYCZNE
% Baza: UR1_ukl_D | Test: UR2_ukl_D
% ----------------------------------------------------------------
%               Parametr               Baza_Referencyjna    Model_Badany    Zmiana_procentowa
%     _____________________________    _________________    ____________    _________________
% 
%     {'Czas narastania (t_r) [s]'}            0.015               0.01          -33.333     
%     {'Czas ustalania (t_us) [s]'}            0.043              0.084           95.349     
%     {'Przeregulowanie [%]'      }           5.7977             48.261           732.42     
%     {'Uchyb ustalony'           }       3.1131e-13         1.3462e-10           133.89     
% 
% ================================================================
% Porównanie dla metody uproszczonej tablicowej
% 
% ================================================================
% ANALIZA PORÓWNAWCZA: TAB 7. WPŁYW DYSKRETYZACJI: OPTIMUM MODUŁU TABLICOWE
% Baza: UR1_tab | Test: UR1_tab_D
% ----------------------------------------------------------------
%               Parametr               Baza_Referencyjna    Model_Badany    Zmiana_procentowa
%     _____________________________    _________________    ____________    _________________
% 
%     {'Czas narastania (t_r) [s]'}        0.015204               0.015          -1.3412     
%     {'Czas ustalania (t_us) [s]'}        0.042188               0.043           1.9249     
%     {'Przeregulowanie [%]'      }           4.321              5.7977           34.176     
%     {'Uchyb ustalony'           }               0          3.1131e-13          0.31131     
% 
% ================================================================
% 
% ================================================================
% ANALIZA PORÓWNAWCZA: TAB 8, WPŁYW DYSKRETYZACJI: SYMETRYCZNE TABLICOWE
% Baza: UR2_tab | Test: UR2_tab_D
% ----------------------------------------------------------------
%               Parametr               Baza_Referencyjna    Model_Badany    Zmiana_procentowa
%     _____________________________    _________________    ____________    _________________
% 
%     {'Czas narastania (t_r) [s]'}        0.010922               0.011          0.71218     
%     {'Czas ustalania (t_us) [s]'}        0.079826               0.082           2.7232     
%     {'Przeregulowanie [%]'      }          37.572              41.853           11.396     
%     {'Uchyb ustalony'           }               0          2.2105e-13          0.22105     
% 
% ================================================================
% 
% ================================================================
% ANALIZA PORÓWNAWCZA: TAB 9. DOKŁADNOŚĆ METODY: MODUŁOWE UKL VS TABLICOWE
% Baza: UR1_ukl_D | Test: UR1_tab_D
% ----------------------------------------------------------------
%               Parametr               Baza_Referencyjna    Model_Badany    Zmiana_procentowa
%     _____________________________    _________________    ____________    _________________
% 
%     {'Czas narastania (t_r) [s]'}            0.015              0.015             0        
%     {'Czas ustalania (t_us) [s]'}            0.043              0.043             0        
%     {'Przeregulowanie [%]'      }           5.7977             5.7977             0        
%     {'Uchyb ustalony'           }       3.1131e-13         3.1131e-13             0        
% 
% ================================================================
% 
% ================================================================
% ANALIZA PORÓWNAWCZA: TAB 10. DOKŁADNOŚĆ METODY: SYMETRYCZNE UKL VS TABLICOWE
% Baza: UR2_ukl_D | Test: UR2_tab_D
% ----------------------------------------------------------------
%               Parametr               Baza_Referencyjna    Model_Badany    Zmiana_procentowa
%     _____________________________    _________________    ____________    _________________
% 
%     {'Czas narastania (t_r) [s]'}             0.01              0.011               10     
%     {'Czas ustalania (t_us) [s]'}            0.084              0.082           -2.381     
%     {'Przeregulowanie [%]'      }           48.261             41.853          -13.277     
%     {'Uchyb ustalony'           }       1.3462e-10         2.2105e-13          -57.284     
% 
% ================================================================
% 
% ===========================================================
% TAB. 11 PARAMETRY REGULATORÓW W PRZEPROWADZONEJ SYMULACJI
% ===========================================================
% Pozycja 1 [Reg1_tab_D_sym]: Dane ze struktury
% Pozycja 2 [Reg2_tab_D_sym]: Dane ze struktury
% Pozycja 3 [Reg1_ukl_D_sym]: Dane ze struktury
% Pozycja 4 [Reg2_ukl_D_sym]: Dane ze struktury
% 
%                Metryka               Reg1_tab_D_sym    Reg2_tab_D_sym    Reg1_ukl_D_sym    Reg2_ukl_D_sym
%     _____________________________    ______________    ______________    ______________    ______________
% 
%     {'Czas narastania (t_r) [s]'}      {[0.0165]}       {[ 0.0165]}        {[0.0165]}       {[ 0.0165]}  
%     {'Czas ustalania (t_us) [s]'}      {[0.0212]}       {[ 0.0891]}        {[0.0212]}       {[ 0.0899]}  
%     {'Przeregulowanie [%]'      }      {[1.8923]}       {[37.5846]}        {[1.8923]}       {[44.4128]}  
%     {'Uchyb ustalony'           }      {'-'     }       {'-'      }        {'-'     }       {'-'      }  
%     {'Stabilny'                 }      {'N/A'   }       {'N/A'    }        {'N/A'   }       {'N/A'    }  
% 
% ===========================================================
% 
% ================================================================
% ANALIZA PORÓWNAWCZA: TAB 12. BEZ OGRANICZENIA VS Z OGRANICZENIEM BEZ WINDUP
% Baza: UR1_tab_D | Test: Reg1_tab_D_sym
% ----------------------------------------------------------------
%               Parametr               Baza_Referencyjna    Model_Badany    Zmiana_procentowa
%     _____________________________    _________________    ____________    _________________
% 
%     {'Czas narastania (t_r) [s]'}            0.015          0.016458            9.7233     
%     {'Czas ustalania (t_us) [s]'}            0.043          0.021202           -50.692     
%     {'Przeregulowanie [%]'      }           5.7977            1.8923           -67.361     
%     {'Uchyb ustalony'           }       3.1131e-13               NaN               NaN     
% 
% ================================================================
% 
% ================================================================
% ANALIZA PORÓWNAWCZA: TAB 13. BEZ OGRANICZENIA VS Z OGRANICZENIEM BEZ WINDUP
% Baza: UR2_tab_D | Test: Reg2_tab_D_sym
% ----------------------------------------------------------------
%               Parametr               Baza_Referencyjna    Model_Badany    Zmiana_procentowa
%     _____________________________    _________________    ____________    _________________
% 
%     {'Czas narastania (t_r) [s]'}            0.011          0.016458           49.623      
%     {'Czas ustalania (t_us) [s]'}            0.082          0.089083           8.6375      
%     {'Przeregulowanie [%]'      }           41.853            37.585            -10.2      
%     {'Uchyb ustalony'           }       2.2105e-13               NaN              NaN      
% 
% ================================================================
% 
% ================================================================
% ANALIZA PORÓWNAWCZA: TAB 14. BEZ OGRANICZENIA VS Z OGRANICZENIEM BEZ WINDUP
% Baza: UR1_ukl_D | Test: Reg1_ukl_D_sym
% ----------------------------------------------------------------
%               Parametr               Baza_Referencyjna    Model_Badany    Zmiana_procentowa
%     _____________________________    _________________    ____________    _________________
% 
%     {'Czas narastania (t_r) [s]'}            0.015          0.016458            9.7233     
%     {'Czas ustalania (t_us) [s]'}            0.043          0.021202           -50.692     
%     {'Przeregulowanie [%]'      }           5.7977            1.8923           -67.361     
%     {'Uchyb ustalony'           }       3.1131e-13               NaN               NaN     
% 
% ================================================================
% 
% ================================================================
% ANALIZA PORÓWNAWCZA: TAB 15. BEZ OGRANICZENIA VS Z OGRANICZENIEM BEZ WINDUP
% Baza: UR2_ukl_D | Test: Reg2_ukl_D_sym
% ----------------------------------------------------------------
%               Parametr               Baza_Referencyjna    Model_Badany    Zmiana_procentowa
%     _____________________________    _________________    ____________    _________________
% 
%     {'Czas narastania (t_r) [s]'}             0.01          0.016458            64.585     
%     {'Czas ustalania (t_us) [s]'}            0.084          0.089914            7.0407     
%     {'Przeregulowanie [%]'      }           48.261            44.413           -7.9738     
%     {'Uchyb ustalony'           }       1.3462e-10               NaN               NaN     
% 
% ================================================================
% 
% ===========================================================
% TAB. 16 PARAMETRY REGULATORÓW W PRZEPROWADZONEJ SYMULACJI Z MODUŁEM WINDUP
% ===========================================================
% Pozycja 1 [Reg1_tab_D_sym_wind]: Dane ze struktury
% Pozycja 2 [Reg2_tab_D_sym_wind]: Dane ze struktury
% Pozycja 3 [Reg1_ukl_D_sym_wind]: Dane ze struktury
% Pozycja 4 [Reg2_ukl_D_sym_wind]: Dane ze struktury
% 
%                Metryka               Reg1_tab_D_sym_wind    Reg2_tab_D_sym_wind    Reg1_ukl_D_sym_wind    Reg2_ukl_D_sym_wind
%     _____________________________    ___________________    ___________________    ___________________    ___________________
% 
%     {'Czas narastania (t_r) [s]'}        {[0.0165]}             {[0.0165]}             {[0.0165]}             {[0.0165]}     
%     {'Czas ustalania (t_us) [s]'}        {[0.2022]}             {[0.0388]}             {[0.2022]}             {[0.0356]}     
%     {'Przeregulowanie [%]'      }        {[     0]}             {[     0]}             {[     0]}             {[0.3710]}     
%     {'Uchyb ustalony'           }        {'-'     }             {'-'     }             {'-'     }             {'-'     }     
%     {'Stabilny'                 }        {'N/A'   }             {'N/A'   }             {'N/A'   }             {'N/A'   }     
% 
% ===========================================================
% 
% ================================================================
% ANALIZA PORÓWNAWCZA: TAB 17. BEZ OGRANICZENIA VS Z OGRANICZENIEM Z WINDUP
% Baza: UR1_tab_D | Test: Reg1_tab_D_sym_wind
% ----------------------------------------------------------------
%               Parametr               Baza_Referencyjna    Model_Badany    Zmiana_procentowa
%     _____________________________    _________________    ____________    _________________
% 
%     {'Czas narastania (t_r) [s]'}            0.015          0.016458           9.7233      
%     {'Czas ustalania (t_us) [s]'}            0.043           0.20217           370.16      
%     {'Przeregulowanie [%]'      }           5.7977                 0             -100      
%     {'Uchyb ustalony'           }       3.1131e-13               NaN              NaN      
% 
% ================================================================
% 
% ================================================================
% ANALIZA PORÓWNAWCZA: TAB 18. BEZ OGRANICZENIA VS Z OGRANICZENIEM Z WINDUP
% Baza: UR2_tab_D | Test: Reg2_tab_D_sym_wind
% ----------------------------------------------------------------
%               Parametr               Baza_Referencyjna    Model_Badany    Zmiana_procentowa
%     _____________________________    _________________    ____________    _________________
% 
%     {'Czas narastania (t_r) [s]'}            0.011          0.016458           49.623      
%     {'Czas ustalania (t_us) [s]'}            0.082          0.038753           -52.74      
%     {'Przeregulowanie [%]'      }           41.853                 0             -100      
%     {'Uchyb ustalony'           }       2.2105e-13               NaN              NaN      
% 
% ================================================================
% 
% ================================================================
% ANALIZA PORÓWNAWCZA: TAB 19. BEZ OGRANICZENIA VS Z OGRANICZENIEM Z WINDUP
% Baza: UR1_ukl_D | Test: Reg1_ukl_D_sym_wind
% ----------------------------------------------------------------
%               Parametr               Baza_Referencyjna    Model_Badany    Zmiana_procentowa
%     _____________________________    _________________    ____________    _________________
% 
%     {'Czas narastania (t_r) [s]'}            0.015          0.016458           9.7233      
%     {'Czas ustalania (t_us) [s]'}            0.043           0.20217           370.16      
%     {'Przeregulowanie [%]'      }           5.7977                 0             -100      
%     {'Uchyb ustalony'           }       3.1131e-13               NaN              NaN      
% 
% ================================================================
% 
% ================================================================
% ANALIZA PORÓWNAWCZA: TAB 20. BEZ OGRANICZENIA VS Z OGRANICZENIEM Z WINDUP
% Baza: UR2_ukl_D | Test: Reg2_ukl_D_sym_wind
% ----------------------------------------------------------------
%               Parametr               Baza_Referencyjna    Model_Badany    Zmiana_procentowa
%     _____________________________    _________________    ____________    _________________
% 
%     {'Czas narastania (t_r) [s]'}             0.01          0.016458            64.585     
%     {'Czas ustalania (t_us) [s]'}            0.084           0.03563           -57.584     
%     {'Przeregulowanie [%]'      }           48.261           0.37095           -99.231     
%     {'Uchyb ustalony'           }       1.3462e-10               NaN               NaN     

% ================================================================
% ANALIZA PORÓWNAWCZA: TAB 21. REGULATOR SYMULOWANY BEZ WINDUP VS  Z WIND UP
% Baza: Reg1_tab_D_sym | Test: Reg1_tab_D_sym_wind
% ----------------------------------------------------------------
%               Parametr               Baza_Referencyjna    Model_Badany    Zmiana_procentowa
%     _____________________________    _________________    ____________    _________________
% 
%     {'Czas narastania (t_r) [s]'}        0.016458           0.016458                0      
%     {'Czas ustalania (t_us) [s]'}        0.021202            0.20217           853.51      
%     {'Przeregulowanie [%]'      }          1.8923                  0             -100      
%     {'Uchyb ustalony'           }             NaN                NaN              NaN      
% 
% ================================================================
% 
% ================================================================
% ANALIZA PORÓWNAWCZA: TAB 22. REGULATOR SYMULOWANY BEZ WINDUP VS  Z WIND UP
% Baza: Reg2_tab_D_sym | Test: Reg2_tab_D_sym_wind
% ----------------------------------------------------------------
%               Parametr               Baza_Referencyjna    Model_Badany    Zmiana_procentowa
%     _____________________________    _________________    ____________    _________________
% 
%     {'Czas narastania (t_r) [s]'}        0.016458           0.016458                 0     
%     {'Czas ustalania (t_us) [s]'}        0.089083           0.038753           -56.498     
%     {'Przeregulowanie [%]'      }          37.585                  0              -100     
%     {'Uchyb ustalony'           }             NaN                NaN               NaN     
% 
% ================================================================
% 
% ================================================================
% ANALIZA PORÓWNAWCZA: TAB 23. REGULATOR SYMULOWANY BEZ WINDUP VS  Z WIND UP
% Baza: Reg1_ukl_D_sym | Test: Reg1_ukl_D_sym_wind
% ----------------------------------------------------------------
%               Parametr               Baza_Referencyjna    Model_Badany    Zmiana_procentowa
%     _____________________________    _________________    ____________    _________________
% 
%     {'Czas narastania (t_r) [s]'}        0.016458           0.016458                0      
%     {'Czas ustalania (t_us) [s]'}        0.021202            0.20217           853.51      
%     {'Przeregulowanie [%]'      }          1.8923                  0             -100      
%     {'Uchyb ustalony'           }             NaN                NaN              NaN      
% 
% ================================================================
% 
% ================================================================
% ANALIZA PORÓWNAWCZA: TAB 24. REGULATOR SYMULOWANY BEZ WINDUP VS  Z WIND UP
% Baza: Reg2_ukl_D_sym | Test: Reg2_ukl_D_sym_wind
% ----------------------------------------------------------------
%               Parametr               Baza_Referencyjna    Model_Badany    Zmiana_procentowa
%     _____________________________    _________________    ____________    _________________
% 
%     {'Czas narastania (t_r) [s]'}        0.016458           0.016458                 0     
%     {'Czas ustalania (t_us) [s]'}        0.089914            0.03563           -60.374     
%     {'Przeregulowanie [%]'      }          44.413            0.37095           -99.165     
%     {'Uchyb ustalony'           }             NaN                NaN               NaN     
% 
% ================================================================


%% === Wnioski jakościowe ===

% a) Poprawa dynamiki układu po zastosowaniu regulatorów Kesslera -
% Zastosowanie regulatorów Kesslera znacząco poprawiło parametry odpowiedzi
% względem modelu wyjściowego G_0. Dla obiektu bez regulatora czas narastania
% wynosił 0.3956s, czas ustalania 0.7093s, a uchyb ustalony 0.6183
% (tab. 1). Po zastosowaniu regulatorów ciągłych czas narastania zmniejszył
% się do zakresu 0.0106-0.0152s, czas ustalania do zakresu 0.0422-0.0828s,
% a uchyb ustalony został sprowadzony do 0 (tab. 2). Świadczy to o tym,
% że regulatory skutecznie przyspieszyły odpowiedź układu oraz usunęły
% błąd statyczny.

% b) Wzrost przeregulowania dla optimum symetrycznego w układzie ciągłym -
% Regulatory dobrane według optimum symetrycznego zapewniły krótszy czas
% narastania, ale jednocześnie zwiększyły przeregulowanie. Dla UR2_ukl
% czas narastania wyniósł 0.0106s, podczas gdy dla UR1_ukl wyniósł
% 0.0152s (tab. 2). Jednocześnie przeregulowanie wzrosło z 4.32%
% dla UR1_ukl do 43.39% dla UR2_ukl. Taki wynik pokazuje, że optimum
% symetryczne poprawia szybkość początkowej reakcji układu, ale powoduje
% znacznie większą oscylacyjność odpowiedzi.

% c) Niewielkie pogorszenie odpowiedzi po dyskretyzacji optimum modułu -
% Dyskretyzacja układu regulacji z regulatorem optimum modułu nie zmieniła
% znacząco czasów odpowiedzi, ale zwiększyła przeregulowanie. Dla UR1_ukl
% czas narastania zmniejszył się z 0.0152s do 0.0150s, natomiast czas
% ustalania wzrósł z 0.0422s do 0.0430s, czyli o 1.92% (tab. 4).
% Przeregulowanie zwiększyło się z 4.32% do 5.80%, czyli o 34.18%.
% W praktyce dyskretyzacja w tym przypadku zachowała podobną szybkość
% odpowiedzi, ale pogorszyła tłumienie przebiegu przejściowego.

% d) Wzrost oscylacyjności po dyskretyzacji optimum symetrycznego -
% Dyskretyzacja układu regulacji z regulatorem optimum symetrycznego
% zwiększyła agresywność odpowiedzi. Dla UR2_ukl czas narastania zmniejszył
% się z 0.0106s do 0.0100s, a czas ustalania wzrósł z 0.0828s do 0.0840s
% (tab. 5). Przeregulowanie wzrosło z 43.39% do 48.26%, czyli o 11.22%.
% Wynik ten potwierdza, że przejście do domeny dyskretnej dodatkowo
% zwiększyło oscylacyjny charakter odpowiedzi regulatora symetrycznego.

% e) Szybsze narastanie kosztem gorszego tłumienia dla optimum
% symetrycznego w domenie dyskretnej - W domenie dyskretnej układ
% z regulatorem optimum symetrycznego szybciej rozpoczyna odpowiedź niż
% układ z regulatorem optimum modułu, ale ma znacznie większe
% przeregulowanie. Czas narastania zmniejszył się z 0.015s dla UR1_ukl_D
% do 0.010s dla UR2_ukl_D, czyli o 33.33% (tab. 6). Jednocześnie czas
% ustalania wzrósł z 0.043s do 0.084s, a przeregulowanie z 5.80%
% do 48.26%. Zależność ta wskazuje, że optimum symetryczne daje szybszą
% reakcję początkową, ale pogarsza stabilizację odpowiedzi.

% f) Zgodność metody tablicowej i metody układu zamkniętego dla optimum
% modułu - Dla optimum modułu metoda tablicowa dała takie same parametry
% jak metoda układu zamkniętego. Dla UR1_ukl_D oraz UR1_tab_D czas
% narastania wyniósł 0.015s, czas ustalania 0.043s, przeregulowanie
% 5.80%, a uchyb ustalony 3.1131e-13 (tab. 9). Potwierdza to, że dla
% optimum modułu metoda tablicowa bardzo dobrze odwzorowuje metodę
% układu zamkniętego i może być stosowana jako prostszy sposób doboru
% nastaw.

% g) Ograniczenie przeregulowania przez metodę tablicową dla optimum
% symetrycznego - Dla optimum symetrycznego metoda tablicowa zmniejszyła
% przeregulowanie względem metody układu zamkniętego. Dla UR2_ukl_D
% przeregulowanie wyniosło 48.26%, natomiast dla UR2_tab_D spadło
% do 41.85%, czyli o 13.28% (tab. 10). Czas ustalania również zmniejszył
% się z 0.084s do 0.082s, przy jednoczesnym wzroście czasu narastania
% z 0.010s do 0.011s. Można więc uznać, że metoda tablicowa dla optimum
% symetrycznego daje odpowiedź nieco wolniejszą na początku, ale lepiej
% tłumioną.

% h) Tłumiący wpływ ograniczenia prądowego dla regulatora optimum modułu
% metodą tablicową - W symulacji z regulatorem Reg1_tab_D_sym ograniczenie
% prądowe zmniejszyło przeregulowanie względem liniowego układu UR1_tab_D.
% Czas narastania wzrósł z 0.015s do 0.0165s, czyli o 9.72% (tab. 12).
% Jednocześnie czas ustalania zmniejszył się z 0.043s do 0.0212s,
% a przeregulowanie spadło z 5.80% do 1.89%, czyli o 67.36%.
% Wskazuje to, że ograniczenie prądowe działa w tym przypadku jak dodatkowe
% tłumienie, ograniczając oscylacje odpowiedzi.

% i) Pogorszenie szybkości narastania dla regulatora symetrycznego metodą
% tablicową po wprowadzeniu ograniczenia prądowego - W symulacji
% z regulatorem Reg2_tab_D_sym ograniczenie prądowe wyraźnie wydłużyło
% czas narastania względem układu UR2_tab_D. Czas narastania wzrósł
% z 0.011s do 0.0165s, czyli o 49.62% (tab. 13). Czas ustalania wzrósł
% z 0.082s do 0.0891s, natomiast przeregulowanie spadło z 41.85%
% do 37.58%. Taki przebieg pokazuje, że ograniczenie prądowe zmniejsza
% agresywność regulatora, ale odbywa się to kosztem wolniejszego początku
% odpowiedzi.

% j) Tłumiący wpływ ograniczenia prądowego dla regulatora optimum modułu
% metodą układu zamkniętego - W symulacji z regulatorem Reg1_ukl_D_sym
% ograniczenie prądowe spowodowało takie same zmiany jak dla metody
% tablicowej optimum modułu. Czas narastania wzrósł z 0.015s do 0.0165s,
% czyli o 9.72% (tab. 14). Czas ustalania zmniejszył się z 0.043s
% do 0.0212s, a przeregulowanie spadło z 5.80% do 1.89%, czyli o 67.36%.
% Na tej podstawie można stwierdzić, że dla optimum modułu sposób
% wyznaczenia regulatora nie wpłynął na zachowanie układu po wprowadzeniu
% ograniczenia prądowego.

% k) Pogorszenie szybkości narastania dla regulatora symetrycznego metodą
% układu zamkniętego po wprowadzeniu ograniczenia prądowego - W symulacji
% z regulatorem Reg2_ukl_D_sym ograniczenie prądowe spowodowało największy
% wzrost czasu narastania spośród analizowanych przypadków. Czas narastania
% wzrósł z 0.010s do 0.0165s, czyli o 64.59% (tab. 15). Czas ustalania
% zwiększył się z 0.084s do 0.0899s, natomiast przeregulowanie spadło
% z 48.26% do 44.41%. Dane te pokazują, że regulator symetryczny wyznaczony
% metodą układu zamkniętego jest szczególnie wrażliwy na ograniczenie
% sygnału sterującego.

% l) Wydłużenie czasu ustalania dla regulatora tablicowego optimum modułu
% po zastosowaniu anti-windup - W symulacji z regulatorem Reg1_tab_D_sym
% zastosowanie mechanizmu anti-windup nie zmieniło czasu narastania,
% który wyniósł 0.0165s zarówno bez anti-windup, jak i z anti-windup
% (tab. 21). Jednocześnie czas ustalania wzrósł z 0.0212s do 0.2022s,
% czyli o 853.51%. Przeregulowanie zostało natomiast całkowicie usunięte,
% ponieważ spadło z 1.89% do 0%. Dla tego regulatora anti-windup poprawia
% tłumienie odpowiedzi, ale znacząco wydłuża czas dojścia do stanu
% ustalonego.

% m) Skrócenie czasu ustalania dla regulatora tablicowego symetrycznego
% po zastosowaniu anti-windup - W symulacji z regulatorem Reg2_tab_D_sym
% mechanizm anti-windup poprawił czas ustalania oraz usunął przeregulowanie.
% Czas narastania pozostał bez zmian i wyniósł 0.0165s (tab. 22). Czas
% ustalania zmniejszył się z 0.0891s do 0.0388s, czyli o 56.50%,
% a przeregulowanie spadło z 37.58% do 0%. Wyniki te wskazują, że dla
% regulatora tablicowego symetrycznego anti-windup skutecznie ogranicza
% stan przejściowy.

% n) Wydłużenie czasu ustalania dla regulatora optimum modułu metodą
% układu zamkniętego po zastosowaniu anti-windup - W symulacji z regulatorem
% Reg1_ukl_D_sym zastosowanie anti-windup dało takie same efekty jak dla
% regulatora tablicowego optimum modułu. Czas narastania nie zmienił się
% i wyniósł 0.0165s (tab. 23). Czas ustalania wzrósł z 0.0212s do 0.2022s,
% czyli o 853.51%, natomiast przeregulowanie spadło z 1.89% do 0%.
% Pokazuje to, że w przypadku optimum modułu anti-windup eliminuje
% przeregulowanie, ale jednocześnie silnie spowalnia stabilizację odpowiedzi.

% o) Skrócenie czasu ustalania dla regulatora symetrycznego metodą układu
% zamkniętego po zastosowaniu anti-windup - W symulacji z regulatorem
% Reg2_ukl_D_sym mechanizm anti-windup wyraźnie poprawił przebieg
% przejściowy. Czas narastania pozostał bez zmian i wyniósł 0.0165s
% (tab. 24). Czas ustalania zmniejszył się z 0.0899s do 0.0356s,
% czyli o 60.37%, a przeregulowanie spadło z 44.41% do 0.37%.
% Takie wartości świadczą o bardzo skutecznym ograniczeniu oscylacji
% przy jednoczesnym skróceniu czasu stabilizacji.

% p) Brak wpływu anti-windup na czas narastania w symulacjach
% nieliniowych - Mechanizm anti-windup nie wpłynął na czas narastania
% w żadnym z analizowanych wariantów. We wszystkich porównaniach czas
% narastania wyniósł 0.0165s zarówno bez anti-windup, jak i z anti-windup
% (tab. 21-24). Z tego wynika, że działanie anti-windup ujawnia się głównie
% po początkowej fazie odpowiedzi, czyli w czasie ustalania oraz
% przeregulowaniu.

% q) Różny wpływ anti-windup w zależności od kryterium doboru regulatora -
% Mechanizm anti-windup działał korzystnie dla regulatorów symetrycznych,
% ale pogarszał czas ustalania regulatorów optimum modułu. Dla Reg1_tab_D_sym
% oraz Reg1_ukl_D_sym czas ustalania wzrósł o 853.51% (tab. 21 i tab. 23).
% Dla Reg2_tab_D_sym czas ustalania zmniejszył się o 56.50%, a dla
% Reg2_ukl_D_sym o 60.37% (tab. 22 i tab. 24). Wskazuje to, że skuteczność
% anti-windup zależy od nastaw regulatora oraz charakteru jego odpowiedzi.

% r) Redukcja przeregulowania przez anti-windup we wszystkich symulacjach -
% Zastosowanie anti-windup zmniejszyło przeregulowanie we wszystkich
% analizowanych wariantach. Dla Reg1_tab_D_sym i Reg1_ukl_D_sym
% przeregulowanie spadło z 1.89% do 0% (tab. 21 i tab. 23). Dla
% Reg2_tab_D_sym spadło z 37.58% do 0% (tab. 22), a dla Reg2_ukl_D_sym
% z 44.41% do 0.37% (tab. 24). Potwierdza to, że anti-windup skutecznie
% ogranicza skutki nasycenia regulatora i zmniejsza oscylacyjność odpowiedzi.

%% === Rekomendacje ===

% Na podstawie przeprowadzonej analizy zaleca się stosowanie mechanizmu
% anti-windup w układzie z ograniczeniem prądowym. Mechanizm ten jest
% szczególnie uzasadniony wtedy, gdy priorytetem jest ograniczenie
% przeregulowania, wygładzenie odpowiedzi oraz zmniejszenie wpływu nasycenia
% regulatora na przebieg stanu przejściowego.

% W kwestii czasu regulacji optymalnym wyborem jest regulator dobrany według
% kryterium symetrycznego z mechanizmem anti-windup. Wariant wyznaczony
% metodą układu zamkniętego można rekomendować, gdy najważniejszy jest
% krótki czas ustalania, natomiast wariant tablicowy jest dobrym wyborem,
% gdy istotna jest prostsza procedura doboru nastaw oraz małe
% przeregulowanie.
