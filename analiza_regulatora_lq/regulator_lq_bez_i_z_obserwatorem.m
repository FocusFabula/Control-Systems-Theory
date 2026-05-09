

function [S] = gen_UR(R, phi, vel, czy_dyskretny, Ts)
    s = tf('s');
    % disp('=== obiekt nominalny ======');

    G = 1/(s^2);                     % transmitancja obiektu
    [A,B,C,D] = ssdata(G);           % model stanu
    Gss = ss(G);                     % obiekt w przestrzeni stanów

    n = size(A,1);                   
    Q = eye(n) * phi;                % macierz wag Q

    K = lqr(A,B,Q,R);                % regulator statyczny LQR
    Acl = A - B*K;                   % macierz układu zamkniętego

    P = eig(Acl);                    % bieguny układu zamkniętego
    L = place(A',C',vel*P)';         % obserwator

    if czy_dyskretny
        % --- DYSKRETYZACJA OBIEKTU ---
        Gss_d = c2d(Gss, Ts, 'zoh');         % Dyskretyzacja modelu Gss (obiekt fizyczny)
        [Ad, Bd, Cd, Dd] = ssdata(Gss_d);

        % --- UKŁAD ZAMKNIĘTY IDEALNY (Dyskretny) ---
        Acl_d = Ad - Bd*K;                   % macierz układu zamkniętego (dyskretna)
        S.Gclz_q = ss(Acl_d, Bd, Cd, Dd, Ts); % układ zamknięty bez obserwatora
        S.F_q = 1/dcgain(S.Gclz_q);          % kompensacja wejścia
        S.Gclz_q = S.Gclz_q * S.F_q; 
        
        S.Gclu_q = ss(Acl_d, Bd, -K, 0, Ts) * S.F_q; % sygnał sterujący bez obserwatora

        % Dyskretyzacja regulatora dynamicznego
        Css = reg(Gss, K, L);                
        Css_d = c2d(Css, Ts, 'zoh');         

        % --- POŁĄCZENIE DYSKRETNEGO OBIEKTU Z DYSKRETNYM REGULATOREM ---
        S.Gclz_o = feedback(Gss_d, Css_d, 1);
        S.F_o = 1/dcgain(S.Gclz_o);          
        S.Gclz_o = S.Gclz_o * S.F_o;         
        
        S.Gclu_o = feedback(Gss_d*Css_d, 1, 1) * S.F_o; 
        
        S.t = (0:Ts:10)';                     % wektor czasu dla dyskretnego
    else
        % --- WERSJA CIĄGŁA ---
        S.Gclz_q = ss(Acl,B,C,0);
        S.F_q = 1/dcgain(S.Gclz_q);
        S.Gclz_q = S.Gclz_q * S.F_q;
        
        S.Gclu_q = ss(Acl,B,-K,0) * S.F_q;

        Css = reg(Gss,K,L);
        S.Gclz_o = feedback(Gss,Css,1);
        S.F_o = 1/dcgain(S.Gclz_o);
        S.Gclz_o = S.Gclz_o * S.F_o;
        
        S.Gclu_o = feedback(Gss*Css, 1, 1) * S.F_o;
        
        S.t = (0:0.01:10)';               
    end

    % --- OBLICZENIE UCHYBÓW ---
    [y_q, ~] = step(S.Gclz_q, S.t);
    [y_o, ~] = step(S.Gclz_o, S.t);

    S.e_q = 1 - y_q;                   
    S.e_o = 1 - y_o;
    
    % --- DODATKOWE INFO ---
    S.param.R = R;
    S.param.phi = phi;
    S.param.vel = vel;
    S.param.Ts = Ts;
    S.param.czy_dyskretny = czy_dyskretny;
end


function gen_wykresy(S_array, naglowek, wlasna_legenda)
    n = length(S_array);
    figure('Name', naglowek, 'NumberTitle', 'off');
    sgtitle(naglowek, 'FontSize', 14, 'FontWeight', 'bold');
    colors = get(gca, 'ColorOrder'); 
   
    
    % Definiujemy wspólny zakres częstotliwości do 10 rad/s
    w_range = logspace(-2, 1, 500); 

    for i = 1:n
        S = S_array(i);
        t = S.t;
        c_idx = mod(i-1, size(colors, 1)) + 1;
        current_color = colors(c_idx, :);
        is_dist = S.param.czy_dyskretny; 
        
        if is_dist; func = @stairs; else; func = @plot; end


        % --- 1. Odpowiedź skokowa y(t) ---
        subplot(3, 2, 1); hold on;
        [y_q, ~] = step(S.Gclz_q, t); [y_o, ~] = step(S.Gclz_o, t);
        func(t, squeeze(y_q), 'Color', current_color, 'LineWidth', 1.5);
        func(t, squeeze(y_o), 'Color', current_color, 'LineStyle', '--');
        title('Odpowiedź skokowa y(t)'); grid on;
        xlabel('czas [s]');
        ylabel('y(t)');

        % --- 2. Sygnał sterujący u(t) ---
        subplot(3, 2, 2); hold on;
        [u_q, ~] = step(S.Gclu_q, t); [u_o, ~] = step(S.Gclu_o, t);
        func(t, squeeze(u_q), 'Color', current_color);
        func(t, squeeze(u_o), 'Color', current_color, 'LineStyle', '--');
        title('Sygnał sterujący u(t)'); grid on;
        xlabel('czas [s]');
        ylabel('u(t)');

        % --- 3. Uchyb regulacji e(t) ---
        subplot(3, 2, 3); hold on;
        func(t, S.e_q, 'Color', current_color);
        func(t, S.e_o, 'Color', current_color, 'LineStyle', '--');
        title('Uchyb e(t)'); grid on;
        xlabel('czas [s]');
        ylabel('e(t)');
        
        % --- 4. Odpowiedź impulsowa h(t) ---
        subplot(3, 2, 4); hold on;
        [yi_q, ~] = impulse(S.Gclz_q, t); [yi_o, ~] = impulse(S.Gclz_o, t);
        func(t, squeeze(yi_q), 'Color', current_color);
        func(t, squeeze(yi_o), 'Color', current_color, 'LineStyle', '--');
        title('Odpowiedź Impulsowa h(t)'); grid on;
        xlabel('czas [s]');
        ylabel('h(t)');

        % --- 5. Magnituda (Bode) - Ideal vs Obs ---
        subplot(3, 2, 5); hold on;
        [mag_q, ~] = bode(S.Gclz_q, w_range);
        [mag_o, ~] = bode(S.Gclz_o, w_range);
        semilogx(w_range, 20*log10(squeeze(mag_q)), 'Color', current_color);
        semilogx(w_range, 20*log10(squeeze(mag_o)), 'Color', current_color, 'LineStyle', '--');
        title('Magnituda [dB]'); grid on; ylabel('dB'); xlim([min(w_range) max(w_range)]);
        xlabel('rad/s')

        % --- 6. Faza (Bode) - Ideal vs Obs ---
        subplot(3, 2, 6); hold on;
        [~, ph_q] = bode(S.Gclz_q, w_range);
        [~, ph_o] = bode(S.Gclz_o, w_range);
        semilogx(w_range, squeeze(ph_q), 'Color', current_color);
        semilogx(w_range, squeeze(ph_o), 'Color', current_color, 'LineStyle', '--');
        title('Faza [deg]'); grid on; ylabel('deg'); xlabel('rad/s'); xlim([min(w_range) max(w_range)]);
    end
    
    subplot(3, 2, 1);
    legend(wlasna_legenda, 'Location', 'best', 'FontSize', 7);
end


function [T_parametry] = gen_tabela_danych(ob, tytul)
    % ob - struktura z danymi (wynik gen_UR)
    % tytul - nazwa testu
    
    Gclz_q = ob.Gclz_q;
    Gclz_o = ob.Gclz_o;
    is_dist = ob.param.czy_dyskretny;

    % 1. Nagłówek raportu
    disp(' ');
    disp('***********************************************************');
    disp(['Parametry czasowe i bieguny Gclz_q oraz Gclz_0 dla : ', upper(tytul)]);
    if is_dist
        disp(['DOMENA: Dyskretna (Z), Ts = ', num2str(ob.param.Ts), ' s']);
    else
        disp('DOMENA: Ciągła (S)');
    end
    disp('***********************************************************');

    % % 2. WYŚWIETLENIE TRANSMITANCJI
    % disp('--- TRANSMITANCJE UKŁADU ZAMKNIĘTEGO (W-Y) ---');
    % disp('Układ Idealny Gclz_q:');
    % 
    % Gclz_q = tf(ob.Gclz_q) 
    % disp('Układ z Obserwatorem Gclz_o:');
    % Gclz_o =tf(Gclz_o)    
    % disp('-----------------------------------------------------------');

    % 3. Parametry czasowe
    iq = stepinfo(Gclz_q);
    io = stepinfo(Gclz_o);
    ess_q = abs(1 - dcgain(Gclz_q));
    ess_o = abs(1 - dcgain(Gclz_o));

    Metryka = {'Czas narastania (t_r)'; 'Czas ustalania (t_us)'; 'Przeregulowanie [%]'; 'Uchyb ustalony'};
    T_parametry = table(Metryka, ...
        [iq.RiseTime; iq.SettlingTime; iq.Overshoot; ess_q], ...
        [io.RiseTime; io.SettlingTime; io.Overshoot; ess_o], ...
        'VariableNames', {'Metryka', 'Gclz_q', 'Gclz_o'});

    disp('--- PARAMETRY CZASOWE ---');
    disp(T_parametry);

    % % 4. Analiza Stabilności (Bieguny)
    % disp('--- ANALIZA BIEGUNÓW (Stabilność) ---');
    % p_q = pole(Gclz_q);
    % p_o = pole(Gclz_o);
    % 
    % if is_dist
    %     T_p_o = table(p_o, abs(p_o), 'VariableNames', {'Bieguny_Obs', 'Modul_abs_Z'});
    %     stabilny = all(abs(p_o) < 1.0);
    % else
    %     T_p_o = table(p_o, real(p_o), 'VariableNames', {'Bieguny_Obs', 'Czesc_Real_S'});
    %     stabilny = all(real(p_o) < 0);
    % end
    % disp(T_p_o);
    % 
    % % 5. Analiza Zer
    % disp('--- ANALIZA ZER ---');
    % z_o = zero(Gclz_o);
    % if isempty(z_o)
    %     disp('Układ nie posiada zer.');
    % else
    %     if is_dist
    %         T_z = table(z_o, abs(z_o), 'VariableNames', {'Zera_Obs', 'Modul_abs_Z'});
    %     else
    %         T_z = table(z_o, real(z_o), 'VariableNames', {'Zera_Obs', 'Czesc_Real_S'});
    %     end
    %     disp(T_z);
    % end
end


function [T_delta] = gen_tabela_porownawcza(tab1, tab2, nazwa1, nazwa2, tryb)
    % tryb: 
    % 1 - Porównanie dwóch różnych badań (np. R1 vs R01)
    % 2 - Porównanie Ideal vs Obserwator (tylko dla tab1)
    
    Metryka = tab1.Metryka;
    eps_val = 1e-10;

    if tryb == 1
        % --- PORÓWNANIE MIĘDZY BADANIAMI ---
        V1_Ideal = tab1.Gclz_q;
        V2_Ideal = tab2.Gclz_q;
        V1_Obs = tab1.Gclz_o;
        V2_Obs = tab2.Gclz_o;

        Delta_Ideal = ((V2_Ideal - V1_Ideal) ./ (abs(V1_Ideal) + eps_val)) * 100;
        Delta_Obs = ((V2_Obs - V1_Obs) ./ (abs(V1_Obs) + eps_val)) * 100;

        T_delta = table(Metryka, Delta_Ideal, Delta_Obs);
        T_delta.Properties.VariableNames{2} = 'Zmiana_Ideal_proc';
        T_delta.Properties.VariableNames{3} = 'Zmiana_Obs_proc';
        
        naglowek = ['PORÓWNANIE PARAMETRÓW: ', nazwa1, ' vs ', nazwa2];

    elseif tryb == 2
        % --- PORÓWNANIE WEWNĘTRZNE (Ideal vs Obserwator) ---
        Wartosc_Ideal = tab1.Gclz_q;
        Wartosc_Obs = tab1.Gclz_o;

        % Liczymy o ile procent obserwator różni się od ideału
        Blad_procentowy = ((Wartosc_Obs - Wartosc_Ideal) ./ (abs(Wartosc_Ideal) + eps_val)) * 100;

        % Tworzymy tabelę z czasami ORAZ procentem
        T_delta = table(Metryka, Wartosc_Ideal, Wartosc_Obs, Blad_procentowy);
        
        % Czytelne nazwy kolumn
        T_delta.Properties.VariableNames{2} = 'Ideal_LQR';
        T_delta.Properties.VariableNames{3} = 'Z_Obserwatorem';
        T_delta.Properties.VariableNames{4} = 'Roznica_proc';
        
        naglowek = ['ANALIZA WPŁYWU OBSERWATORA: ', nazwa1];
    end
    disp(' ');
    disp('================================================================');
    disp(['ANALIZA JAKOŚCIOWA: ', nazwa1, ' vs ', nazwa2]);
    disp(T_delta);
    disp('================================================================');
end



% === Przyjęte Oznaczenie ===
% parametr odpowiadający za pomnożenie biegunów obserwatora w tym
% sprawozdaniu został oznaczony Vs

% === Sprawdzenie wpływu R na układ regulacji przy Vs =3===
% Do sprawdzenia przyjęto wartości R = [1, 0.1, 0.01] oraz dla szybkości
% obserwatora Vs = 3
[UR_R1] = gen_UR(1, 1, 3, false, 0); % w dalszej części sprawozdania te układy regulacji będą traktowane jako układy odniesienia
[UR_R01] = gen_UR(0.1, 1, 3, false, 0);
[UR_R001] = gen_UR(0.01, 1, 3, false, 0);

% === Wygenerowanie odpowiednich wykresów
gen_wykresy([UR_R1, UR_R01, UR_R001],'Wykresy pokazujące wpływ R =[1, 0.1, 0.01] Vs=3', [{'Gclz_q R=1', 'Gclz_o R=1'}, {'Gclz_q R=0.1', 'Gclz_o R=0.1'}, {'Gclz_q R=0.01', 'Gclz_o R=0.01'}])

% === Wygenerowanie odpowiednich danych do poźniejszej analizy ===
UR_R1_param = gen_tabela_danych(UR_R1, 'Tab 1.  R = 1 Q=1 Vs=3');
UR_R01_param = gen_tabela_danych(UR_R01, 'Tab 2.R = 0.1 Q=1 Vs=3');
UR_R001_param = gen_tabela_danych(UR_R001, 'Tab3. R = 0.01 Q=1 Vs=3');


% === Sprawdzenie wplywu Q na odpowiedzi czasowe układu regulacji przy Vs =3 === 
% Do sprawdzenia przyjęto wartości Q=[1, 10, 100] oraz szybkość obserwatora
% równą Vs = 3 
[UR_Q1] = gen_UR(1, 1, 3, false, 0);
[UR_Q10] = gen_UR(1, 10, 3, false, 0);
[UR_Q100] = gen_UR(1, 100, 3, false, 0);

% === Wygenerowanie odpowiednich wykresów
gen_wykresy([UR_Q1, UR_Q10, UR_Q100],'Wykresy pokazujące wpływ Q =[1, 10, 100] Vs=3', [{'Gclz_q Q=1', 'Gclz_o Q=1'}; {'Gclz_q Q=10', 'Gclz_o Q=10'}; {'Gclz_q Q=100', 'Gclz_o Q=100'}])

% === Wygenerowanie odpowiednich danych do poźniejszej analizy ===
UR_Q1_param = gen_tabela_danych(UR_Q1, 'Tab 4. R = 1 Q=1 Vs=3');
UR_Q10_param = gen_tabela_danych(UR_Q10, 'Tab 5. R = 1 Q=1 Vs=3');
UR_Q100_param = gen_tabela_danych(UR_Q100, 'Tab 6. R = 1 Q=100 Vs=3');


% == Sprawdzenie odpowiedzi układu regulaci względem Vs =3 i Vs =5 ===
% Do sprawdzenia przyjęto wartości Vs =[3, 5] oraz R = 1 oraz Q = 1 
[UR_VS3] = gen_UR(1, 1, 3, false, 0); % w dalszej części sprawozdania te układy regulacji będą traktowane jako układy odniesienia
[UR_VS5] = gen_UR(1, 1, 5, false, 0);

% === Wygenerowanie odpowiednich wykresów ===
gen_wykresy([UR_VS3, UR_VS5],'Wykresy pokazujące wpływ Vs=[3, 5] przy R=1 oraz Q=1 ', [{'Gclz_q VS=3', 'Gclz_o Vs=3'}, {'Gclz_q VS=5', 'Gclz_o Vs=5'}])


% === Wygenerowanie odpowiednich danych do poźniejszej analizy ===
UR_VS3_param = gen_tabela_danych(UR_VS3, 'Tab 7. R=1 Q=1 Vs=3');
UR_VS5_param = gen_tabela_danych(UR_VS5, 'Tab 8. R=1 Q=1 Vs=5');

disp('=== Tabele porównawcze Q=? R=? Vs=?')

% === Wygenerowanie tabel porównawczych Względem R=? Q=1 i Vs=3 ===
gen_tabela_porownawcza(UR_VS3_param, UR_R01_param, 'Tab 9. R=1', 'R=0.1', 1); 
gen_tabela_porownawcza(UR_VS3_param, UR_R001_param, 'Tab 10. R=1', 'R=0.01', 1);


% === Wygenerowanie tabel porównawczych Względem  Q=? R=1 i Vs=3 ===
gen_tabela_porownawcza(UR_VS3_param, UR_Q10_param, 'Tab 11 Q=1', 'Q=10', 1); 
gen_tabela_porownawcza(UR_VS3_param, UR_Q100_param, 'Tab 12. Q=1', 'Q=100', 1);


% === Wygenerowanie tabel porównawczych układu bez obserwatora i z obserwatorem ===
gen_tabela_porownawcza(UR_VS3_param, UR_VS3_param , 'Tab 12. bez_obs', 'z_obs Vs=3', 2);
gen_tabela_porownawcza(UR_VS5_param, UR_VS5_param , 'TAb 13. bez_obs', 'z_obs Vs=5', 2);

%

disp('=== Tabele zawierające dane odnośnie układów wyregulowanego ciągłego')

% === Regulacja regulatora ===
% Wymagania regulacji tr -0.5s ~> tr tus - 0.4  <= t_us  uchyb ~0 z obserwatorem 
[UR_ureg] = gen_UR(0.1, 5, 3, false, 0); 
UR_ureg_param = gen_tabela_danych(UR_ureg, 'Tab 13. R=0.1 Q=5 Vs=3');

% === wygenerowanie tabeli porównawczej względem układu odniesienia R=1 Q=1
% Vs=1
gen_tabela_porownawcza(UR_R01_param, UR_ureg_param, 'Tab 14. Do układu odniesienia', 'UR_ureg R=1, Q=5 Vs=3', 1);

[UR_ureg5] = gen_UR(0.1, 5, 5, false, 0); 
UR_ureg5_param = gen_tabela_danych(UR_ureg5, 'Tab 15. R=0.1 Q=5 Vs=5');


%na podstawie tabeli nr 15 przyjęto że dalej w sprawozdanie Vs=5

% === wygenerowanie tabeli porównawczej względem układu odniesienia R=1 Q=1
gen_tabela_porownawcza(UR_R01_param, UR_ureg5_param, 'Tab 16. do układu odniesienia','UR_ureg5', 1);

% === Wygenerowanie tabel porównawczych układu bez obserwatora i z obserwatorem ===
gen_tabela_porownawcza(UR_ureg5_param, UR_ureg5_param , 'Tab 17. dla R=0.1 oraz Q=5 bez_obs', 'z_obs Vs=5', 2);



%% Dla wartości R=0.1 Q=5 Vs=5 osiągnięto wymagania przyjęte regulacji
gen_wykresy(UR_ureg5, 'Tab 18. Wyregulowany układ dla wartości R=0.1 Q=5 oraz Vs=5', [{'Gclz_q', 'Gclz_o'}]);

% % === Przeprowadzenie dyskretyzacji ===

%% === Dobór parametrów dyskretyzacji ===
% Czas próbkowania Tp wyznaczono w oparciu o charakterystykę częstotliwościową układu.
% 1. Z wykresu Bodego odczytano pulsację graniczną (pasmo -3 dB): w_b ≈ 1 rad/s.
% 2. Częstotliwość graniczna pasma: f_b = w_b / 2π ≈ 0,16 Hz.
% 3. Przyjęto inżynierskie kryterium częstotliwości próbkowania f_p = 20 * f_b.
% 4. Obliczone parametry: f_p = 3,2 Hz => Tp = 1/f_p ≈ 0,3125 s.
% Przyjęto ostatecznie: Tp = 0.3 s.

% === Wygenerowanie układu wyregulowanego ===
[UR_ureg5_d] = gen_UR(0.1, 5, 5, true, 0.3);
gen_wykresy(UR_ureg5_d, 'Tab 19. Wyregulowany układ dyskretny dla wartości Tp=0.32', [{'Gclz_q', 'Gclz_o'}]);
UR_ureg5_d_param = gen_tabela_danych(UR_ureg5_d, 'UR_ureg5_d_new dla tp=0.32');

% Na przedstawionym wykresie zaobserwowano bardzo mocne oscylacje jak i
% błąd dyskretyzacji wynika to z tego że przyjęta metoda aproksymacji Tp
% nie bierze pod uwagę dodatkowej dynamiki obserwatora żeby wziąc dynamikę
% pod uwagę za aproksymowany Tp przyjęto Tp/5= 0.06s

[UR_ureg5_d_new] = gen_UR(0.1, 5, 5, true, 0.06);
gen_wykresy(UR_ureg5_d_new, 'Wyregulowany układ dyskretny dla wartości Tp=0.06', [{'Gclz_q', 'Gclz_o'}]);
UR_ureg5_d_new_param = gen_tabela_danych(UR_ureg5_d_new, 'Tab 21 UR_ureg5_d_new dla Tp=0.06');

% Porównanie wyregulowanego modelu ciągłego UR_ureg5 do modelu dyskretnego
% UR_ureg5_d dla Tp=0.06
gen_wykresy([UR_ureg5,UR_ureg5_d_new], 'Wyregulowany układ dyskretny oraz układ liniowy', [{'Gclz_q', 'Gclz_o'}, {'Gclz_qd', 'Gclz_od'}]);
gen_tabela_porownawcza(UR_ureg5_param, UR_ureg5_d_new_param , 'Tab 22. UR_ureg5_d_param', 'UR_ureg5_d_new', 1);

% W przypadku wykonania dyskretyzacji na układnie nie posiadającym
% obserwatora układ ma bardzo zbliżone parametry czasowe jednak w przypadku
% układy posiadajacego obserwator układ jest 23% szybszy żeby zbliżyć układ
% z obserwatora do takich samych parametrów czasowych zmniejszono dynamikę
% obserwatora do Vs=3 oraz Q=2

[UR_ureg5_d_new] = gen_UR(0.1, 2, 3, true, 0.06);
gen_wykresy([UR_ureg5,UR_ureg5_d_new], 'Wyregulowany układ dyskretny oraz układ liniowy', [{'Gclz_q', 'Gclz_o'}, {'Gclz_qd', 'Gclz_od'}]);
UR_ureg5_d_new_param = gen_tabela_danych(UR_ureg5_d_new, 'Tab 22. UR_ureg5_d_new dla Tp=0.06');
gen_tabela_porownawcza(UR_ureg5_param, UR_ureg5_d_new_param , 'Tab 23. UR_ureg5_d_param', 'UR_ureg5_d_new', 1);

% Dla wyregulowanego układu z obserwatorem R=0.1 Q=2 i Vs=3 parametry
% czasowe są o 5% szybsze od wyregulowanego układu.

%ze względu na to że Dynamika układu obserwatora się zmieniła zmieniono
%czas próbkowania Tp /3 = 0.32s /3s. ze względu na zmianę szybkości
%próbkowania różnica pomiędzy paametrami czasowymi zwiększyła się do -21%
% więc żeby zbliżyć te parametry zmieniono Q= 1

[UR_ureg5_d_new] = gen_UR(0.1, 0.2, 3, true, 0.1);
gen_wykresy([UR_ureg5,UR_ureg5_d_new], 'Wyregulowany układ dyskretny oraz układ liniowy', [{'Gclz_q', 'Gclz_o'}, {'Gclz_qd', 'Gclz_od'}]);
UR_ureg5_d_new_param = gen_tabela_danych(UR_ureg5_d_new, 'Tab 24. UR_ureg5_d_new dla Tp=0.1');
gen_tabela_porownawcza(UR_ureg5_param, UR_ureg5_d_new_param , 'Tab 25. UR_ureg5_d_param', 'UR_ureg5_d_new', 1);

% Gotowe

%% === Wnioski ilościowe ===
% ***********************************************************
% Parametry czasowe i bieguny Gclz_q oraz Gclz_0 dla : TAB 1.  R = 1 Q=1 VS=3
% DOMENA: Ciągła (S)
% ***********************************************************
% --- PARAMETRY CZASOWE ---
%              Metryka               Gclz_q        Gclz_o  
%     _________________________    __________    __________
% 
%     {'Czas narastania (t_r)'}        2.7344        2.8394
%     {'Czas ustalania (t_us)'}        4.3456        4.6108
%     {'Przeregulowanie [%]'  }       0.43333       0.41593
%     {'Uchyb ustalony'       }    1.1102e-16    4.4409e-16
% 
% 
% ***********************************************************
% Parametry czasowe i bieguny Gclz_q oraz Gclz_0 dla : TAB 2.R = 0.1 Q=1 VS=3
% DOMENA: Ciągła (S)
% ***********************************************************
% --- PARAMETRY CZASOWE ---
%              Metryka             Gclz_q      Gclz_o  
%     _________________________    ______    __________
% 
%     {'Czas narastania (t_r)'}    2.2754         2.429
%     {'Czas ustalania (t_us)'}    4.1001        4.3919
%     {'Przeregulowanie [%]'  }         0             0
%     {'Uchyb ustalony'       }         0    2.2204e-16
% 
% 
% ***********************************************************
% Parametry czasowe i bieguny Gclz_q oraz Gclz_0 dla : TAB3. R = 0.01 Q=1 VS=3
% DOMENA: Ciągła (S)
% ***********************************************************
% --- PARAMETRY CZASOWE ---
%              Metryka             Gclz_q      Gclz_o  
%     _________________________    ______    __________
% 
%     {'Czas narastania (t_r)'}    2.2042        2.3892
%     {'Czas ustalania (t_us)'}    3.9982        4.3379
%     {'Przeregulowanie [%]'  }         0             0
%     {'Uchyb ustalony'       }         0    1.1102e-16
% 
% 
% ***********************************************************
% Parametry czasowe i bieguny Gclz_q oraz Gclz_0 dla : TAB 4. R = 1 Q=1 VS=3
% DOMENA: Ciągła (S)
% ***********************************************************
% --- PARAMETRY CZASOWE ---
%              Metryka               Gclz_q        Gclz_o  
%     _________________________    __________    __________
% 
%     {'Czas narastania (t_r)'}        2.7344        2.8394
%     {'Czas ustalania (t_us)'}        4.3456        4.6108
%     {'Przeregulowanie [%]'  }       0.43333       0.41593
%     {'Uchyb ustalony'       }    1.1102e-16    4.4409e-16
% 
% 
% ***********************************************************
% Parametry czasowe i bieguny Gclz_q oraz Gclz_0 dla : TAB 5. R = 0.1 Q=1 VS=3
% DOMENA: Ciągła (S)
% ***********************************************************
% --- PARAMETRY CZASOWE ---
%              Metryka             Gclz_q      Gclz_o  
%     _________________________    ______    __________
% 
%     {'Czas narastania (t_r)'}    2.2754         2.429
%     {'Czas ustalania (t_us)'}    4.1001        4.3919
%     {'Przeregulowanie [%]'  }         0             0
%     {'Uchyb ustalony'       }         0    2.2204e-16
% 
% 
% ***********************************************************
% Parametry czasowe i bieguny Gclz_q oraz Gclz_0 dla : TAB 6. R = 0.01 Q=1 VS=3
% DOMENA: Ciągła (S)
% ***********************************************************
% --- PARAMETRY CZASOWE ---
%              Metryka             Gclz_q    Gclz_o
%     _________________________    ______    ______
% 
%     {'Czas narastania (t_r)'}    2.2042    2.3892
%     {'Czas ustalania (t_us)'}    3.9982    4.3379
%     {'Przeregulowanie [%]'  }         0         0
%     {'Uchyb ustalony'       }         0         0
% 
% 
% ***********************************************************
% Parametry czasowe i bieguny Gclz_q oraz Gclz_0 dla : TAB 7. R=1 Q=1 VS=3
% DOMENA: Ciągła (S)
% ***********************************************************
% --- PARAMETRY CZASOWE ---
%              Metryka               Gclz_q        Gclz_o  
%     _________________________    __________    __________
% 
%     {'Czas narastania (t_r)'}        2.7344        2.8394
%     {'Czas ustalania (t_us)'}        4.3456        4.6108
%     {'Przeregulowanie [%]'  }       0.43333       0.41593
%     {'Uchyb ustalony'       }    1.1102e-16    4.4409e-16
% 
% 
% ***********************************************************
% Parametry czasowe i bieguny Gclz_q oraz Gclz_0 dla : TAB 8. R=1 Q=1 VS=5
% DOMENA: Ciągła (S)
% ***********************************************************
% --- PARAMETRY CZASOWE ---
%              Metryka               Gclz_q        Gclz_o  
%     _________________________    __________    __________
% 
%     {'Czas narastania (t_r)'}        2.7344        2.7637
%     {'Czas ustalania (t_us)'}        4.3456         4.454
%     {'Przeregulowanie [%]'  }       0.43333       0.42764
%     {'Uchyb ustalony'       }    1.1102e-16    2.2204e-16
% 
% === Tabele porównawcze Q=? R=? Vs=?
% 
% ================================================================
% ANALIZA JAKOŚCIOWA: Tab 9. R=1 vs R=0.1
%              Metryka             Zmiana_Ideal_proc    Zmiana_Obs_proc
%     _________________________    _________________    _______________
% 
%     {'Czas narastania (t_r)'}           -16.785             -14.451  
%     {'Czas ustalania (t_us)'}           -5.6492             -4.7486  
%     {'Przeregulowanie [%]'  }              -100                -100  
%     {'Uchyb ustalony'       }       -0.00011102         -0.00022204  
% 
% ================================================================
% 
% ================================================================
% ANALIZA JAKOŚCIOWA: Tab 10. R=1 vs R=0.01
%              Metryka             Zmiana_Ideal_proc    Zmiana_Obs_proc
%     _________________________    _________________    _______________
% 
%     {'Czas narastania (t_r)'}           -19.389             -15.855  
%     {'Czas ustalania (t_us)'}           -7.9942             -5.9194  
%     {'Przeregulowanie [%]'  }              -100                -100  
%     {'Uchyb ustalony'       }       -0.00011102         -0.00033307  
% 
% ================================================================
% 
% ================================================================
% ANALIZA JAKOŚCIOWA: Tab 11 Q=1 vs Q=10
%              Metryka             Zmiana_Ideal_proc    Zmiana_Obs_proc
%     _________________________    _________________    _______________
% 
%     {'Czas narastania (t_r)'}           -16.785             -14.451  
%     {'Czas ustalania (t_us)'}           -5.6492             -4.7486  
%     {'Przeregulowanie [%]'  }              -100                -100  
%     {'Uchyb ustalony'       }       -0.00011102         -0.00022204  
% 
% ================================================================
% 
% ================================================================
% ANALIZA JAKOŚCIOWA: Tab 12. Q=1 vs Q=100
%              Metryka             Zmiana_Ideal_proc    Zmiana_Obs_proc
%     _________________________    _________________    _______________
% 
%     {'Czas narastania (t_r)'}           -19.389             -15.855  
%     {'Czas ustalania (t_us)'}           -7.9942             -5.9194  
%     {'Przeregulowanie [%]'  }              -100                -100  
%     {'Uchyb ustalony'       }       -0.00011102         -0.00044409  
% 
% ================================================================
% 
% ================================================================
% ANALIZA JAKOŚCIOWA: Tab 12. bez_obs vs z_obs Vs=3
%              Metryka             Ideal_LQR     Z_Obserwatorem    Roznica_proc
%     _________________________    __________    ______________    ____________
% 
%     {'Czas narastania (t_r)'}        2.7344          2.8394           3.8392 
%     {'Czas ustalania (t_us)'}        4.3456          4.6108           6.1041 
%     {'Przeregulowanie [%]'  }       0.43333         0.41593          -4.0147 
%     {'Uchyb ustalony'       }    1.1102e-16      4.4409e-16       0.00033307 
% 
% ================================================================
% 
% ================================================================
% ANALIZA JAKOŚCIOWA: TAb 13. bez_obs vs z_obs Vs=5
%              Metryka             Ideal_LQR     Z_Obserwatorem    Roznica_proc
%     _________________________    __________    ______________    ____________
% 
%     {'Czas narastania (t_r)'}        2.7344          2.7637           1.0709 
%     {'Czas ustalania (t_us)'}        4.3456           4.454           2.4955 
%     {'Przeregulowanie [%]'  }       0.43333         0.42764          -1.3133 
%     {'Uchyb ustalony'       }    1.1102e-16      2.2204e-16       0.00011102 
% 
% ================================================================
% === Tabele zawierające dane odnośnie układów wyregulowanego ciągłego
% 
% ***********************************************************
% Parametry czasowe i bieguny Gclz_q oraz Gclz_0 dla : TAB 13. R=0.1 Q=5 VS=3
% DOMENA: Ciągła (S)
% ***********************************************************
% --- PARAMETRY CZASOWE ---
%              Metryka             Gclz_q      Gclz_o  
%     _________________________    ______    __________
% 
%     {'Czas narastania (t_r)'}    2.2139        2.3915
%     {'Czas ustalania (t_us)'}    4.0262        4.3508
%     {'Przeregulowanie [%]'  }         0             0
%     {'Uchyb ustalony'       }         0    2.2204e-16
% 
% 
% ================================================================
% ANALIZA JAKOŚCIOWA: Tab 14. Do układu odniesienia vs UR_ureg R=1, Q=5 Vs=3
%              Metryka             Zmiana_Ideal_proc    Zmiana_Obs_proc
%     _________________________    _________________    _______________
% 
%     {'Czas narastania (t_r)'}         -2.7047             -1.5454    
%     {'Czas ustalania (t_us)'}         -1.8017            -0.93575    
%     {'Przeregulowanie [%]'  }               0                   0    
%     {'Uchyb ustalony'       }               0                   0    
% 
% ================================================================
% 
% ***********************************************************
% Parametry czasowe i bieguny Gclz_q oraz Gclz_0 dla : TAB 15. R=0.1 Q=5 Q=10 VS=5
% DOMENA: Ciągła (S)
% ***********************************************************
% --- PARAMETRY CZASOWE ---
%              Metryka             Gclz_q    Gclz_o
%     _________________________    ______    ______
% 
%     {'Czas narastania (t_r)'}    2.2139    2.2784
%     {'Czas ustalania (t_us)'}    4.0262    4.1781
%     {'Przeregulowanie [%]'  }         0         0
%     {'Uchyb ustalony'       }         0         0
% 
% 
% ================================================================
% ANALIZA JAKOŚCIOWA: Tab 16. do układu odniesienia vs UR_ureg5
%              Metryka             Zmiana_Ideal_proc    Zmiana_Obs_proc
%     _________________________    _________________    _______________
% 
%     {'Czas narastania (t_r)'}         -2.7047                -6.199  
%     {'Czas ustalania (t_us)'}         -1.8017               -4.8665  
%     {'Przeregulowanie [%]'  }               0                     0  
%     {'Uchyb ustalony'       }               0           -0.00022204  
% 
% ================================================================
% 
% ================================================================
% ANALIZA JAKOŚCIOWA: Tab 17. dla R=0.1 oraz Q=5 bez_obs vs z_obs Vs=5
%              Metryka             Ideal_LQR    Z_Obserwatorem    Roznica_proc
%     _________________________    _________    ______________    ____________
% 
%     {'Czas narastania (t_r)'}     2.2139          2.2784           2.9168   
%     {'Czas ustalania (t_us)'}     4.0262          4.1781           3.7739   
%     {'Przeregulowanie [%]'  }          0               0                0   
%     {'Uchyb ustalony'       }          0               0                0   
% 
% ================================================================
% 
% ***********************************************************
% Parametry czasowe i bieguny Gclz_q oraz Gclz_0 dla : UR_UREG5_D_NEW DLA TP=0.32
% DOMENA: Dyskretna (Z), Ts = 0.3 s
% ***********************************************************

% --- PARAMETRY CZASOWE ---
%              Metryka             Gclz_q      Gclz_o  
%     _________________________    ______    __________
% 
%     {'Czas narastania (t_r)'}     NaN             NaN
%     {'Czas ustalania (t_us)'}     NaN             NaN
%     {'Przeregulowanie [%]'  }     NaN             NaN
%     {'Uchyb ustalony'       }       0      2.2204e-16
% 
% 
% ***********************************************************
% Parametry czasowe i bieguny Gclz_q oraz Gclz_0 dla : UR_UREG5_D_NEW DLA TP=0.06
% DOMENA: Dyskretna (Z), Ts = 0.06 s
% ***********************************************************
% --- PARAMETRY CZASOWE ---
%              Metryka             Gclz_q      Gclz_o  
%     _________________________    ______    __________
% 
%     {'Czas narastania (t_r)'}     2.22           1.74
%     {'Czas ustalania (t_us)'}     4.02            3.3
%     {'Przeregulowanie [%]'  }        0              0
%     {'Uchyb ustalony'       }        0     1.1102e-16
% 
% 
% ================================================================
% ANALIZA JAKOŚCIOWA: Tab 22. UR_ureg5_d_param vs UR_ureg5_d_new
%              Metryka             Zmiana_Ideal_proc    Zmiana_Obs_proc
%     _________________________    _________________    _______________
% 
%     {'Czas narastania (t_r)'}         0.27662              -23.632   
%     {'Czas ustalania (t_us)'}        -0.15407              -21.018   
%     {'Przeregulowanie [%]'  }               0                    0   
%     {'Uchyb ustalony'       }               0           0.00011102   
% 
% ================================================================
% 
% ***********************************************************
% Parametry czasowe i bieguny Gclz_q oraz Gclz_0 dla : TAB 22. UR_UREG5_D_NEW DLA TP=0.06
% DOMENA: Dyskretna (Z), Ts = 0.06 s
% ***********************************************************
% --- PARAMETRY CZASOWE ---
%              Metryka               Gclz_q        Gclz_o  
%     _________________________    __________    __________
% 
%     {'Czas narastania (t_r)'}          2.22          2.16
%     {'Czas ustalania (t_us)'}          4.08          4.08
%     {'Przeregulowanie [%]'  }             0             0
%     {'Uchyb ustalony'       }    1.1102e-16    1.1102e-16
% 
% 
% ================================================================
% ANALIZA JAKOŚCIOWA: Tab 23. UR_ureg5_d_param vs UR_ureg5_d_new
%              Metryka             Zmiana_Ideal_proc    Zmiana_Obs_proc
%     _________________________    _________________    _______________
% 
%     {'Czas narastania (t_r)'}          0.27662             -5.1987   
%     {'Czas ustalania (t_us)'}           1.3362             -2.3491   
%     {'Przeregulowanie [%]'  }                0                   0   
%     {'Uchyb ustalony'       }       0.00011102          0.00011102   
% 
% ================================================================
% 
% ***********************************************************
% Parametry czasowe i bieguny Gclz_q oraz Gclz_0 dla : TAB 24. UR_UREG5_D_NEW DLA TP=0.06
% DOMENA: Dyskretna (Z), Ts = 0.1 s
% ***********************************************************
% --- PARAMETRY CZASOWE ---
%              Metryka              Gclz_q       Gclz_o  
%     _________________________    ________    __________
% 
%     {'Czas narastania (t_r)'}         2.4           2.2
%     {'Czas ustalania (t_us)'}         4.2             4
%     {'Przeregulowanie [%]'  }    0.009084             0
%     {'Uchyb ustalony'       }           0    3.3307e-16
% 
% 
% ================================================================
% ANALIZA JAKOŚCIOWA: Tab 25. UR_ureg5_d_param vs UR_ureg5_d_new
%              Metryka             Zmiana_Ideal_proc    Zmiana_Obs_proc
%     _________________________    _________________    _______________
% 
%     {'Czas narastania (t_r)'}           8.4072             -3.4431   
%     {'Czas ustalania (t_us)'}           4.3166             -4.2638   
%     {'Przeregulowanie [%]'  }        9.084e+09                   0   
%     {'Uchyb ustalony'       }                0          0.00033307   
% ================================================================


%% === Wnioski ilościowe ===
% a) Wpływ wartości R na parametry czasowe układu regulacji bez obserwatora - Wraz ze
% zmniejszaniem wartości R=1 poprzez 0.1 do wartości 0.01 odnotowano spadek
% czasu ustalania z wartości 2.73s (tab 2) do odpowiednio 2.27s (tab 2)
% oraz 2.2s (tab 3) jest to spadek o 16.&% (tab 9) dla układu gdzie R=0.1
% oraz 19.4% (tab 10) w przypadku układu regulacji gdzie R=0.01. Wynika z
% tego że wraz ze spadkiem wartości R układ zwiększa swoją dynamikę. 
%
% b) Wpływ wartości R na parametry czasowe układu regulacji z obserwatorem - 
% Wraz ze zmniejszaniem wartości R=1 poprzez 0.1 do wartości 0.01 odnotowano spadek
% czasu ustalania z wartości 2.83s (tab 2) do odpowiednio 2.42s (tab 2)
% oraz 2.4s (tab 3) jest to spadek o 14.4% (tab 9) dla układu gdzie R=0.1
% oraz 15.9% (tab 10) w przypadku układu regulacji gdzie R=0.01. Wynika z
% tego że wraz ze spadkiem wartości R układ zwiększa swoją dynamikę. 
% 
% c) Wpływ wartości Q na parametry czasowe układu regulacji bez obserwatora - 
% Wraz ze zwiększaniem wartości Q=1 poprzez 10 do wartości 100 odnotowano
% spadek
% czasu ustalania z wartości 2.73s (tab 2) do odpowiednio 2.27s (tab 5)
% oraz 2.2s (tab 6) jest to spadek o 16.&% (tab 10) dlaukładu gdzie Q=10
% oraz 19.4% (tab 12) w przypadku układu regulacji gdzie Q=100. Wynika z
% tego że wraz ze wzrostem wartości Q układ zwiększa swoją dynamikę. 

% d) Wpływ wartości Q na parametry czasowe układu regulacji z obserwatorem - 
% Wraz ze zwiększaniem wartości Q=1 poprzez 10 do wartości 100 odnotowano spadek
% czasu ustalania z wartości 2.83s (tab 2) do odpowiednio 2.42s (tab 6)
% oraz 2.4s (tab 6) jest to spadek o 14.4% (tab 10) dla układu gdzie Q=10
% oraz 15.9% (tab 12) w przypadku układu regulacji gdzie Q=100. Wynika z
% tego że wraz ze wzrostem wartości Q układ zwiększa swoją dynamikę. 
%
% e) Relacja wag Q i R - Zaobserwowano, że w przypadkach, gdy wartość R jest 
% odwrotnością Q (np. R=0.1 przy Q=1 oraz R=1 przy Q=10), uzyskano identyczne 
% parametry czasowe i te same wartości spadków czasu ustalania. Wynika z tego, 
% że dla dynamiki układu kluczowy jest stosunek wag Q/R. Potwierdza to teoretyczną 
% właściwość algorytmu LQR, w którym macierz wzmocnień sprzężenia zwrotnego zależy 
% od relatywnego stosunku kar za uchyb stanu (Q) do kar za wysiłek sterowania (R), 
% a nie od ich bezwzględnych wartości.
%
% f) Wpływ dołączenia obserwatora ( Vs=3 ) - Porównanie układu idealnego 
% z układem z obserwatorem wykazało pogorszenie parametrów czasowych. 
% Odnotowano wzrost czasu narastania o 3.8% (tab 12) oraz czasu ustalania o 6.1%.(tab 12) 
% Jednocześnie zaobserwowano spadek przeregulowania o 4% (tab 12). 
% Wynika z tego, że dołączenie obserwatora o dynamice 3-krotnie większej 
% od układu wprowadza niewielkie opóźnienie, ale nie wpływa znacząco 
% na stabilność i dokładność uchyb pozostaje bliski zeru.
%
% g) Wpływ obserwatora (Vs=3) - Dołączenie obserwatora spowodowało 
% nieznaczny wzrost czasu narastania o 1.1% oraz czasu ustalania o 2.5%. 
% Przeregulowanie spadło o 1.3% (tab 13). Uchyb ustalony pozostał bliski 
% zeru. Wyniki są bardziej zbliżone do układu idealnego niż w przypadku 
% wolniejszego obserwatora (x3).  
%
% h) Parametry wyregulowanego układu ciągłego - W celu spełnienia założonych 
% kryteriów jakości regulacji (tr ~ -tr-0.5s, tus ~ tus-0.4s, uchyb ≈ 0) przyjęto 
% wartości R=0.1, Q=5 oraz Vs=5. W układzie bez obserwatora uzyskano spadek 
% czasu narastania z 2.73s (tab 2) do 2.21s (tab 15). W układzie z obserwatorem 
% odnotowano spadek czasu narastania z 2.73s (tab 1) do 2.21s (tab 15), co 
% potwierdza spełnienie pierwszego kryterium. 
% W przypadku czasu ustalania uzyskano spadek z 4.35s (tab 1) do 4s (tab 15) 
% dla układu bez obserwatora, natomiast w układzie z obserwatorem osiągnięto 
% redukcję z 4.61s (tab 1) do 4.17s (tab 15).
%
% i) Zależność czasu próbkowania od dynamiki obserwatora - Pierwotnie dobrany 
% czas próbkowania Tp = 0.32s (wynikający z pasma obiektu ok. 1 rad/s) okazał 
% się niewystarczający po zwiększeniu dynamiki obserwatora. Zwiększenie 
% szybkości obserwatora (3x oraz 5x) wymusiło proporcjonalne skrócenie czasu 
% próbkowania (odpowiednio Tp/3 oraz Tp/5). Wynika z tego, że w układach 
% z obserwatorem czas próbkowania musi być dostosowany do najszybszego 
% podukładu (najdalszych biegunów), aby uniknąć utraty stabilności i zapewnić 
% poprawną pracę estymatora w wersji dyskretnej.
% 
% j) Wpływ dyskretyzacji na parametry układu wyregulowanego - Porównanie 
% modelu ciągłego (tab 15) z modelem dyskretnym dla Tp=0.06s (tab 21) 
% wykazało istotne różnice w dynamice. W przypadku układu bez obserwatora 
% parametry pozostały niemal identyczne (tr=2.22s, tus=4.02s). 
% Jednak dla układu z obserwatorem odnotowano gwałtowne 
% przyspieszenie odpowiedzi: czas narastania spadł z 2.28s (tab 15) 
% do 1.74s (tab 21), co stanowi spadek o 23.6%. Czas ustalania uległ 
% skróceniu z 4.1781s (tab 15) do 3.3s (tab 21), czyli o 21%. 
% Wynika z tego, że dyskretyzacja przy małym Tp i wysokiej dynamice 
% obserwatora (Vs=5) powoduje przemieszczenie biegunów skutkujące 
% znacznie agresywniejszą pracą układu w domenie Z.
%
% k) Optymalizacja parametrów układu dyskretnego - Aby zniwelować efekt 
% nadmiernego przyspieszenia odpowiedzi w domenie Z, dokonano redukcji 
% dynamiki obserwatora do Vs=3 oraz zmniejszenia wagi Q do poziomu 0.2 
% (przy Tp=0.1s). Porównanie z modelem ciągłym (tab 15) wykazało wysoką 
% zbieżność wyników: dla układu z obserwatorem czas narastania wyniósł 
% 2.2s (tab 24) względem 2.27s (tab 15), a czas ustalania osiągnął 
% 4s (tab 24) względem 4.17s (tab 15). % Wynika z tego, że w 
% przypadku dyskretnej postaci regulatora należy stosować łagodniejsze 
% nastawy (mniejsze Q, mniejsze Vs) niż w modelu ciągłym. Wynika to z 
% faktu, że dyskretyzacja naturalnie zwiększa agresywność sterowania, 
% co przy zachowaniu ciągłych nastaw prowadziłoby  do przeregulowań i 
% skrócenia czasów poza założone kryteria jakości.
% 
% l) Wrażliwość na zmiany Tp przy wysokiej dynamice - Zaobserwowano, że dla 
% szybkiego obserwatora (Vs=5) układ wykazuje dużą wrażliwość numeryczną. 
% Zbyt małe Tp w połączeniu z dużymi wzmocnieniami obserwatora może prowadzić 
% do gwałtownej niestabilności. Większe wartości Tp, mimo wprowadzania 
% opóźnień, zapewniają szerszy margines stabilności numerycznej, filtrując 
% szybkozmienne błędy dyskretyzacji.
%
% m) Zależność stabilności od zapasu fazy - Zaobserwowano, że przy wysokiej 
% dynamice obserwatora układ staje się bardzo wrażliwy na opóźnienia 
% fazowe wprowadzane przez ZOH. Zmniejszanie Tp przy szybkim obserwatorze 
% gwałtownie redukuje zapas fazy, co przy minimalnym błędzie doboru 
% parametrów prowadzi do utraty stabilności. Większe Tp ogranicza pasmo 
% układu, co paradoksalnie pozwala na zachowanie większego marginesu 
% fazy kosztem ogólnej szybkości odpowiedzi.

%% === Rekomendacja ===
% 1. Projekt w dziedzinie ciągłej (S): Jako punkt wyjścia należy przyjąć 
%    nastawy R=0.1 oraz Q=5 przy Vs=5 (tab 15). Pozwala to na wyznaczenie 
%    docelowej, idealnej trajektorii układu (tr=2.2s, tus=4.17s) bez wpływu 
%    opóźnień próbkowania.
% 
% 2. Dobór czasu próbkowania: Przy przejściu na model cyfrowy należy skrócić 
%    Tp co najmniej do wartości Tp = 0.1s dla Vs=3 lub Tp = 0.06s dla Vs=5. 
%    Stosowanie Tp wynikającego tylko z pasma obiektu (0.32s) jest 
%    niedopuszczalne przy szybkim obserwatorze, gdyż prowadzi do 
%    natychmiastowej niestabilności (tab 19).
% 
% 3. Korekta parametrów w dziedzinie Z: Implementator musi uwzględnić, że 
%    układ dyskretny z obserwatorem wykazuje o ok. 23% szybszą dynamikę 
%    (tab 21). Aby zachować parametry czasowe z modelu ciągłego, należy 
%    dokonać redukcji dynamiki obserwatora do Vs=3 oraz drastycznego 
%    zmniejszenia wagi Q (z 5 na 0.2).
% 
% 4. Weryfikacja końcowa: Ostateczny układ dyskretny (Tp=0.1, Q=0.2, Vs=3) 
%    należy uznać za optymalny, ponieważ najdokładniej odwzorowuje założone 
%    kryteria jakościowe modelu ciągłego przy zachowaniu 
%    wysokiej stabilności numerycznej i braku przeregulowań (tab 24).
