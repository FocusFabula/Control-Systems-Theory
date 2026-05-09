# 📉 Control Systems Theory & Digital Signal Processing

Repozytorium zawiera projekty z zakresu klasycznej i nowoczesnej teorii sterowania, przetwarzania sygnałów oraz optymalizacji układów regulacji. Modele i obliczenia zostały zrealizowane w środowiskach **MATLAB** oraz **Simulink**.

---

## 📑 Zawartość

### 1. Synteza cyfrowych filtrów dolnoprzepustowych (IIR)
Projekt i analiza numeryczna filtrów o nieskończonej odpowiedzi impulsowej (NOI) wyznaczonych metodą biliniowej.

* **Charakterystyka:** Filtr dolnoprzepustowy o niskiej częstotliwości odcięcia ($f_c = 4,8\text{ Hz}$), zaprojektowany do eliminacji zakłóceń niskoczęstotliwościowych i wygładzania sygnałów pomiarowych.
* **Analiza:** Porównanie implementacji II oraz VIII rzędu pod kątem stabilności i nieliniowości fazowej.
* **Cel:** Wyznaczenie równań różnicowych do implementacji w systemach czasu rzeczywistego.

### 2. Optymalizacja nastaw regulatorów (Kryteria Kesslera)
Analiza porównawcza klasycznych metod doboru nastaw regulatorów PI/PID dla obiektów inercyjnych.

* **Zastosowane kryteria:** Optimum Modułowe (OM) oraz Optimum Symetryczne (OS).
* **Metodologia:** Porównanie wyników wyliczonych analitycznie (**metoda tablicowa**) z rzeczywistym zachowaniem modelu w pętli sprzężenia zwrotnego (**układ zamknięty**).
* **Weryfikacja:** Analiza porównawcza obu kryteriów pod kątem czasu regulacji, przeregulowania oraz odporności na zakłócenia w środowisku Simulink.


### 3. Sterowanie optymalne i estymacja stanu (LQR + Obserwator)
Zaawansowane sterowanie w przestrzeni stanów z uwzględnieniem ograniczeń pomiarowych obiektu.

* **Regulator LQR (Linear Quadratic Regulator):** Projektowanie macierzy wagowych $Q$ i $R$ w celu minimalizacji funkcji kosztu sterowania.
* **Obserwator Stanu (Luenberger Observer):** Konstrukcja układu estymacji zmiennych stanu w przypadku braku bezpośredniego dostępu pomiarowego do wszystkich składowych wektora stanu.
* **Analiza porównawcza:** * Porównanie dynamiki układu z idealnym sprzężeniem od stanu oraz układu z obserwatorem.
    * Badanie wpływu szybkości obserwatora na jakość regulacji i odporność na zakłócenia.



---

## 🛠️ Wykorzystane narzędzia
* **MATLAB / Filter Design Toolbox:** Synteza filtrów i obliczenia macierzowe (LQR, bieguny obserwatora).
* **Simulink:** Porównawcze symulacje czasowe przebiegów z obserwatorem i bez.

---
*Projekty stanowią część opracowań realizowanych na Wydziale Mechatroniki, Uzbrojenia i Lotnictwa WAT.*
