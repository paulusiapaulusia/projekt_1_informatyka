Projekt 1: TRANSFORMACJE

Program służy do transformacji współrzędnych.

Dostępne transformacje:

XYZ (geocentryczne) --> BLH (elipsoidalne)
BLH --> XYZ
XYZ --> NEUp (topocentryczne)
BL --> PL2000
BL --> PL1992

Obsługiwane elipsoidy:
- GRS80
- WGS84
- Krasowski

Wymagania:
Python 3.11 lub 3.12
biblioteka numpy
biblioteka math
biblioteka argparse

Program został napisany dla systemu operacyjnego Windows 11

Opis programu:
Program oczekuje trzech rodzajów informacji podanych za pomocą flag:
-plik: Tutaj należy podać nazwę pliku, który zawiera dane potrzebne do przeprowadzenia transformacji. Ważne, aby plik miał odpowiednie rozszerzenie.
-elip: W tym polu należy wpisać nazwę wybranej elipsoidy, na której chcemy przeprowadzić transformację. Dostępne opcje to: 'WGS84', 'GRS80' oraz 'Elipsoida Krasowskiego'.
-funkcja: W tej sekcji należy wpisać nazwę transformacji, którą chcemy zastosować. Dostępne możliwości to: 'XYZ2BLH', 'BLH2XYZ', 'XYZ2NEU', 'BL92' oraz 'BL200'.

Po wyborze parametrów i załadowaniu pliku z danymi utworzy się plik tekstowy zawierający wyniki wykonanych obliczeń...
