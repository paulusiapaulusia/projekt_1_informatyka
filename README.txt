Projekt 1: TRANSFORMACJE

Program służy do transformacji współrzędnych.


DOSTĘPNE TRANSFORMACJE:
XYZ (geocentryczne) --> BLH (elipsoidalne)
BLH --> XYZ
XYZ --> NEUp (topocentryczne)
BL --> PL2000
BL --> PL1992

OBSŁUGIWANE ELIPSOIDY:
- GRS80
- WGS84
- Krasowski

WYMAGANIA:
Python 3.11 lub 3.12
biblioteka numpy
biblioteka math
biblioteka argparse

Program został napisany dla systemu operacyjnego Windows 11

OPIS PROGRAMU:
Program oczekuje trzech rodzajów informacji podanych za pomocą flag:
-p: Tutaj należy podać nazwę pliku, który zawiera dane potrzebne do przeprowadzenia transformacji. Ważne, aby plik miał odpowiednie rozszerzenie (txt).
-el: W tym polu należy wpisać nazwę wybranej elipsoidy, na której chcemy przeprowadzić transformację. Dostępne opcje to: 'WGS84', 'GRS80' oraz 'Elipsoida Krasowskiego'.
-t: W tej sekcji należy wpisać nazwę transformacji, którą chcemy zastosować. Dostępne możliwości to: 'XYZ2BLH', 'BLH2XYZ', 'XYZ2NEU', 'BL92' oraz 'BL200'.

Po wyborze parametrów i załadowaniu pliku z danymi utworzy się plik tekstowy zawierający wyniki wykonanych obliczeń.


PRZYKŁADOWE WYWOŁANIE:

geo_v1.py -p wsp_inp.txt -el grs80 -t xyz2blh

--->W folderze w którym znajduje się program powinien się utworzyć plik wynikowy.Wyświetli się następujący komunikat:

Plik wynikowy zostal utworzony.

--->Program dopuszcza zrobienie kolejnej transformacji bez ponownego uruchamiania. Odrazu po zapisaniu pliku tekstowego dostaniemy komunikat:

Jeżeli chcesz wykonać kolejną transformację, wpisz TAK. Wciśnij ENTER aby zakończyć: 

--->Wpisanie "TAK" pozwoli nam na ponowne wywołanie programu. Natomiast jeśli użytkownik zdecyduje się zakończyć działanie i wciśnie "ENTER", pojawi się poniższy komunikat:

Koniec transformacji

--->Całość wygląda w następujący sposób:


Podaj nazwę elipsoidy: GRS80
Wklej scieżkę do pliku z danymi (format musi być w formie txt): wsp_inp.txt
Jaką transformację program ma wykonać?: hirvonen
Plik wynikowy zostal utworzony.
Jeżeli chcesz wykonać kolejną transformację, wpisz TAK. Wciśnij ENTER aby zakończyć: TAK
Koniec transformacji
Podaj nazwę elipsoidy: WGS84
Wklej scieżkę do pliku z danymi (format musi być w formie txt): wsp_inp.txt
Jaką transformację program ma wykonać?: bl2000
Plik wynikowy zostal utworzony.
Jeżeli chcesz wykonać kolejną transformację, wpisz TAK. Wciśnij ENTER aby zakończyć: 
Koniec transformacji

--->Jeśli użytkownik przy wywoływaniu programu flagą nie poda którejś wartości, program samoistnie zapyta o daną wartość. Przykład:
geo_v1.py -el grs80 -t xyz2blh
Wklej scieżkę do pliku z danymi (format musi być w formie txt):


OPIS DANYCH DO TRANSFORMACJI:

--------BL --> PL2000/PL1992
51.3456 21.2340
50.2309 20.3456
50.0976 19.3467

Pierwsza kolumna to szerokość geodezyjna[stopnie], a druga to długość geodezyjna[stopnie]. Plik wynikowy prezentuje się w następujący sposób:

WYNIK_bl200.txt
5690122.096 7516302.661
5566306.665 7453311.168
5552144.313 6596348.084

Pierwsza wartość to współrzędna X[m], a druga to współrzędna Y[m].

-----------BLH --> XYZ
plik_dane_BLH2XYZ.txt
51.3489 21.0122 100.000
50.2983 20.3490 110.000
49.3874 19.6266 111.000

Pierwsza wartość to szerokość geodezyjna[stopnie], druga to długość geodezyjna[stopnie], a trzecia to wysokość punktu[m].

WYNIK_plh2xyz.txt
3726411.578 1431345.843 4957958.223
3827682.529 1419625.341 4884135.412
3918212.663 1397262.818 4818798.584

Pierwsza wartość to współrzędna X[m], druga to współrzędna Y[m], a trzecia to współrzędna Z[m].

---------XYZ (geocentryczne) --> BLH (elipsoidalne)
3654515.7218380477 1403730.0433548905 5018560.040857761
3856423.941196924 1399432.680121372 4867579.706638583
3528330.6974752024 1190061.3100085442 5160990.220680703
3836640.8048792393 1175798.397398429 4941181.911009777
3730268.373914433 1317661.09184508 4986406.785425414

Pierwsza wartość to współrzędna X[m], druga to współrzędna Y[m], a trzecia to współrzędna Z[m].

WYNIK_xyz2plh.txt
52.2297000000 21.0122000000 78.000
50.0647000000 19.9450000000 219.000
54.3722000000 18.6386000000 5.000
51.1079000000 17.0385000000 118.000
51.7592000000 19.4550000000 182.000

Pierwsza wartość to szerokość geodezyjna[stopnie], druga to długość geodezyjna[stopnie], a trzecia to wysokość punktu[m].

XYZ --> NEUp
3673602.0707659787,1410163.7139986877,5002803.345368741
22348783.04399999976,14703137.28000000119,38477.139999999942
-23911309.4600000089,-10891803.3950000142,-2229962.48600000034
12567758.58699999936,12885129.25399999879,19449249.87600000203
19821597.2880000250,3882086.66999999925,17313733.10199999809
-5921297.01099999940,-24454287.61199999973,-8484889.633000001311

W pierwszym wierszu znajdują się kolejno współrzędne X0, Y0, Z0 odbiornika na Ziemi. W pozostałych wierszach pierwsza wartość to współrzędna X[m], druga to współrzędna Y[m], a trzecia to współrzędna Z[m].

WYNIK_xyz2neu.txt
0.000 0.000 0.000
-20549047.679 5717473.618 9754855.210
19314544.051 -1599327.517 -24268687.235
-889561.130 7525422.589 19027805.115
-4998318.921 -3479185.048 19527908.835
6058904.131 -20708041.190 -21849862.779

W wierszach znajdują się kolejno współrzędne X[m], Y[m], Z[m] w układzie topocentrycznym. W pierwszym wierszu są zera, które należy zignorować.


UWAGI:
->Program zwraca błędne wartości współrzędnych geodezyjnych PL2000 i PL1992 w odniesieniu do elipsoidy Krasowskiego
->Pliki muszą być w formacie .txt
->Pliki nadpisują się. Jeśli wykonamy transformację współrzędnych geodezyjnych na PL2000 w elipsoidzie WGS84, musimy zmienić nazwę utworzonego pliku by program nie nadpisał nam w przypadku wykonania tej samej transformacji na innej elipsoidzie.
->Dane w plikach muszą być oddzielone PRZECINKIEM, nie spacją.