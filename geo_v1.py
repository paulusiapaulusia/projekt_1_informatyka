from math import sin, cos, sqrt, atan
import numpy as np
import argparse


class Transformacje:
    def __init__(self, elipsoida):
        self.a, self.ecc2 = elip[elipsoida]

    def Npu(self, fi):
        N = self.a / np.sqrt(1 - self.ecc2 * np.sin(fi)**2)
        return N
    
    def hirvonen(self, X, Y, Z):
        wyniki = []
        for X, Y, Z in zip(X, Y, Z):
            p = np.sqrt(X**2 + Y**2)
            fi = np.arctan(Z / (p * (1 - self.ecc2)))
            while True:
                N = self.Npu(fi)
                h = (p / np.cos(fi)) - N
                fi_poprzednia = fi
                fi = np.arctan((Z / p) / (1 - ((N * self.ecc2) / (N + h))))
                if abs(fi_poprzednia - fi) < (0.000001 / 206265):
                    break
            N = self.Npu(fi)
            h = p / np.cos(fi) - N
            lam = np.arctan2(Y, X)  # Zmiana z arctan(Y/X) na arctan2(Y, X)
            wyniki.append([np.rad2deg(fi), np.rad2deg(lam), h])
        return wyniki

    def odwrocony_hirvonen(self, fi_list, lam_list, h_list):
        wyniki = []
        for fi, lam, h in zip(fi_list, lam_list, h_list):
            N = self.Npu(fi)
            Xk = (N + h) * np.cos(fi) * np.cos(lam)
            Yk = (N + h) * np.cos(fi) * np.sin(lam)
            Zk = (N * (1 - self.ecc2) + h) * np.sin(fi)
            wyniki.append([Xk, Yk, Zk])
        return wyniki

    def Rneu(self, fi, lam): 
        Rneu = np.array([[-np.sin(fi) * np.cos(lam), -np.sin(lam), np.cos(fi) * np.cos(lam)],
                         [-np.sin(fi) * np.sin(lam),  np.cos(lam), np.cos(fi) * np.sin(lam)],
                         [np.cos(fi), 0, np.sin(fi)]])
        return Rneu
    
    def xyz2neu(self, X, Y, Z, X0, Y0, Z0):
        neu = []
        p = np.sqrt(X0**2 + Y0**2)
        fi = np.arctan(Z0 / (p * (1 - self.ecc2)))
        while True:
            N = self.Npu(fi)
            h = (p / np.cos(fi)) - N
            fi_poprzednia = fi
            fi = np.arctan((Z0 / p) / (1 - ((N * self.ecc2) / (N + h))))
            if abs(fi_poprzednia - fi) < (0.000001 / 206265):
                break 
        N = self.Npu(fi)
        h = p / np.cos(fi) - N
        lam = np.arctan2(Y0, X0)  # Zmiana z arctan(Y0/X0) na arctan2(Y0, X0)

        R_neu = self.Rneu(fi, lam)
        for x, y, z in zip(X, Y, Z):
            X_sr = [x - X0, y - Y0, z - Z0]
            X_rneu = R_neu.T @ X_sr
            neu.append(X_rneu.T)

        return neu

    def bl92(self, fi, lam):
        lam0 = (19 * np.pi) / 180
        m = 0.9993
        wsp = []
        for fi, lam in zip(fi, lam):
            b2 = (self.a ** 2) * (1 - self.ecc2)
            e2p = (self.a ** 2 - b2) / b2
            dlam = lam - lam0
            t = np.tan(fi)
            ni = np.sqrt(e2p * (np.cos(fi)) ** 2)
            N = self.Npu(fi)
            A0 = 1 - (self.ecc2 / 4) - (3 * self.ecc2 ** 2 / 64) - (5 * self.ecc2 ** 3 / 256)
            A2 = (3 / 8) * (self.ecc2 + (self.ecc2 ** 2 / 4) + (15 * self.ecc2 ** 3 / 128))
            A4 = (15 / 256) * (self.ecc2 ** 2 + ((3 * self.ecc2 ** 3) / 4))
            A6 = (35 * self.ecc2 ** 3) / 3072
            sigma = self.a * (A0 * fi - A2 * np.sin(2 * fi) + A4 * np.sin(4 * fi) - A6 * np.sin(6 * fi))

            xgk = sigma + ((dlam ** 2) / 2) * N * np.sin(fi) * np.cos(fi) * (1 + ((dlam ** 2) / 12) * (
                    np.cos(fi)) ** 2 * (5 - (t ** 2) + 9 * (ni ** 2) + 4 * (ni ** 4)) + ((dlam ** 4) / 360) * (
                                             np.cos(fi) ** 4) * (61 - 58 * (t ** 2) + (t ** 4) + 270 * (ni ** 2) - 330 * (ni ** 2) * (
                            t ** 2)))
            ygk = (dlam * N * np.cos(fi)) * (1 + (((dlam) ** 2 / 6) * (np.cos(fi) ** 2) * (1 - (t ** 2) + ni ** 2)) + (
                    ((dlam ** 4 / 120) * (np.cos(fi) ** 4)) * (
                            5 - 18 * (t ** 2) + (t ** 4) + (14 * ni ** 2) - (58 * ni ** 2) * (t ** 2))))

            x92 = xgk * m - 5300000
            y92 = ygk * m + 500000
            wsp.append([x92, y92])

        return wsp
    def bl200(self, fi, lama, m=0.999923):
        wyniki = []
        for fi, lama in zip(fi, lama):
            lama0 = 0
            strefa = 0
            if lama > np.deg2rad(13.5) and lama < np.deg2rad(16.5):
                strefa = 5
                lama0 = np.deg2rad(15)
            elif lama > np.deg2rad(16.5) and lama < np.deg2rad(19.5):
                strefa = 6
                lama0 = np.deg2rad(18)
            elif lama > np.deg2rad(19.5) and lama < np.deg2rad(22.5):
                strefa = 7
                lama0 = np.deg2rad(21)
            elif lama > np.deg2rad(22.5) and lama < np.deg2rad(25.5):
                strefa = 8
                lama0 = np.deg2rad(24)
            b2 = self.a ** 2 * (1 - self.ecc2)
            ep2 = (self.a ** 2 - b2) / b2
            dellama = lama - lama0
            t = np.tan(fi)
            ni2 = ep2 * (np.cos(fi) ** 2)
            N = self.Npu(fi)
            A0 = 1 - (self.ecc2 / 4) - (3 * self.ecc2 ** 2 / 64) - (5 * self.ecc2 ** 3 / 256)
            A2 = (3 / 8) * (self.ecc2 + (self.ecc2 ** 2 / 4) + (15 * self.ecc2 ** 3 / 128))
            A4 = (15 / 256) * (self.ecc2 ** 2 + ((3 * self.ecc2 ** 3) / 4))
            A6 = (35 * self.ecc2 ** 3) / 3072
            sigma = self.a * (A0 * fi - A2 * np.sin(2 * fi) + A4 * np.sin(4 * fi) - A6 * np.sin(6 * fi))
            xgk = sigma + (((dellama ** 2 / 2) * N * np.sin(fi) * np.cos(fi)) * (
                        1 + ((dellama ** 2 / 12) * (np.cos(fi) ** 2) * (5 - t ** 2 + 9 * ni2 + 4 * ni2 ** 2)) + (
                            (dellama ** 4 / 360) * (np.cos(fi) ** 4) * (
                            61 - 58 * t ** 2 + t ** 4 + 270 * ni2 - 330 * ni2 * t ** 2))))
            ygk = (dellama * N * np.cos(fi)) * (
                        1 + (((dellama ** 2 / 6) * (np.cos(fi) ** 2) * (1 - t ** 2 + ni2)) + (
                            (((dellama ** 4 / 120) * (np.cos(fi) ** 4)) * (
                            5 - (18 * t ** 2) + t ** 4 + (14 * ni2) - (58 * ni2 * t ** 2))))))
            x2000 = xgk * m
            y2000 = ygk * m + (strefa * 1000000) + 500000
            wyniki.append([x2000, y2000])

        return wyniki

def wczytywanie_pliku(plik, funkcja, transformacje):
    poprawne_dane = []
    with open(plik, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                if any(char.isalpha() for char in line):
                    continue
                poprawne_dane.append(line.strip().split(','))
    data = np.array(poprawne_dane, dtype=float)
    if funkcja == 'xyz2plh':
        X = data[:, 0]
        Y = data[:, 1]
        Z = data[:, 2]
        blh = transformacje.hirvonen(X, Y, Z)
        np.savetxt(f"WYNIK_{funkcja}.txt", blh, delimiter=";", fmt='%0.10f %0.10f %0.3f')

    elif funkcja == 'plh2xyz':
        fi = np.deg2rad(data[:, 0])
        lam = np.deg2rad(data[:, 1])  
        h = data[:, 2]
        XYZ = transformacje.odwrocony_hirvonen(fi, lam, h)
        np.savetxt(f"WYNIK_{funkcja}.txt", XYZ, delimiter=";", fmt='%0.3f %0.3f %0.3f')

    elif funkcja == 'xyz2neu':
        X0 = data[0, 0]
        Y0 = data[0, 1]
        Z0 = data[0, 2]
        X = data[:, 0]  # Wczytaj wszystkie wartości X z pliku
        Y = data[:, 1]  # Wczytaj wszystkie wartości Y z pliku
        Z = data[:, 2]  # Wczytaj wszystkie wartości Z z pliku
        neu = transformacje.xyz2neu(X, Y, Z, X0, Y0, Z0)
        np.savetxt(f"WYNIK_{funkcja}.txt", neu, delimiter=";", fmt='%0.3f %0.3f %0.3f')

    elif funkcja == 'bl292':
        fi = np.deg2rad(data[:, 0])
        lam = np.deg2rad(data[:, 1])
        wsp92 = transformacje.bl92(fi, lam)
        np.savetxt(f"WYNIK_{funkcja}.txt", wsp92, delimiter=";", fmt='%0.3f %0.3f')

    elif funkcja == 'bl200':
        fi = np.deg2rad(data[:, 0])
        lam = np.deg2rad(data[:, 1])
        wsp00 = transformacje.bl200(fi, lam)
        np.savetxt(f"WYNIK_{funkcja}.txt", wsp00, delimiter=";", fmt='%0.3f %0.3f')

if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser(description="Podaj plik")
        parser.add_argument("-p", type=str, help="Podaj nazwę pliku, w którym znajdują się dane wejsciowe (ps. oprócz nazwy podaj rozszerzenie:)")
        parser.add_argument("-el", type=str, help="Wybierz elipsoidę, na której ma wykonać się transformacja, wpisz jedną: 'wgs84', 'grs80', 'KRASOWSKI'")
        parser.add_argument("-t", type=str, help="Wybierz transformację jaką chcesz obliczyć:  'xyz2plh', 'bl200', 'xyz2neu', 'xyz2plh', 'plh2xyz', 'rneu', 'bl292'")
        args = parser.parse_args()
    except SyntaxError:
        print(f"Niestety nie ma takiego pliku. Spróbuj podać pełną scieżkę do pliku lub upewnij się że wpisujesz dobrą nazwę")

elip = {'WGS84':[6378137.000, 0.00669438002290], 'GRS80':[6378137.000, 0.00669438002290], 'KRASOWSKI':[6378245.000, 0.00669342162296]}
funkcja = {'xyz2plh': 'xyz2plh', 'plh2xyz': 'plh2xyz', 'xyz2neu': 'xyz2neu', 'bl292': 'bl292', 'bl200': 'bl200'}

while True:
    try:
        if args.el is None:
            args.el = input(str('Podaj nazwę elipsoidy: '))
        if args.p is None:
            args.p = input(str('Wklej scieżkę do pliku z danymi (format musi być w formie txt): '))
        if args.t is None:
            args.t = input(str('Jaką transformację program ma wykonać?: '))

        geo = Transformacje(args.el.upper())
        wczytywanie_pliku(args.p, args.t, geo)
        
        print('Plik wynikowy został utworzony.')

        wybor = input(str("Jeżeli chcesz wykonać kolejną transformację, wpisz TAK. Wciśnij ENTER aby zakończyć: ")).upper()
        args.el = None
        args.p = None
        args.t = None

        if wybor != 'TAK':
            break

    except FileNotFoundError:
        print('Nie można znaleźć danego pliku.')
    except KeyError:
        print('Źle podana elipsoida lub transformacja.')
    except IndexError:
        print('Zły format danych w pliku.')
    except ValueError:
        print('Zły format danych w pliku.')
    finally:
        print('Koniec transformacji')