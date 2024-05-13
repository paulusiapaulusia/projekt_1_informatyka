from math import sin, cos, sqrt, atan, atan2, degrees, radians
import sys
import numpy as np
from argparse import ArgumentParser #to trzeba dodac do pliku tekstowego
import argparse


class Transformacje:
    def __init__(self, elipsoida):
        self.a, self.ecc2 = elip[elipsoida]

    def Npu(self,fi):
        N = self.a / np.sqrt(1 - self.ecc2 * np.sin(fi)**2)
        
        return N
    
    def plh2xyz(self, fi, lam, h):     
       for fi, lam, h in zip(fi, lam, h):
            N = self.Npu(fi)
            x = (N+h)*np.cos(fi)*np.cos(lam)
            y = (N+h)*np.cos(fi)*np.sin(lam)
            z = (N*(1-self.ecc2)+h)*np.sin(fi)  
            return(x, y, z)

    def hirvonen(self, X, Y, Z, output='dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny.
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.
        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie orto-kartezjańskim,

        Returns
        -------
        lat
            [stopnie dziesiętne] - szerokość geodezyjna
        lon
            [stopnie dziesiętne] - długośc geodezyjna.
        h : TYPE
            [metry] - wysokość elipsoidalna
        output [STR] - optional, defoulf
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """
        r = sqrt(X ** 2 + Y ** 2)  # promień
        lat_prev = atan(Z / (r * (1 - self.ecc2)))  # pierwsze przybliilizenie
        lat = 0
        while abs(lat_prev - lat) > 0.000001 / 206265:
            lat_prev = lat
            N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev) ** 2)
            h = r / cos(lat_prev) - N
            lat = atan((Z / r) * (((1 - self.ecc2 * N / (N + h)) ** (-1))))
        lon = atan(Y / X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat)) ** 2)
        h = r / cos(lat) - N
        if output == "dec_degree":
            return degrees(lat), degrees(lon), h
        elif output == "dms":
            lat = self.deg2dms(degrees(lat))
            lon = self.deg2dms(degrees(lon))
            return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - output format not defined")

   # XYZ ---> NEU
    def Rneu(self, fi, lam): #fi, lam--->radians
       '''
       Obliczenie macierzy Rneu

       Parameters
       ----------
       fi : float
       lam : float

       Returns: rneu
       -------
       None.

       '''
       Rneu = np.array([[-np.sin(fi)*np.cos(lam), -np.sin(lam), np.cos(fi)*np.cos(lam)],
                        [-np.sin(fi)*np.sin(lam),  np.cos(lam), np.cos(fi)*np.sin(lam)],
                        [             np.cos(fi),            0,             np.sin(fi)]])
       return(Rneu)
 
    def xyz2neu(self, X, Y, Z, X0, Y0, Z0):
       '''
       Przeliczenie wsp XYZ na neu
       
       Parameters
       ----------
       X : TYPE
           DESCRIPTION.
       Y : TYPE
           DESCRIPTION.
       Z : TYPE
           DESCRIPTION.
       X0 : TYPE
           DESCRIPTION.
       Y0 : TYPE
           DESCRIPTION.
       Z0 : TYPE
           DESCRIPTION.

       Returns
       -------
       None.

       '''
       neu = []
       p = np.sqrt(X0**2 + Y0**2)
       fi = np.arctan(Z0 / (p*(1 - self.ecc2)))
       while True:
           N = self.Npu(fi)
           h = (p / np.cos(fi)) - N
           fi_poprzednia = fi
           fi = np.arctan((Z0 / p)/(1-((N * self.ecc2)/(N + h))))
           if abs(fi_poprzednia - fi) < (0.000001/206265):
               break 
       N = self.Npu(fi)
       h = p/np.cos(fi) - N
       lam = np.arctan(Y0 / X0)
       
       R_neu = self.Rneu(fi, lam)
       X_sr = [X - X0, Y - Y0, Z - Z0] 
       X_rneu = R_neu.T@X_sr
       neu.append(X_rneu.T)
           
       return(neu)

    def bl92(self, fi, lam):
       '''
     Algorytm przelicza współrzędne geodezyjne (BL) na współrzędne w układzie 1992 (XY)      

       Parameters
       ----------
       fi : TYPE
           DESCRIPTION.
       lam : TYPE
           DESCRIPTION.

       Returns
       -------
       None.

       '''
       lam0 = (19*np.pi)/180
       m = 0.9993
       wsp = []
       for fi,lam in zip(fi,lam):
           b2 = (self.a**2) * (1-self.ecc2)   #krotsza polowa
           e2p = (self.a**2 - b2 ) / b2   #drugi mimosrod elipsy
           dlam = lam - lam0
           t = np.tan(fi)
           ni = np.sqrt(e2p * (np.cos(fi))**2)
           N = self.Npu(fi)
           A0 = 1 - (self.ecc2/4)-(3*self.ecc2**2/64)-(5*self.ecc2**3/256)
           A2 = (3/8)*(self.ecc2+(self.ecc2**2/4)+(15*self.ecc2**3/128))
           A4 = (15/256)*(self.ecc2**2+((3*self.ecc2**3)/4))
           A6 = (35*self.ecc2**3)/3072
           sigma = self.a *(A0*fi-A2*np.sin(2*fi)+A4*np.sin(4*fi)-A6*np.sin(6*fi))
           
           xgk = sigma + ((dlam**2)/2)*N*np.sin(fi)*np.cos(fi) * ( 1+ ((dlam**2)/12)*(np.cos(fi))**2 * ( 5 - (t**2)+9*(ni**2) + 4*(ni**4)     )  + ((dlam**4)/360)*(np.cos(fi)**4) * (61-58*(t**2)+(t**4) + 270*(ni**2) - 330*(ni**2)*(t**2))  )
           ygk = (dlam*N* np.cos(fi)) * (1+(((dlam)**2/6)*(np.cos(fi))**2) *(1-(t**2)+(ni**2))+((dlam**4)/120)*(np.cos(fi)**4)*(5-18*(t**2)+(t**4)+14*(ni**2)-58*(ni**2)*(t**2)) )
                       
           x92 = xgk*m - 5300000
           y92 = ygk*m + 500000
           wsp.append([x92, y92]) 
           
       return(wsp)

    def bl200(self, fi, lam):
       '''
       Następujący algorytm umożliwia przeliczenie współrzędnych geodezyjnych (BLH) na współrzędne w układzie 2000 (XY)

       Parameters
       ----------
       fi : TYPE
           DESCRIPTION.
       lam : TYPE
           DESCRIPTION.

       Returns
       -------
       None.

       '''
       m=0.999923
       print(fi, lam)
       wsp = []
       for fi,lam in zip(fi,lam):
           lam0=0 
           strefa = 0
           if lam >np.deg2rad(13.5) and lam < np.deg2rad(16.5):
               strefa = 5
               lam0 = np.deg2rad(15)
           elif lam >np.deg2rad(16.5) and lam < np.deg2rad(19.5):
               strefa = 6
               lam0 = np.deg2rad(18)
           elif lam >np.deg2rad(19.5) and lam < np.deg2rad(22.5):
               strefa =7
               lam0 = np.deg2rad(21)
           elif lam >np.deg2rad(22.5) and lam < np.deg2rad(25.5):
               strefa = 8
               lam0 = np.deg2rad(24)
           else:
               print("Punkt poza strefami odwzorowawczymi układu PL-2000")        
                    
           b2 = (self.a**2) * (1-self.ecc2)   #krotsza polos
           e2p = ( self.a**2 - b2 ) / b2   #drugi mimosrod elipsy
           dlam = lam - lam0
           t = np.tan(fi)
           ni = np.sqrt(e2p * (np.cos(fi))**2)
           N = self.Npu(fi)
           A0 = 1 - (self.ecc2/4)-(3*self.ecc2**2/64)-(5*self.ecc2**3/256)
           A2 = (3/8)*(self.ecc2+(self.ecc2**2/4)+(15*self.ecc2**3/128))
           A4 = (15/256)*(self.ecc2**2+((3*self.ecc2**3)/4))
           A6 = (35*self.ecc2**3)/3072
           sigma = self.a *(A0*fi-A2*np.sin(2*fi)+A4*np.sin(4*fi)-A6*np.sin(6*fi))
       
           xgk = sigma + ((dlam**2)/2)*N*np.sin(fi)*np.cos(fi) * ( 1+ ((dlam**2)/12)*(np.cos(fi))**2 * ( 5 - (t**2)+9*(ni**2) + 4*(ni**4)     )  + ((dlam**4)/360)*(np.cos(fi)**4) * (61-58*(t**2)+(t**4) + 270*(ni**2) - 330*(ni**2)*(t**2))  )
           ygk = (dlam*N* np.cos(fi)) * (1+(((dlam)**2/6)*(np.cos(fi))**2) *(1-(t**2)+(ni**2))+((dlam**4)/120)*(np.cos(fi)**4)*(5-18*(t**2)+(t**4)+14*(ni**2)-58*(ni**2)*(t**2)) )
                    
           x00 = xgk * m
           y00 = ygk * m + strefa*1000000 + 500000
           wsp.append([x00, y00])
       return(wsp)



def wczytywanie_pliku(plik, funkcja, transformacje):
    poprawne_dane = []
    with open(plik, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                if any(char.isalpha() for char in line):
                    continue
                poprawne_dane.append(line.strip().split(','))
    # Konwersja poprawnych danych na tablicę numpy
    data = np.array(poprawne_dane, dtype=float)
    if funkcja == 'xyz2plh':
        X = data[0, 0]
        Y = data[0, 1]
        Z = data[0, 2]
        # to zmienic e starej wersji, bedzie latwiej
        blh = transformacje.hirvonen(X, Y, Z)
        np.savetxt(f"WYNIK_{funkcja}.txt", blh, delimiter=";")

    elif funkcja == 'plh2xyz':
        fi = np.deg2rad(data[:, 0])
        lam = np.deg2rad(data[:, 1])
        h = data[:, 2]
        XYZ = transformacje.plh2xyz(fi, lam, h)
        np.savetxt(f"WYNIK_{funkcja}.txt", XYZ, delimiter=";")

    elif funkcja == 'xyz2neu':
        X0 = data[0, 0]
        Y0 = data[0, 1]
        Z0 = data[0, 2]
        X = data[1, 0]
        Y = data[1, 1]
        Z = data[1, 2]
        neu = transformacje.xyz2neu(X, Y, Z, X0, Y0, Z0)
        np.savetxt(f"WYNIK_{funkcja}.txt", neu, delimiter=";")

    ## tutaj brakuje funkcji w klasie 
    elif funkcja == 'bl292':
        fi = np.deg2rad(data[:, 0])
        lam = np.deg2rad(data[:, 1])
        wsp92 = transformacje.bl292(fi, lam)
        np.savetxt(f"WYNIK_{funkcja}.txt", wsp92, delimiter=";")

    elif funkcja == 'bl200':
        fi = np.deg2rad(data[:, 0])
        lam = np.deg2rad(data[:, 1])
        wsp00 = transformacje.bl200(fi, lam)
        np.savetxt(f"WYNIK_{funkcja}.txt", wsp00, delimiter=";")


if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser(description="Podaj plik")
        parser.add_argument("-p", type = str, help = "Podaj nazwę pliku, w którym znajdują się dane wejsciowe (ps. oprócz nazwy podaj rozszerzenie:)")
        parser.add_argument("-el", type = str, help = "Wybierz elipsoidę, na której ma wykonać się transformacja, wpisz jedną: 'wgs84', 'grs80', 'KRASOWSKI' ")
        parser.add_argument("-t", type = str, help = "Wybierz transformację jaką chcesz obliczyć:  'xyz2plh', 'bl200', 'xyz2neu', 'xyz2plh', 'plh2xyz', 'rneu', 'bl292' ")
        args = parser.parse_args()
    except SyntaxError:
        print(f"Niestety nie ma takiego pliku. Spróbuj podać pełną scieżkę do pliku lub upewnij się że wpisujesz dobrą nazwę")


elip = {'WGS84':[6378137.000, 0.00669438002290], 'GRS80':[6378137.000, 0.00669438002290], 'KRASOWSKI':[6378245.000, 0.00669342162296]}
funkcja = { 'xyz2plh' :  'xyz2plh', 'plh2xyz' : 'plh2xyz', 'xyz2neu' : 'xyz2neu', 'bl292' : 'bl292', 'bl200' : 'bl200'}

try:
    if args.el == None:
        args.el = input(str('Podaj nazwe elipsoidy: '))
    if args.p == None:
        args.p = input(str('Wklej sciezke do pliku txt z danymi: '))
    if args.t == None:
        args.t = input(str('Jaka transformacje wykonac?: '))

    geo = Transformacje(args.el.upper())
    wczytywanie_pliku(args.p, args.t, geo)
    
    print('Plik wynikowy zostal utworzony.')

    wybor = input(str("Jezeli chcesz wykonac kolejna transformacje wpisz TAK jesli chcesz zakonczyc ENTER: ")).upper()
    args.el = None
    args.p = None
    args.t = None

except FileNotFoundError:
    print('Podany plik nie istnieje.')
except KeyError:
    print('Zle podana elipsoida lub transformacja.')
except IndexError:
    print('Zly format danych w pliku.')
except ValueError:
    print('Zly format danych w pliku.')
finally:
    print('Koniec programu')