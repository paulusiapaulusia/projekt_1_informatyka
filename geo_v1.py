from math import sin, cos, sqrt, atan, atan2, degrees, radians
import sys
o = object()
print(sys.arg)
import numpy as np
    

class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2


    
    def xyz2plh(self, X, Y, Z, output = 'dec_degree'):
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
        r   = sqrt(X**2 + Y**2)           # promień
        lat_prev = atan(Z / (r * (1 - self.ecc2)))    # pierwsze przybliilizenie
        lat = 0
        while abs(lat_prev - lat) > 0.000001/206265:    
            lat_prev = lat
            N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
            h = r / cos(lat_prev) - N
            lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lon = atan(Y/X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2);
        h = r / cos(lat) - N       
        if output == "dec_degree":
            return degrees(lat), degrees(lon), h 
        elif output == "dms":
            lat = self.deg2dms(degrees(lat))
            lon = self.deg2dms(degrees(lon))
            return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - output format not defined")
            

if __name__ == "__main__":
    # utworzenie obiektu
    geo = Transformacje(model = "wgs84")
    # dane XYZ geocentryczne
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    phi, lam, h = geo.xyz2plh(X, Y, Z)
    print(phi, lam, h)
    # phi, lam, h = geo.xyz2plh2(X, Y, Z)
    # print(phi, lam, h)

    def pliczek(self, plik, funkcja):
        data = np.genfromtxt(plik,  delimiter = " ")
        if funkcja == "XYZ_BLH":
            X = data[:,0]
            Y = data[:,1]
            Z = data[:,2]
                # to zmienic e starej wersji, bedzie latwiej
            blh = self.hirvonen(X, Y, Z)
            np.savetxt(f"WYNIK_{funkcja}.txt", blh, delimiter=";")

        elif funkcja == "BLH_XYZ":
            fi = np.deg2rad(data[:,0])
            lam = np.deg2rad(data[:,1])
            h = data[:,2]
            XYZ = self.filh2XYZ(fi, lam, h)
            np.savetxt(f"WYNIK_{funkcja}.txt",XYZ, delimiter=";")
            
        elif funkcja == "XYZ_NEU":
            X0 = data[0,0]
            Y0 = data[0,1]
            Z0 = data[0,2]
            X = data[1,0]
            Y = data[1,1]
            Z = data[1,2]
            neu = self.xyz2neup(X, Y, Z, X0, Y0, Z0)
            np.savetxt(f"WYNIK_{funkcja}.txt", neu, delimiter=";")
            
        elif funkcja == "BL_PL1992":
            fi = np.deg2rad(data[:,0])
            lam = np.deg2rad(data[:,1])
            wsp92 = self.cale92(fi, lam)
            np.savetxt(f"WYNIK_{funkcja}.txt", wsp92, delimiter=";")
            
        elif funkcja == "BL_PL2000":
            fi = np.deg2rad(data[:,0])
            lam = np.deg2rad(data[:,1])
            wsp00 = self.cale00(fi, lam)
            np.savetxt(f"WYNIK_{funkcja}.txt", wsp00, delimiter=";")
    elip = {'WGS84':[6378137.000, 0.00669438002290], 'GRS80':[6378137.000, 0.00669438002290], 'Elipsoida Krasowskiego':[6378245.000, 0.00669342162296]}
    funkcja = {'XYZ_BLH' : 'hirvonen', 'BLH_XYZ' : 'filh2XYZ', 'XYZ_NEU' : 'xyz2neup', 'BL_PL1992' : 'cale92', 'BL_PL2000' : 'cale00'}
        
    try:
        geo = Transformations(elip[args.elip.upper()])
        finito = geo.pliczek(args.plik, args.funkcja.upper())
        print("Zapisano")
    except KeyError():
        print(f"Podana funkcja/elipsoida nie istnieją, proszę upewnij się, że korzystasz z istniejących elipsoid")
    except AttributeError:
        print("Podana funkcja/elipsoida nie istnieje, proszę wprowadzić dostępne wartosci.")
    except FileNotFoundError:
        print("Nie znaleziono takiego pliku. Proszę spróbować wprowadzić inny plik.")
    except IndexError:
        print("Nieodpowiedni format pliku. Proszę wprowadzić dane do pliku tak jak pokazano w przyładzie.")
    except ValueError:
        print("Nieodpowiedni format pliku. Proszę wprowadzić dane do pliku tak jak pokazano w przyładzie.")

if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser(description="Podaj plik")
        parser.add_argument("-plik", type = str, help = "Podaj nazwę pliku, w którym znajdują się dane wejsciowe (ps. oprócz nazwy podaj rozszerzenie:)")
        parser.add_argument("-elip", type = str, help = "Wybierz elipsoidę, na której ma wykonać się transformacja, wpisz jedną: 'WGS84', 'GRS80', 'Elipsoida Krasowskiego' ")
        parser.add_argument("-funkcja", type = str, help = "Wybierz transformację jaką chcesz obliczyć: 'XYZ_BLH', 'BLH_XYZ', 'XYZ_NEU' ")
        args = parser.parse_args()
    except SyntaxError:
        print(f"Niestety nie ma takiego pliku. Spróbuj podać pełną scieżkę do pliku lub upewnij się że wpisujesz dobrą nazwę")

if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser(description="Podaj plik")
        parser.add_argument("-plik", type = str, help = "Podaj nazwę pliku, w którym znajdują się dane wejsciowe (ps. oprócz nazwy podaj rozszerzenie:)")
        parser.add_argument("-elip", type = str, help = "Wybierz elipsoidę, na której ma wykonać się transformacja, wpisz jedną: 'WGS84', 'GRS80', 'Elipsoida Krasowskiego' ")
        parser.add_argument("-funkcja", type = str, help = "Wybierz transformację jaką chcesz obliczyć: 'XYZ_BLH', 'BLH_XYZ', 'XYZ_NEU' ")
        args = parser.parse_args()
    except SyntaxError:
        print(f"Niestety nie ma takiego pliku. Spróbuj podać pełną scieżkę do pliku lub upewnij się że wpisujesz dobrą nazwę")


coords_plh= []        
with open('przyklad/wsp_inp.txt') as f:
    lines = f.readlines()
    coords_lines = lines[4:]
    for coords_lines in coords_lines:
        coords_lines= coords_lines('\n')
        x_str, y_str, z_str= coords_lines.split(',')
        x, y, z= (float( x_str),float( y_str), float(z_str))
        phi, lam, h = geo.xyzwplh(x,y,z)
        coords_plh.append([phi,lam,h])

input_file_path= sys.args[-1]        
with open(input_file_path, 'w') as f:
    f.write('phi[deg], lam[deg], h[m]\n')
    for coords_list in coords_plh:
        line= ','.join([str(coord) for coord in coords_list])
        f.writelines(line + '\n')
        
#kras-->grs80--->2000/92
#2 tygodnie---13 maj--ma pobierac dane tylko z pliku
if '--plh2zyx' in sys.args:      
    def plh2xyz(self, phi, lam, h):     
        Rm= self.a / sqrt(1 - self.acc) * sin(phi)**2
        x= (Rm +h)*cos(phi) *cos(lam)
        y= (Rm +h) *cos(phi) *sin(lam)
        z= (Rm +h) * sin(pho) - q
    
#lalal
        
if '--xyz2plh' in sys.args:        
    def hirvonen(self, X, Y, Z):
        '''
        Następujący algorytm przelicza współrzędne z układu ortokartezjańskiego na współrzędne geodezyjne.

        Parameters
        ----------
        X : float
        Y : float
        Z : float

        Returns: fi, lamda, h 
        -------
        None.

        '''
        flh = []
        for X,Y,Z in zip(X,Y,Z):
            p = np.sqrt(X**2 + Y**2)
            fi = np.arctan(Z / (p * (1 - self.e2)))
            while True:
                N = self.Npu(fi)
                h = p / np.cos(fi) - N
                fip = fi     #fip - fi poprzednie, fi - fi nowe
                fi = np.arctan(Z / (p * (1 - N * self.e2 / (N + h))))
                if abs(fip - fi) < (0.000001/206265):
                    break
        
            lam = np.arctan2(Y, X)
            flh.extend([np.rad2deg(fi), np.rad2deg(lam), h])
        return(flh)

        
if '--rneu' in sys.args:         # XYZ ---> NEU
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
   
 
if '--xyz2neup' in sys.args:
   def xyz2neup(self, X, Y, Z, X0, Y0, Z0):
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
       fi = np.arctan(Z0 / (p*(1 - self.e2)))
       while True:
           N = self.Npu(fi)
           h = (p / np.cos(fi)) - N
           fi_poprzednia = fi
           fi = np.arctan((Z0 / p)/(1-((N * self.e2)/(N + h))))
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

       # TRANSFORMACJA WSP BL ---> 1992
if '--bl292' in sys.args:
   def cale92(self, fi, lam):
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
           b2 = (self.a**2) * (1-self.e2)   #krotsza polowa
           e2p = (self.a**2 - b2 ) / b2   #drugi mimosrod elipsy
           dlam = lam - lam0
           t = np.tan(fi)
           ni = np.sqrt(e2p * (np.cos(fi))**2)
           N = self.Npu(fi)
           sigma = self.Sigma(fi)
           
           xgk = sigma + ((dlam**2)/2)*N*np.sin(fi)*np.cos(fi) * ( 1+ ((dlam**2)/12)*(np.cos(fi))**2 * ( 5 - (t**2)+9*(ni**2) + 4*(ni**4)     )  + ((dlam**4)/360)*(np.cos(fi)**4) * (61-58*(t**2)+(t**4) + 270*(ni**2) - 330*(ni**2)*(t**2))  )
           ygk = (dlam*N* np.cos(fi)) * (1+(((dlam)**2/6)*(np.cos(fi))**2) *(1-(t**2)+(ni**2))+((dlam**4)/120)*(np.cos(fi)**4)*(5-18*(t**2)+(t**4)+14*(ni**2)-58*(ni**2)*(t**2)) )
                       
           x92 = xgk*m - 5300000
           y92 = ygk*m + 500000
           wsp.append([x92, y92]) 
           
       return(wsp)
           
                     
       # TRANSFORMACJA WSP BL ---> 2000
if '--bl200' in sys.args:
   def cale00(self, fi, lam):
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
                    
           b2 = (self.a**2) * (1-self.e2)   #krotsza polos
           e2p = ( self.a**2 - b2 ) / b2   #drugi mimosrod elipsy
           dlam = lam - lam0
           t = np.tan(fi)
           ni = np.sqrt(e2p * (np.cos(fi))**2)
           N = self.Npu(fi)
           sigma = self.Sigma(fi)
       
           xgk = sigma + ((dlam**2)/2)*N*np.sin(fi)*np.cos(fi) * ( 1+ ((dlam**2)/12)*(np.cos(fi))**2 * ( 5 - (t**2)+9*(ni**2) + 4*(ni**4)     )  + ((dlam**4)/360)*(np.cos(fi)**4) * (61-58*(t**2)+(t**4) + 270*(ni**2) - 330*(ni**2)*(t**2))  )
           ygk = (dlam*N* np.cos(fi)) * (1+(((dlam)**2/6)*(np.cos(fi))**2) *(1-(t**2)+(ni**2))+((dlam**4)/120)*(np.cos(fi)**4)*(5-18*(t**2)+(t**4)+14*(ni**2)-58*(ni**2)*(t**2)) )
                    
           x00 = xgk * m
           y00 = ygk * m + strefa*1000000 + 500000
           wsp.append([x00, y00])
       return(wsp)
   
    

    
