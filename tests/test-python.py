from math import exp, log10, floor


def PbAlMarhounCorrelation(Yg, Yo, Rsb, T):        
    Pb = (5.38088 * 10 ** -3) * (Rsb ** 0.715082) * (Yg ** -1.87784) * (Yo ** 3.1437) * (T ** 1.32657)
        
    if Pb < 0:
        Pb = 0
        
    PbAlMarhoun = Pb
        
    return PbAlMarhoun

def PbAlShammasiCorrelation(Yg, Yo, Rsb, T):
    c1 = 5.527215
    c2 = -1.841408
    c3 = 0.783716
        
    Pb = (Yo ** c1) * exp(c2 * floor(Yo * Yg)) * (Rsb * (T + 460) * Yg) ** c3 
        
    PbAlShammasi = Pb
        
    return PbAlShammasi

def PbDeGhettoCorrelation(Yg, Rsb, API, T):
    if API <= 10: # Extra-heavy oil
        Pb = (Rsb / Yg) ** (1 / 1.1128) * 10.7025 / (10 ** ((0.0169 * API) - (0.00156 * T)))
    else:
        Pb = 15.7286*(((Rsb/Yg)**0.7885)*((10**(0.0020*T))/(10**(0.0142*API))))
        
    PbDeGhetto = Pb
        
    return PbDeGhetto



def PbDindorukChristmanCorrelation(Yg, Rsb, API, T):
    a1 = 1.42828 * 10 ** -10
    a2 = 2.844591797
    a3 = -6.74896 * 10 ** -4
    a4 = 1.225226436
    a5 = 0.033383304
    a6 = - 0.272945957
    a7 = -0.084226069
    a8 = 1.869979257
    a9 = 1.221486524
    a10 = 1.370508349
    a11 = 0.011688308
        
    A = ((a1 * T**a2) + (a3 * API**a4)) / (a5 + (2 * Rsb**a6 / Yg**a7))**2
       
    Pb = a8 * ((Rsb**a9 * 10**A / Yg**a10) + a11)
            
    PbDindorukChristman = Pb
        
    return PbDindorukChristman

def PbDoklaOsmanCorrelation(Yg, Rsb, API, T):
    Yo = 141.5 / ( API + 131.5 ) # Leandro
       
    Pb = (0.836386 * 10 ** 4) * (Yg ** -1.01049) * (Yo ** 0.107991) * (T ** -0.952584) * (Rsb ** 0.724047)
        
    PbDoklayOsman = Pb
        
    return PbDoklayOsman

    def PbGlasoCorrelation(Yg, Rsb, API, T, n2Concentration, co2Concentration, h2sConcentration):
        """ GLASO CORRELATION, CALCULATION OF BUBBLE POINT PRESSURE  
        
        DATA BANK:
        Based on 26 samples from the North Sea (collected from wells in the region 56 to 62°N) and 19 samples from the Middle East, Algeria, and several areas in the U.S.
        
        :see: Oistein Glaso. "Generalized Pressure-Volume-Temperature Correlations," Journal of Petroleum Technology , 1980.
        
        :precondition: T: 80 - 280 [°F]
        :precondition: API: 22.3 - 48.1 [°API]
        :precondition: Yg: 0.650 - 1.276 [ratio air=1]
        :precondition: Rs: 90 - 2637 [scf/STB]
        :precondition: n2Concentration: 0 - 26 [mol percent]
        :precondition: co2Concentration: 0 - 26 [mol percent]
        :precondition: h2sConcentration: 0 - 50 [mol percent]
        :precondition: Pb: 150 - 7127 [psig]
        
        :type Yg: number
        :param Yg: Gas specific gravity [ratio air=1]
        :type Rsb: number
        :param Rsb: Solution GOR at Pb [scf/STB]
        :type API: number
        :param API: API oil gravity [°API]
        :type T: number
        :param T: Temperature [°F]
        :type n2Concentration: number
        :param n2Concentration: molar fraction of nitrogen in the gas [gas mol/mixture mol]
        :type co2Concentration: number
        :param co2Concentration: molar fraction of carbon dioxide in the gas [gas mol/mixture mol]
        :type h2sConcentration: number
        :param h2sConcentration: molar fraction of hydrogen sulfide in the gas [gas mol/mixture mol]
        
        :return: Pb = Bubble point pressure [psia]
        """
        
        X = (Rsb / Yg) ** 0.816 * T ** 0.172 / API ** 0.989
        
        PbHC = 10 ** (1.7669 + (1.7447 * log10(X)) - (0.30218 * (log10(X)) ** 2))
        
        # Effects of nonhydrocarbons on bubble point pressure
        n2Effect = 1 + ((-2.65 * 10 **-4 * API + 5.5 * 10 **-3) * T + (0.0931 * API - 0.8295)) * n2Concentration + ((1.954 * 10 **-11 * API **4.699) * T + (0.027 * API - 2.366)) * n2Concentration **2
        co2Effect = 1 - 693.8 * co2Concentration * T **-1.553
        h2sEffect = 1 - (0.9035 + 0.0015 * API) * h2sConcentration + 0.019 * (45 - API) *  h2sConcentration **2
  
        Pb = PbHC * n2Effect * co2Effect * h2sEffect
        
        if Pb < 0:
            
            Pb = 0
        
        PbGlaso = Pb
        
        return PbGlaso

def PbHanafyCorrelation(Rsb):
    # Total flash gas-oil ratio
    Rsft = 69 + 1.071 * Rsb
        
    # Initial differential gas-oil ratio
    Rsi = 23.94 + 1.101 * Rsft
        
    Pb = (3.205 * Rsi) + 157.27
        
    PbHanafy = Pb
        
    return PbHanafy

def PbKartoatmodjoSchmidtCorrelation(Yg, Rsb, API, T, Tsep, Psep):
    C = 1 + 0.1595 * (API ** 0.4078) * (Tsep ** -0.2466) * log10(Psep / 114.7)
    Ygcorr = C * Yg
    if API > 30:
        C1 = 0.0315
        C2 = 0.7587
        C3 = 11.289
        C4 = 0.9143
    else:
        C1 = 0.05958
        C2 = 0.7972
        C3 = 13.1405
        C4 = 0.9986
    Pb = (Rsb / (C1 * (Ygcorr ** C2) * (10 ** (C3 * API / (T + 460) )))) ** C4
    if Pb < 0:
        Pb = 0
    PbKartoatmodjoSchmidt = Pb
    return PbKartoatmodjoSchmidt

def PbLasaterCorrelation(Yg, Rsb, API, T):
    Yo = 141.5 / (API + 131.5)
    # Effective molecular weight
    if API <= 40:
        Mo = 630 - (10 * API)
    else:
        Mo = 73110 * (API ** -1.562) 
            
    # Mol fraction of gas
    Ygfactor = (Rsb / 379.3) / ((Rsb / 379.3) + (350 * Yo / Mo))
        
    # Bubble point pressure factor
    if Ygfactor <= 0.6:
        Pbfactor = (0.679 * exp (2.786 * Ygfactor)) - 0.323 
    else:
        Pbfactor = (8.26 * (Ygfactor ** 3.56)) + 1.95
            
    Pb = ((Pbfactor) * (T + 459.6)) / Yg
        
    PbLasater = Pb
        
    return PbLasater

def PbPetroskyFarshadCorrelation(Yg, Rsb, API, T):
    X = 4.561 * 10 ** -5 * T ** 1.3911 - (7.916 * 10 ** -4 * API ** 1.541)
        
    Pb = 112.727 * (((Rsb ** 0.5774 / Yg ** 0.8439) * (10 ** X)) - 12.34)
        
    PbPetroskyFarshad = Pb
        
    return PbPetroskyFarshad

def PbStandingCorrelation(Yg, Rsb, API, T):
    Pb = 18.2 * (((Rsb / Yg) ** 0.83 * (10 ** ((0.00091 * T) - (0.0125 * API)))) - 1.4)
    return Pb

def PbTotalCFPCorrelation(Yg, Rsb, API, T):
    if API <= 10:
        c1 = 12.847
        c2 = 0.9636
        c3 = 0.000993
        c4 = 0.03417
    
    if 10 < API <= 35:
        c1 = 25.2755
        c2 = 0.7617
        c3 = 0.000835
        c4 = 0.011292
        
    if 35 < API <= 45:
        c1 = 216.4711
        c2 = 0.6922
        c3 = -0.000427
        c4 = 0.02314
            
    Pb = c1 * ((Rsb / Yg) ** c2) * (10 ** ((c3 * T) - (c4 * API)))
        
    PbTotalCFP = Pb
        
    return PbTotalCFP


def PbVasquezBeggsCorrelation(Yg, Rsb, API, T, Tsep, Psep):
    if API <= 30:
        c1 = 0.0362
        c2 = 1.0937
        c3 = 25.724
    else:
        c1 = 0.0178
        c2 = 1.187
        c3 = 23.931
        
    # Gas gravity that would result from separator conditions of 100 psig (approximately 114.7 psia)
    Ygs = Yg*(1.+5.912*10**-5*(API)*(Tsep)*log10(Psep/114.7))
    
    Pb = (Rsb / (c1 * Ygs * exp((c3 * API) / (T + 460)))) ** (1 / c2)
            
    PbVasquezBegg = Pb
    
    return PbVasquezBegg

def PbVelardeCorrelation(Yg, Rsb, API, T):
    x = (0.013098 * T ** 0.282372) - (8.2 * 10 ** -6 * API ** 2.176124)
    Pb = 1091.47 * ((Rsb ** 0.081465 * Yg ** -0.161488 * 10 ** x) - 0.740152) ** 5.354891
    PbVelarde = Pb
    return PbVelarde

def PbCegarraCorrelation(Yg, Rsb, API, T):
    if API < 29.9:
        l1 = 154.158
        l2 = 0.4577
        l3 = 0.0006680
        l4 = 0.000514
        l5 = 4.70257
    elif API >= 29.9:
        l1 = 809.238
        l2 = 0.32
        l3 = 0.00061
        l4 = 0.011 
        l5 = 1.1142
            
    Pb = l1 * ((((Rsb / Yg) ** l2) * (10 ** ((l3 * T) - (l4 * API)))) - l5)
        
    PbCegarra = Pb
        
    return PbCegarra

def PbPerezMLCorrelation(Yg, Rsb, API, T):
    X = (0.0002573 * T) - (0.0253643 * API)
    Pb = 10.667657 * ((Rsb / Yg) ** 1.00139101) * (10 ** X)
    PbPerezML = Pb
    return PbPerezML

def PbMillanArciaCorrelation(Yg, Rsb, API, T):
    X = (0.00091 * T) - (0.0125 * API)
    Pb = 25.3302 * ((((Rsb / Yg) ** 0.8303) * (10 ** X)) ** 0.9433)
    PbMillanArcia = Pb
    return PbMillanArcia


