from math import exp, log10, floor, log

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

def PbManucciRosalesCorrelation(Yg, Rsb, API, T):
    X = (0.000922 * T) - (0.0072 * API)
    Pb = 84.88 * (((Rsb / Yg) ** 0.53) * (10 ** X))
    PbManucciRosales = Pb
    return PbManucciRosales

def RsAlShammasiCorrelation(Yg, Pb, P, Yo, T, Rsb):
    # Transformation from °API to oil specific gravity
    #Yo = 141.5 / ( API + 131.5 ) # conversion Leandro
    print Yg
    print Pb
    print P
    print Yo
    print T
    print Rsb
    c1 = 5.527215
    c2 = -1.841408
    c3 = 0.783716
    if P >= Pb: # Logical condition
        Rs = Rsb
    else:
        Rs = ((P / (Yo ** c1 * exp(c2 * floor(Yo * Yg)))) ** (1 / c3)) / ((T + 460) * Yg)
    RsAlShammasi = Rs
    return RsAlShammasi


def RsAlMarhounCorrelation(Yg, Pb, P, Yo, T, Rsb):
    c1 = 5.38088 * 10 ** -3
    c2 = -1.87784
    c3 = 3.1437
    c4 = 1.32657
    c5 = 1 / 0.715082
        
    if P >= Pb: # Logical condition
        Rs = Rsb
    else:
        Rs = (P / (c1 * Yg ** c2 * Yo ** c3 * T ** c4)) ** c5
    RsAlMarhoun = Rs
    return RsAlMarhoun

def RsDeGhettoCorrelation(Yg, Pb, P, API, T, Tsep, Psep, Rsb):
    if P >= Pb: # Logical condition
        Rs = Rsb
    else:
        if API <= 10: # Extra-heavy oil
            Rs = Yg * ((P / 10.7025) * 10 **(0.0169 * API - 0.00156 * T)) ** 1.1128
        else:
            YgCorr = Yg * (1 + 0.5912 * API * Tsep * log10(Psep / 114.7) * 10 ** -4)
            Rs = ((YgCorr * P ** 1.2057) / 56.434) * 10 ** (10.9267 * API / (T + 460))
    RsDeGhetto = Rs
    return RsDeGhetto
        

def RsDindorukChristmanCorrelation(Yg, Pb, P, API, T, Rsb):
    if P >= Pb: # Logical condition
        Rs = Rsb
    else:    
        a1 = 4.86996 * 10 ** -6
        a2 = 5.730982539
        a3 = 9.9251 * 10 ** -3
        a4 = 1.776179364
        a5 = 44.2500268
        a6 = 2.702889206
        a7 = 0.744335673
        a8 = 3.35975497
        a9 = 28.10133245
        a10 = 1.57905016
        a11 = 0.928131344
        A = (a1 * API ** a2 + a3 * T ** a4) / (a5 + (2 * API ** a6) / (Pb ** a7)) ** 2
        Rs = ((P / a8 + a9) * Yg ** a10 * 10 ** A) ** a11
    RsDindorukChristman = Rs
    return  RsDindorukChristman

def RsDoklaOsmanCorrelation(Yg, Pb, P, Yo, T, Rsb):
    if P >= Pb: # Logical condition
        Rs = Rsb
    else: 
        Rs = ((1/0.836386 * 10 **4) * P * Yg ** 1.01049 * Yo ** -0.107991 * T ** 0.952584) ** (1/0.724047)
    RsDoklaOsman = Rs
    return RsDoklaOsman       

def RsGlasoCorrelation (Yg, Pb, P, API, T, Rsb):
    if P >= Pb: # Logical condition
        Rs = Rsb
    else:
        a = -0.30218
        b = 1.7447
        c = 1.7669 - log10(P)
        if (b ** 2 - 4 * a * c) < 0: # To avoid calculating the root of a negative number
            R = -b / (2 * a)
        else:
            R = (-b + (b ** 2 - 4 * a * c) ** 0.5) / (2 * a)
        F = 10 ** R
        Rs = Yg*((F * API ** 0.989) / (T ** 0.172)) ** (1 / 0.816)
        RsGlaso = Rs
        return RsGlaso


def RsKartoatmodjoSchmidtCorrelation(Yg, Pb, P, API, T, Tsep, Psep, Rsb):
    if P >= Pb: # Logical condition
        Rs = Rsb
    else:
        c = 1 + 0.1595 * API ** 0.4078 * Tsep ** -0.2466 * log10(Psep / 114.7)
        YgCorr = c * Yg 
        if API > 30:
            c1 = 0.0315
            c2 = 0.7587
            c3 = 1.0937
            c4 = 11.289
        else:
            c1 = 0.05958
            c2 = 0.7972
            c3 = 1.0014
            c4 = 13.1405
    Rs = c1 * (YgCorr) **c2 * (P) ** c3 * (10) ** (c4 * API / (T + 460) )
    RsKartoatmodjoSchmidt = Rs
    return RsKartoatmodjoSchmidt
    
def RsLasaterCorrelation(Yg, Pb, P, Yo, T, Rsb):
    if P >= Pb: # Logical condition
        Rs = Rsb
    else:
        API = 141.5 / Yo - 131.5
        if API <= 40:
            Mo = 630 - 10 * API
        else:
            Mo = 73110 * API ** -1.562
        # Calculation of Xpb = Bubble point pressure factor (P * Yg) / T
        Xpb = P * Yg / T
        if Xpb < 3.29:
            c1 = 0.359
            c2 = ((1.473 * P * Yg) / T)
            c3 = 0.476
            Sy = c1 * log(c2 + c3)
        else:
            c1 = ((0.121 * P * Yg) / T)
            c2 = -0.236
            c3 = 0.281
            Sy = (c1 + c2) ** c3
                
        if Sy>=1:
            Sy=0.99999999
        Rs = (132755 * Yo * Sy) / (Mo * (1 - Sy))    
    
        if Rs < 0: # Logical condition            
            Rs = 0
    
        RsLasater = Rs
    
        return RsLasater


def RsPetroskyFarshadCorrelation(Yg, Pb, P, API, T, Rsb):
    if P >= Pb: # Logical condition
        Rs = Rsb
    else:
        X = 7.916 * 10 ** -4 * API ** 1.541 - 4.561 * 10 ** -5 * T ** 1.3911
        Rs = ((P / 112.727 + 12.34) * Yg ** 0.8439 * 10 ** X) ** 1.73184
    RsPetroskyFarshad = Rs
    return RsPetroskyFarshad


def RsStandingCorrelation(Yg, Pb, P, API, T, Rsb):
    if P >= Pb: # Logical condition
        Rs = Rsb
    else:
        Rs = Yg * (((P / 18.2) + 1.4) * 10 ** (0.0125 * API - 0.00091 * T)) ** 1.2048  
    if Rs < 0: # Logical condition
        Rs = 0
    RsStanding = Rs
    return RsStanding

def RsTotalCFPCorrelation(Yg, Pb, P, API, T, Rsb):
    if P >= Pb: # Logical condition
        Rs = Rsb
    else:
        if API <= 10:
            c1 = 12.2651
            c2 = 0.030405
            c3 = 0
            c4 = 0.9669
            
        if 10 < API <= 35:
            c1 = 15.0057
            c2 = 0.0152
            c3 = 4.484 * 10 ** -4
            c4 = 1.095
        
        if 35 < API <= 45:
            c1 = 112.925
            c2 = 0.0248
            c3 = -1.469 * 10 ** -3
            c4 = 1.129
        
        Rs = Yg * ((P / c1) * (10) ** ((c2 * API) - (c3 * T))) ** c4
        
    if Rs < 0: # Logical condition
        Rs = 0 
    RsTotalCFP = Rs
    return RsTotalCFP

def RsVasquezBeggsCorrelation(Yg, Pb, P, API, T, Tsep, Psep, Rsb):
    if P >= Pb: # Logical condition
        Rs = Rsb
    else:
        if API <= 30:
            c1 = 0.0362
            c2 = 1.0937
            c3 = 25.724
        else:
            c1 = 0.0178
            c2 = 1.187
            c3 = 23.931
        Ygs= Yg*(1.+5.912*10**-5*(API)*(Tsep)*log10(Psep/114.7))
        Rs = c1 * Ygs * P ** c2 * exp((c3 * API) / (T + 460))
    
    if Rs < 0:
        Rs = 0 
    
    RsVasquezBeggs = Rs
    
    return RsVasquezBeggs


def RsVelardeCorrelation(Yg, Pb, P, API, T, Rsb):
    if P >= Pb: # Logical condition
        Rs = Rsb
    else:
        a0 = 9.73 * 10 ** -7
        a1 = 1.672608
        a2 = 0.92987
        a3 = 0.247235
        a4 = 1.056052
        
        b0 = 0.022339
        b1 = -1.00475
        b2 = 0.337711
        b3 = 0.132795
        b4 = 0.302065
        
        c0 = 0.725167
        c1 = -1.48548
        c2 = -0.164741
        c3 = -0.09133
        c4 = 0.047094
        
        S1 = a0 * Yg ** a1 * API ** a2 * T ** a3 * Pb ** a4
        S2 = b0 * Yg ** b1 * API ** b2 * T ** b3 * Pb ** b4
        S3 = c0 * Yg ** c1 * API ** c2 * T ** c3 * Pb ** c4
        
        Pr = 1.0*P / Pb
        Rsr = S1 * Pr ** S2 + (1 - S1) * Pr ** S3
        Rs = Rsr * Rsb    
    
    if Rs < 0: # Logical condition            
        Rs = 0
    
    RsVelarde = Rs
    
    return RsVelarde

def RsCegarraCorrelation(Yg, Pb, P, API, T, Rsb):
    if P >= Pb: # Logical condition
        Rs = Rsb
    else:
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
        Rs = Yg * ((((P/l1) + l5) * (10 ** ((l4 * API) - (l3 * T)))) ** (1/l2))
    if Rs < 0: # Logical condition
        Rs = 0 
    RsCegarra = Rs
    return RsCegarra

def RsPerezMLCorrelation(Yg, Pb, P, API, T, Rsb):
    if P >= Pb: # Logical condition
        Rs = Rsb
    else:
        Rs = Rsb * ((1.0*P/Pb) ** 0.881)
    if Rs < 0: # Logical condition
        Rs = 0 
    RsPerezML = Rs
    return RsPerezML

def RsMillanArciaCorrelation(Yg, Pb, P, API, T, Rsb):
    if P >= Pb: # Logical condition
        Rs = Rsb
    else:
        Rs = (Rsb/1.031) * ((1.0*P/Pb) ** 0.83)
    if Rs < 0: # Logical condition
        Rs = 0 
    RsMillanArcia = Rs
    return RsMillanArcia
        

def RsManucciRosalesCorrelation(Yg, Pb, P, API, T, Rsb):
    if P >= Pb: # Logical condition
        Rs = Rsb
    else:
        X = (0.000922 * T) - (0.0072 * API)
        Rs = Yg * ((P/(84.88 * (10 ** X))) ** (1/0.53))
        if Rs < 0: # Logical condition
            Rs = 0 
    RsManucciRosales = Rs
    return RsManucciRosales

def BoAlmarhounCorrelation(Yg, Yo, Rs, Rsb, T, P, Pb, Co):
    if P < Pb: # Saturated oil
        F = Rs ** 0.74239 * Yg ** 0.323294 * Yo ** -1.20204
        Bo = 0.497069 + 0.862963 * 10 ** -3 * T + 0.182594 * 10 ** -2 * F + 0.318099 * 10 ** -5 * F ** 2 
    else: # Undersaturated oil
        Fb = Rsb ** 0.74239 * Yg ** 0.323294 * Yo ** -1.20204
        Bob = 0.497069 + 0.862963 * 10 ** -3 * T + 0.182594 * 10 ** -2 * Fb + 0.318099 * 10 ** -5 * Fb ** 2
        Bo = Bob * exp(Co * (Pb - P))
        if Bo < 1:
            Bo = 1            
        
    BoAlmarhoun = Bo
    return BoAlmarhoun 


def BoAlShammasiCorrelation(Yg, Yo, Rs, Rsb, T, P, Pb, Co):
    if P < Pb: # Saturated oil
        Bo = 1 + (5.53 * 10 ** -7) * (Rs * (T - 60)) + 0.000181 * (Rs / Yo) + 0.000449 * ((T - 60) / Yo) + 0.000206 * (Rs * Yg / Yo)
    else: # Undersaturated oil
        Bob = 1 + (5.53 * 10 ** -7) * (Rsb * (T - 60)) + 0.000181 * (Rsb / Yo) + 0.000449 * ((T - 60) / Yo) + 0.000206 * (Rsb * Yg / Yo)
        Bo = Bob * exp(Co * (Pb - P))
    if Bo < 1:
        Bo = 1
    BoAlShammasi = Bo
    return BoAlShammasi 
    
def BoAlShammasiCorrelation(Yg, Yo, Rs, Rsb, T, P, Pb, Co):
    if P < Pb: # Saturated oil
        Bo = 1 + (5.53 * 10 ** -7) * (Rs * (T - 60)) + 0.000181 * (Rs / Yo) + 0.000449 * ((T - 60) / Yo) + 0.000206 * (Rs * Yg / Yo)
    else: # Undersaturated oil
        Bob = 1 + (5.53 * 10 ** -7) * (Rsb * (T - 60)) + 0.000181 * (Rsb / Yo) + 0.000449 * ((T - 60) / Yo) + 0.000206 * (Rsb * Yg / Yo)
        Bo = Bob * exp(Co * (Pb - P))
            
    if Bo < 1:
        Bo = 1
           
    BoAlShammasi = Bo
    return BoAlShammasi 

def BoDeGhettoCorrelation(Yg, API, Rs, Rsb, T, Tsep, P, Pb, Psep, Co):
    c1 = 4.677 * 10 ** -4
    c2 = 1.751 * 10 ** -5
    c3 = -1.811 * 10 ** -8
        
    Ygs = Yg * (1. + 5.912 * (10**-5) * (API) * (Tsep) * log10(Psep/114.7))
                         
    if P < Pb: # Saturated oil
        Bo = 1 + c1 * Rs + c2 * (T - 60) * (API / Ygs) + c3 * Rs * (T - 60) * (API / Ygs)
    else: # Undersaturated oil
        Bob = 1 + c1 * Rsb + c2 * (T - 60) * (API / Ygs) + c3 * Rsb * (T - 60) * (API / Ygs)
            
        Bo = Bob * exp(Co * (Pb - P))
                    
    if Bo < 1:
        Bo = 1
    BoDeGhetto = Bo
    return BoDeGhetto


def BoDindorukChristmanCorrelation(Yg, API, Rs, Rsb, T, Tsep, P, Pb, Co):
    a1 = 2.510755 * 10 ** 0
    a2 = -4.852538 * 10 ** 0
    a3 = 1.1835 * 10 ** 1
    a4 = 1.365428 * 10 ** 5
    a5 = 2.25288 * 10 ** 0
    a6 = 1.00719 * 10 ** 1
    a7 = 4.450849 * 10 ** -1
    a8 = 5.352624 * 10 ** 0
    a9 = -6.309052 * 10 ** -1
    a10 = 9.000749 * 10 ** -1
    a11 = 9.871766 * 10 ** -1
    a12 = 7.865146 * 10 ** -4
    a13 = 2.689173 * 10 ** -6
    a14 = 1.100001 * 10 ** -5
    
    b1 = 4.236114474
    b2 = 24.316998249
    b3 = 0.958319868
    b4 = 0.924700438
    b5 = 1.287177430
    b6 = 1.353868836
    b7 = 12.581487761
    b8 = 9.828286832
    
    Yo = 141.5 / (131.5 + API)
    
    if P < Pb: # Saturated oil
        A = (((Rs ** a1 * Yg ** a2) / Yo ** a3) + a4 * (T - 60) ** a5 + a6 * Rs) ** a7 / (a8 + (2 * Rs ** a9 / Yg ** a10) * (T - 60)) ** 2
        Bo = a11 + a12 * A + a13 * A ** 2 + a14 * (T - 60) * API / Yg
        if Bo > 2:
            BoDL = Bo
            n = (T - Tsep) ** b1 * (log10(BoDL) * tanh(BoDL)) ** b2 + b3 * (BoDL - 1) ** b4
            d = (1 + BoDL ** b5 * (T - Tsep) ** b6 * (log10(BoDL)) ** b7 ) ** b8
            Bo = 1 + (n / d)
    else: # Undersaturated oil
        Ab = (((Rsb ** a1 * Yg ** a2) / Yo ** a3) + a4 * (T - 60) ** a5 + a6 * Rsb) ** a7 / (a8 + (2 * Rsb ** a9 / Yg ** a10) * (T - 60)) ** 2
        Bob = a11 + a12 * Ab + a13 * Ab ** 2 + a14 * (T - 60) * API / Yg
        if Bob > 2:
            BoDL = Bob
            n = (T - Tsep) ** b1 * (log10(BoDL) * tanh(BoDL)) ** b2 + b3 * (BoDL - 1) ** b4
            d = (1 + BoDL ** b5 * (T - Tsep) ** b6 * (log10(BoDL)) ** b7 ) ** b8
            Bob = 1 + (n / d)
        Bo = Bob * exp(Co * (Pb - P))
    if Bo < 1:
        Bo = 1
    BoDindorukChristman = Bo
    return BoDindorukChristman
    
def BoDoklaOsmanCorrelation(Yg, Yo, Rs, Rsb, T, P, Pb, Co):
    if P < Pb: # Saturated oil
        M = Rs ** 0.773572 * Yg ** 0.40402 * Yo ** -0.882605
        Bo = 0.431935 * 10 ** -1 + 0.156667 * 10 ** -2 * T + 0.139775 * 10 ** -2 * M + 0.380525 * 10 ** -5 * M ** 2
    else: # Undersaturated oil
        Mb = Rsb ** 0.773572 * Yg ** 0.40402 * Yo ** -0.882605
        Bob = 0.431935 * 10 ** -1 + 0.156667 * 10 ** -2 * T + 0.139775 * 10 ** -2 * Mb + 0.380525 * 10 ** -5 * Mb ** 2
        Bo = Bob * exp(Co * (Pb - P))
    BoDoklaOsman = Bo
    return BoDoklaOsman 

def BoGlasoCorrelation(Yg, Yo, Rs, Rsb, T, P, Pb, Co):
    if P < Pb: # Saturated oil
        F = Rs * (Yg / Yo) ** 0.526 + 0.968 * T
        Bo = 1 + 10 ** (-6.58511 + 2.91329 * log10(F) - 0.27683 * (log10(F)) ** 2)
    else: # Undersaturated oil
        Fb = Rsb * (Yg / Yo) ** 0.526 + 0.968 * T
        Bob = 1 + 10 ** (-6.58511 + 2.91329 * log10(Fb) - 0.27683 * (log10(Fb)) ** 2)
    Bo = Bob * exp(Co * (Pb - P))
    BoGlaso = Bo
        
    return BoGlaso

def BoHanafyCorrelation(Rs, Rsb, P, Pb, Co):
    if P < Pb: # Saturated oil  
        Bo = 0.0006 * Rs + 1.079
    else: # Undersaturated oil
        Bob = 0.0006 * Rsb + 1.079
        Bo = Bob * exp(Co * (Pb - P))
    BoHanafy = Bo
    return BoHanafy

def BoKartoatmodjoCorrelation(Yg, Yo, Rs, Rsb, T, Tsep, P, Pb, Psep, Co):
    API = 141.5 / Yo - 131.5
    c = 1 + 0.1595 * (API ** 0.4078) * (Tsep ** -0.2466) * log10(Psep / 114.7)
    Ygcorr = c * Yg
    if P < Pb: # Saturated oil
        F = (Rs ** 0.755) * (Ygcorr ** 0.25) * (Yo ** -1.5) + 0.45 * T
        Bo = 0.98496 + 0.0001 * F ** 1.5
    else: # Undersaturated oil
        Fb = (Rsb ** 0.755) * (Ygcorr ** 0.25) * (Yo ** -1.5) + 0.45 * T
        Bob = 0.98496 + 0.0001 * Fb ** 1.5
        Bo = Bob * exp(Co * (Pb - P))
        
    BoKartoatmodjoSchmidt = Bo
    return BoKartoatmodjoSchmidt

def BoLasaterCorrelation(Yg, Yo, Rs, Rsb, T, P, Pb, Co):
    if P < Pb: # Saturated oil
        Bo = 0.972 + 0.000147 * (Rs * (Yg / Yo) ** 0.5 + 1.25 * T) ** 1.175 
    else: # Undersaturated oil
        Bob = 0.972 + 0.000147 * (Rsb * (Yg / Yo) ** 0.5 + 1.25 * T) ** 1.175 
        Bo = Bob * exp(Co * (Pb - P))
    BoLasater = Bo
    return BoLasater

def BoPetroskyFarshadCorrelation(Yg, Yo, Rs, Rsb, T, P, Pb, Co):
    if P < Pb: # Saturated oil
        F = (Rs ** 0.3738 * ((Yg ** 0.2914) / (Yo ** 0.6265)) + 0.24626 * T ** 0.5371)
        Bo = 1.0113 + 7.2046 * 10 ** -5 * F ** 3.0936
    else:  # Undersaturated oil
        Fb = (Rsb ** 0.3738 * ((Yg ** 0.2914) / (Yo ** 0.6265)) + 0.24626 * T ** 0.5371)
        Bob = 1.0113 + 7.2046 * 10 ** -5 * Fb ** 3.0936
        Bo = Bob * exp(Co * (Pb - P))
    BoPetroskyFarshad = Bo
    return BoPetroskyFarshad

def BoStandingCorrelation(Yg, Yo, Rs, Rsb, T, P, Pb, Co):
    if P < Pb: # Saturated oil
        Bo = 0.972 + 0.000147 * (Rs * (Yg / Yo) ** 0.5 + 1.25 * T) ** 1.175 
    else: # Undersaturated oil
        Bob = 0.972 + 0.000147 * (Rsb * (Yg / Yo) ** 0.5 + 1.25 * T) ** 1.175 
        Bo = Bob * exp(Co * (Pb - P))
    BoStanding = Bo
    return BoStanding

def BoTotalCFPCorrelation(Yg, API, Rs, Rsb, T, P, Pb, Co):
    if P < Pb: # Saturated oil
        Bo = 1.022 + 4.857 * 10 ** -4 * Rs - 2.009 * 10 ** -6 * (T - 60) * (API / Yg) + 17.569 * 10 ** -9 * Rs * (T - 60) * (API / Yg)
    else: # Undersaturated oil
        Bob = 1.022 + 4.857 * 10 ** -4 * Rsb - 2.009 * 10 ** -6 * (T - 60) * (API / Yg) + 17.569 * 10 ** -9 * Rsb * (T - 60) * (API / Yg)
        Bo = Bob * exp(Co * (Pb - P))
    BoTotalCFP = Bo
    return BoTotalCFP

def BoVasquezBeggsCorrelation(Yg, API, Rs, Rsb, T, Tsep, P, Pb, Psep, Co):
    if API <= 30:
        c1 = 4.677 * 10 ** -4
        c2 = 1.751 * 10 ** -5
        c3 = -1.811 * 10 ** -8
    else:
        c1 = 4.67 * 10 ** -4
        c2 = 1.1 * 10 ** -5
        c3 = 1.337 * 10 ** -9
    Ygs = Yg * (1. + 5.912 * (10 ** -5) * (API) * (Tsep) * log10(Psep/114.7))
    if P < Pb: # Saturated oil
        Bo = 1 + c1 * Rs + c2 * (T - 60) * (API / Ygs) + c3 * Rs * (T - 60) * (API / Ygs)
    else: # Undersaturated oil
        Bob = 1 + c1 * Rsb + c2 * (T - 60) * (API / Ygs) + c3 * Rsb * (T - 60) * (API / Ygs)
        Bo = Bob * exp(Co * (Pb - P))
    BoVasquezBeggs = Bo
    return BoVasquezBeggs

def BoCegarraCorrelation(Yg, Yo, Rs, Rsb, T, P, Pb, Co):
    if P < Pb: # Saturated oil
        Bo = 0.972 + 0.000147 * (Rs * (Yg / Yo) ** 0.5 + 1.25 * T) ** 1.175 
    else: # Undersaturated oil
        Bob = 0.972 + 0.000147 * (Rsb * (Yg / Yo) ** 0.5 + 1.25 * T) ** 1.175 
        Bo = Bob * exp(Co * (Pb - P))
    BoCegarra = Bo
    return BoCegarra


def BoPerezMLCorrelation(Yg, Yo, Rsb, T, P, Pb, Co):
    C = (10 ** (-4)) * (((Rsb) * ((Yg/Yo) ** 0.5)) + (1.25 * T))
    Bob = 0.974322 + (4.35048 * C) - (2.85869 * (C ** 2))
    if P < Pb: # Saturated oil
        Bo = Bob * ((1) - ((8.801 + (0.1547 * T)) * (1e-3) * (1 - (1.0*P/Pb))))
    else: # Undersaturated oil
        Bo = Bob * exp(Co * (Pb - P))
    BoPerezML = Bo
    return BoPerezML

def BoVelardeMcCainCorrelation(Yg, API, Rs, Rsb, T, P, Pb, Co):
    pwater = 62.4 # The specific weight of water is 62.4 lb/ft³ or lb/scf
    psto = pwater * Yo
    A = 0 
    B = 52.8 - 0.01 * Rsb 
    while (fabs(A - B) > 0.00001):
        A = B
        ppo = A
        pa = -49.893 + 85.0149 * Yg - 3.70373 * Yg * ppo + 0.047981 * Yg * ppo ** 2 + 2.98914 * ppo- 0.035688 * ppo ** 2
        ppof = (Rs * Yg + 4600 * Yo) / (73.71 + Rs * Yg / pa)
        B = ppof     
        
    ppo = ppof
        
    if P < Pb: # Saturated oil
        pbs = ppo + (0.167 + 16.181 * (10 ** (-0.0425 * ppo))) * (P / 1000) - 0.01 * (0.299 + 263 * (10 ** (-0.0603 * ppo))) * (P / 1000) ** 2
        poR = pbs - (0.00302 + 1.505 * pbs ** -0.951) * (T - 60) ** 0.938 + (0.0233 * (10 ** (-0.0161 * pbs))) * (T - 60) ** 0.475
        Bo = (psto + 0.01357 * Rs * Yg) / (poR)
            
    else: # Undersaturated oil
        ppob = ppo
        pbsb = ppob + (0.167 + 16.181 * (10 ** (-0.0425 * ppob))) * (Pb / 1000) - 0.01 * (0.299 + 263 * (10 ** (-0.0603 * ppob))) * (Pb / 1000) ** 2
        poRb = pbsb - (0.00302 + 1.505 * pbsb ** -0.951) * (T - 60) ** 0.938 + (0.0233 * (10 ** (-0.0161 * pbsb))) * (T - 60) ** 0.475
        Bob = (psto + 0.01357 * Rsb * Yg) / (poRb)
        Bo = Bob * exp(Co * (Pb - P))
    BoVelarde = Bo
    return BoVelarde

def BoMillanArciaCorrelation(Yg, API, Rs, Rsb, T, P, Pb, Co):
    Bob = 1.3537 * (Rsb ** 0.1505) * (Pb ** (-0.1239)) * exp(-0.00405 * API)
    if P < Pb: # Saturated oil
        Bo = Bob * (0.9419 + (0.0608 * (1.0*P/Pb)))
    else: # Undersaturated oil
        Bo = Bob * exp(Co * (Pb - P))
    BoMillanArcia = Bo
    return BoMillanArcia

def BoManucciRosalesCorrelation(API, Rs, Rsb, T, P, Pb, Co):
    if P < Pb: # Saturated oil
        X = (3.27 * 1e-4 * T) + (0.00321 * API)
        Bo = 0.2378 * (Rs ** 0.221) * (10 ** X)
    else: # Undersaturated oil
        X = (3.27 * 1e-4 * T) + (0.00321 * API)
        Bob = 0.2378 * (Rsb ** 0.221) * (10 ** X)
        Bo = Bob * exp(Co * (Pb - P))
    BoManucciRosales = Bo
    return BoManucciRosales