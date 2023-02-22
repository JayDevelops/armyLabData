import math


def InverseDistance():
    ten_div_natural_log10 = 4.342944824679639
    Pinv = 2 * ten_div_natural_log10 * math.log(D3 / D4) #Inverse distance loss from D3 to D4

    for i in range(24):
        S8[i] = Pinv
        if (S8[i] == 0):
            s8[i] = -0.001
        S3[i] = S8[i]

def AnsiHumidity():
    T = Too + T1
    Pstar = 1
    Acr1 = 0.000000016 * (1.377 * T / (T + 110.4)) / Pstar
    Thto = 2239.1 / T
    Exthto = math.exp(-Thto)
    Rto = Thto / (1 - Exthto)
    Mumaxo = 0.20948 * 4 * math.pi * Rto * Rto * Exthto / 35
    Thtn = 3352 / T #exclamation point after 3352 should just be to delcare it as a single but not sure
    Exthtn = math.exp(-Thtn)
    Rtn = Thtn / (1 - Exthtn)
    Mumaxn = 0.78084 * 4 * math.pi * Rtn * Rtn * Exthtn / 35
    Cs = Co * math.sqrt(T / Too)
    H = H1 * math.exp(math.log10 * (20.5318 - 2939 / T - 2.13759744 * math.log(T))) / Pstar
    Fox = (24 + 44100 * H * (0.05 + H) / (0.391 + H)) * math.sqrt(293 / T)    # Is 293 used as 273 + 20degC?
    Fni = (9 + 350 * H) * (293 / T)
    D1 = R * 0.01

    for i in range(24):
        Freq2 = S4(I)
        Freq2 = Freq2 * Freq2
        Acr = Acr1 * Freq2
        X2o = Freq2 / (Fox * Fox)
        Avibo = 868.5 * Mumaxo * Fox / Cs * X2o / (1 + X2o)
        X2n = Freq2 / (Fni * Fni)
        Avibn = 868.5 * Mumaxn * Fni / Cs * X2n / (1 + X2n)
        AtmAbs[I] = -D1 * (Acr + Avibo + Avibn)


def Atmosphere():
    AnsiHumidity()

    for i in range(24):
        S8[i] = AtmAbs[i] - AtmAbsRef[i]

        if (S8[i] == 0):
            S8[i] = -0.001
        S3[i] = S3[i] + S8[i]
