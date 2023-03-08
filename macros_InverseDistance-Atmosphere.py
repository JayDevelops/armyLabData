import math

#global variables 
atm_abs = []
atm_abs_ref = []
s8 = []
s3 = []

# def inverse_distance():
#     TEN_DIV_NATURAL_LOG10 = 4.342944824679639
#     #D4 should be mic distance source m. changes depending on what macro is used. Can also be detect distance
#     p_inv = 2 * TEN_DIV_NATURAL_LOG10 * math.log(d3 / d4) #Inverse distance loss from D3 to D4. D3 is the mic distance from the target from excel.

#     for i in range(24):
#         s8[i] = p_inv
#         if (s8[i] == 0):
#             s8[i] = -0.001
#         s3[i] = s8[i]

def ansi_humidity():
    TOO = 273.15
    CO = 331.32
    t = TOO + t1 #T1 is Temp deg C from excel.
    PSTAR = 1
    acr1 = 0.000000016 * (1.377 * t / (t + 110.4)) / PSTAR
    thto = 2239.1 / t
    exthto = math.exp(-thto)
    rto = thto / (1 - exthto)
    mumaxo = 0.20948 * 4 * math.pi * rto * rto * exthto / 35
    thtn = 3352 / t
    exthtn = math.exp(-thtn)
    rtn = thtn / (1 - exthtn)
    mumaxn = 0.78084 * 4 * math.pi * rtn * rtn * exthtn / 35
    cs = CO * math.sqrt(t / TOO) 
    h = h1 * math.exp(math.log10 * (20.5318 - 2939 / t - 2.13759744 * math.log(t))) / PSTAR #H1 is RH% from excel
    fox = (24 + 44100 * h * (0.05 + h) / (0.391 + h)) * math.sqrt(293 / t)    # Is 293 used as 273 + 20degC?
    fni = (9 + 350 * h) * (293 / t)
    d1 = r * 0.01 #R can be a different value depending on what helper function is called

    for i in range(24):
        freq2 = s4[i] #S4 is from the Freq. Hz collumn
        freq2 = freq2 * freq2
        acr = acr1 * freq2
        x2o = freq2 / (fox * fox)
        avibo = 868.5 * mumaxo * fox / cs * x2o / (1 + x2o)
        x2n = freq2 / (fni * fni)
        avibn = 868.5 * mumaxn * fni / cs * x2n / (1 + x2n)
        atm_abs[i] = -d1 * (acr + avibo + avibn) #AtmAbs, atmosphere absorption, a global variable this is where it is initialized


def atmosphere():
    ansi_humidity()

    for i in range(24):
        s8[i] = atm_abs[i] - atm_abs_ref[i] #only place I see AtmAbsRef initialized is in Reference Calc(), but is just iniliatized to AtmAbs[i]

        if (s8[i] == 0):
            s8[i] = -0.001
        s3[i] = s3[i] + s8[i]