import pandas as pd  # Install pandas in python or use Anaconda environment
import data  # python module including helper functions (In our case, translated Macros from ADM by Jeol)
import math

def Detection():
    # data_df = data.read_data()
    # print(data_df.iloc[:,29:34])

    data_df = pd.read_excel('ADM - from Joel - Sept-2013.xls', sheet_name='Data')

    # Instantiating specific columns to be pulled from the data sheet in ADM by Joel
    Trg_name = 'Drone'  # Drone db column
    Bkg_name = 'Urban'  # Ambient setting background noise
    Hth_name = 'ISO Std'  # Hearing Threshold based on ISO Standard

    freq_df = data_df['Freq Hz']
    target_df = data_df[Trg_name]
    bkg_df = data_df[Bkg_name]
    ht_df = data_df[Hth_name]
    awt_weights_df = data_df['Awt weights']
    ai_weights_df = data_df['A.I. weights']
    third_octave_bands_df = data_df.iloc[:, 29:34].copy()

    # #Find and set all necessary dataframes for further processing
    # freq_df, target_df, bkg_noise_df, \
    #     hear_thresh_df, awt_weights, ai_weights_df, \
    #     third_obands_df = data.get_data(Trg_name, Bkg_name, Hth_name, data_df)

    # print(third_obands_df)

# Extract relevant reference values to be used in calcuating detectability. InitMacros()
def get_data(target, background_noise, hearing_th, df):
    mic_distance = df.at[24, target]

    # Drop Extra Rows at the bottom
    df = df.drop(range(24, 27))

    # Extract colums from main df
    freq_df = df['Freq Hz']
    target_df = df[target]
    bkg_df = df[background_noise]
    ht_df = df[hearing_th]
    awt_weights_df = df['Awt weights']
    ai_weights_df = df['A.I. weights']

    # print(freq_df)
    # print(target_df)
    # print(bkg_df)
    # print(ht_df)
    # print(awt_weights_df)
    # print(ai_weights_df)
    # Drop extra rows at the bottom of each df
    # freq_df = freq_df.drop(range(24,27))
    # target_df = target_df.drop(range(24,27))
    # bkg_df = bkg_df.drop(range(24,27))
    # ht_df = ht_df.drop(range(24,27))
    # awt_weights_df = awt_weights_df.drop(range(24,27))
    # ai_weights_df = ai_weights_df.drop(range(24,27))

def atmosphere(s3, s8, atm_abs, atm_abs_ref):
    ansi_humidity()

    for i in range(24):
        s8[i] = atm_abs[i] - atm_abs_ref[i] #only place I see AtmAbsRef initialized is in Reference Calc(), but is just iniliatized to AtmAbs[i]

        if (s8[i] == 0):
            s8[i] = -0.001
        s3[i] = s3[i] + s8[i]

def Ingard():
  HspHr = Hs + Hr
  HsmHr = Hs - Hr
  R1 = math.sqrt(R * R + HsmHr * HsmHr)
  R2 = math.sqrt(R * R + HspHr * HspHr)
  R12 = R1 / R2
  Dr = R2 - R1
  SinTheta1 = HspHr / R2
  CosTheta1 = R / R2
  TanTheta1 = HspHr / R
  Tk = zerDegCelsius_inKelvin + celsius_degrees_det #NOT SURE IF WE NEED TO CHANGE THESE VARIABLE NAMES
  Cs = Co * math.sqrt(Tk / zerDegCelsius_inKelvin)
  Bef = math.exp(math.log(2) / 6)
  Mu = (Bef - 1 / Bef) / 2
  Eta = (Bef + 1 / Bef) / 2
  Rhoa = 0.747
  Rhod = 0.747
  L = 1.1
  I1c = math.sqrt(Pi) * Em2 * R1 * L

  for I in range(24):
    freq = freq[I]
    K1 = 2 * Pi * freq / Cs
    K1dr = K1 * Dr
    Del = R1 / (K1 * L * L)
    Omega = math.sqrt(1 + 1 / (Del * Del)) - 1
    Dom1 = Del * Omega
    Dom2 = Del * math.sqrt(2 * Omega)
    At = math.atan(Dom1 / (1 - Dom2)) - math.atan(Dom1 / (1 + Dom2))
    L1 = 0.5 * Dom1 * math.log((1 + Dom2) / Abs(1 - Dom2)) + At
    I1 = I1c * K1 * K1
    I2 = 0.5 * I1 * L1 / ((Dom1 + Del) * Dom2)
    Sd2 = 0.5 * (I1 + I2)
    X = 0.5 * (I1 - I2)
    if(X > 1):
      Ea2 = 0.27 * math.exp(math.log(X) * 0.33)
    else:
      Ea2 = X / (1 + 11 * X / 4)
    #Functions we need for this portion: ComplexDiv, ComplexMul, Polar
    Lnfos = math.log(freq / Sigma)
    Rrc = 1 + 9.08 * math.exp(-0.75 * Lnfos)
    Xrc = 11.9 * math.exp(-0.73 * Lnfos)
    Done = ComplexDiv(1, 0, Rrc, Xrc, Rzr, Rzi)
    Done = ComplexDiv(SinTheta1 - Rzr, -Rzi, SinTheta1 + Rzr, Rzi, Rpr, Rpi)
    Done = ComplexMul(SinTheta1 + Rzr, Rzi, SinTheta1 + Rzr, Rzi, Srzr, Srzi)
    Done = ComplexDiv(Srzr, Srzi, 1 + SinTheta1 * Rzr, SinTheta1 * Rzi, Czr, Czi)
    Done = ComplexMul(0, 0.5 * K1 * R2, Czr, Czi, Wr, Wi)
    Done = Polar(Wr, Wi, Wm, Wp)
    if(Wm < 6):
      Intg = 3
      Fact = 1
      Wsr = Wr
      Wsi = Wi
      W1r = Wr
      W1i = Wi
      for J in range(1, 12):
        Fact = Fact * J
        Cons = Fact * Intg
        Done = ComplexMul(Wr, Wi, W1r, W1i, W2r, W2i)
        Wsr = Wsr + W2r / Cons
        Wsi = Wsi + W2i / Cons
        W1r = W2r
        W1i = W2i
        Intg = Intg + 2
      Ewr = math.exp(-Wr)
      Ewrr = Ewr * math.cos(Wi)
      Ewri = -Ewr * math.sin(Wi)
      Whm = math.sqr(Pi * Wm)
      Whp = Wp / 2
      Whr = Whm * math.cos(Whp + Pi / 2)
      Whi = Whm * math.sin(Whp + Pi / 2)
      Done = ComplexMul(Ewrr, Ewri, Whr - 2 * Wsr, Whi - 2 * Wsi, Ws1r, Ws1i)
      Fr = 1 + Ws1r
      Fi = Ws1i
    else:
        Done = ComplexDiv(1, 0, 2 * Wr, 2 * Wi, W1r, W1i)
        Done = ComplexMul(W1r, W1i, W1r, W1i, W2r, W2i)
        Wsr = W1r + 3 * W2r
        Wsi = W1i + 3 * W2i
        Done = ComplexMul(W1r, W1i, W2r, W2i, W3r, W3i)
        Wsr = Wsr + 15 * W3r
        Wsi = Wsi + 15 * W3i
        Fr = -Wsr
        Fi = -Wsi
    Done = ComplexMul(Fr, Fi, 1 - Rpr, -Rpi, Q1r, Q1i)
    'Done = Polar(Fr, Fi, Fm, Fp)'
    Qr = Rpr + Q1r
    Qi = Rpi + Q1i
    Done = Polar(Qr, Qi, Qm, Qp)
    Mkdr = Mu * K1dr
    Sinc = math.sin(Mkdr) / Mkdr
    Cost = math.cos(Eta * K1dr + Qp)
    Qmr = Qm * R12
    I3 = (1 + Ea2) * (1 + Qmr * Qmr) + 2 * Qmr * (1 + Ea2 * Rhoa) * Cost * Sinc * math.exp(-Sd2 * (1 - Rhod))
    ground_effect[I] = TenDivLog10 * math.log(I3)
    

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


def inverse_distance(s8, s3, d3, d4):
    _NATURAL_LOG10= (10 / math.log(10, 10))
    #D4 should be mic distance source m. changes depending on what macro is used. Can also be detect distance
    p_inv = 2 * TEN_DIV_NATURAL_LOG10 * math.log(d3 / d4) #Inverse distance loss from D3 to D4. D3 is the mic distance from the target from excel.

    for i in range(24):
        s8[i] = p_inv
        if (s8[i] == 0):
            s8[i] = -0.001
        s3[i] = s8[i]
        
 def ground_effect():
     height_source = source_height_det
     height_listener = listener_height_det
     R = detection_dist
     Sigma = SigmaDet
     Em2 = Em2Det
     TrgHgt = height_source
     DetHgt = height_listener
     windspeed = wind_speed
     Iwthr1 = 0
     if (Iwthr1 == 0):
         # calls Ingard()
         if (windspeed >= 0):
             for i in range(24):
                 prop_loss_indiv[i] = ground_effect_ingard[i] - gournd_effect_initial[i]
                 if (prop_loss_indiv[i] == 0):
                     prop_loss_indiv[i] = -0.01
                 prop_loss_cum[i] = prop_loss_cum[i] + prop_loss_indiv[i]
         else:
             for i in range(24):
                 if (prop_loss_indiv[i] == 0):
                     prop_loss_indiv[i] = -0.001
                 prop_loss_indiv[i] = ground_effect_ingard[i] - ground_effect_initial[i]

"""The variables that are used in this function are described as:

barrier_attenuation -> Global array variable (declared in line 35 of VBA Macros). Inititalized for first time inside Barrier()
detection_dist -> Global 'double' variable (declared in line 13 of VBA Macros). It is also used in GroundEffect() on line 186
distance_from_source -> initialized from excel cell "I5" (line 579)
barrier_height_det -> initialized from excel cell "J5" (line 580)
source_height_det -> initialized from excel cell "H3" (line 564)
listener_height_det -> initialized from excel cell "I3" (line 565)
celsius_degrees_det -> initialized from excel cell "J3" (line 566)
TenDivLog10 -> Static number (described in lines 546 and 547)
freq -> Global integer array declared in line 21. Looks like initialized in InitMacros as 'A' column values of Data Sheet (aka freq Hz)
prop_loss_indiv -> Global 'double' array declared in line 25. Is used in many functions, and changes based on the actions of those functions
prop_loss_cum -> Global 'double' array declared in line 20. Is used to hold value after all changes to prop_loss_indiv throughout the program.  
"""

def Barrier():
    for i in range(24):
        barrier_attenuation(i) = -0.001
    if barrier_number > 0:
        if detection_dist >= distance_from_source:
            Hba = barrier_height_det - (source_height_det + distance_from_source * (listener_height_det) / detection_dist)
            X1 = math.sqrt(distance_from_source ** 2 + ((barrier_height_det - source_height_det) **2))
            X1 = X1 + math.sqrt(((detection_dist - distance_from_source) ** 2) + ((barrier_height_det - listener_height_det) ** 2))
            X1 = X1 - math.sqrt((detection_dist **2) + ((listener_height_det - source_height_det)**2))
            Cs = 331.4 * math.sqrt(1 + celsius_degrees_det / 273.15)
            
            for i in range(24):
                if(Hba < 0):
                    X2 = 2 * X1 * freq(2) / Cs
                    X3 = math.sqrt(2 * Pi * abs(X2))
                    X2 = -(X2)
                    if(X2 <= -0.1916):
                        barrier_attenuation(i) = -0.01
                    else: 
                        barrier_attenuation(i) = -(5 + 2 * TenDivLog10 * math.log(X3 / math.tan(X3))) - 0.01
                else:
                    X2 = 2 * X1 * freq(i) / Cs
                    X3 = math.sqrt(2 * Pi * abs(X2))
                    X2 = math.exp(2 * X3)
                    if(X2 == 1):
                        barrier_attenuation(i) = -5 -0.01
                    else:
                        barrier_attenuation(i) = -(5 + 2 * TenDivLog10 * math.log(X3 * (X2 + 1) / (X2 - 1))) - 0.01
                        if(barrier_attenuation(i) < -20):
                            barrier_attenuation(i) = -20
    
    for i in range(24):
        prop_loss_indiv(i) = barrier_attenuation(i)
        prop_loss_cum(i) = prop_loss_cum(i) + prop_loss_indiv(i)



def binary_search(m_meas_distance, D5, D6, M2, precision_fraction):
    Z9 = -1
    detection_dist = m_meas_distance * 25

    while M2 < 0:
        Z9 = Z9 + 1
        D5 = detection_dist
        detection_dist = 2 * detection_dist
        D6 = detection_dist

    if Z9 == 0:
        while abs(D6 - D5) < precision_fraction * detection_dist:
            detection_dist = (D5 + D6) / 2

            if M2 > 0:
                D5 = detection_dist
            else:
                D6 = detection_dist


"""
ground_effect_reference(26) is 'Reference Ground Effect during measurement' from .vbs file
"""
def reference_calc(ground_effect_ref):
    for I in range(0,23):
        ground_effect_ref[I] = ground_effect(I)
        AtmAbsRef(I) = atmos_absorption(I)


def calculate_measure_dist(detection_dist: float):
    return detection_dist * 25.0


def binary_search_A(dBa, m_meas_distance, precision_fraction, D6, D5):
    z9 = -1
    detection_dist = m_meas_distance * 25
    prop_loss = prop_loss_cum(23)
    Ea = 0

    while prop_loss < dBa:
        D5 = detection_dist
        detection_dist = 2 * detection_dist
        D6 = detection_dist

        for x in range(0, 23):
            Ea += Log10Div10 * target_spectrum(x) + prop_loss_cum(x) + A_weight_levels(x)

    if z9 == 0:
        while prop_loss < dBa:
            D5 = detection_dist
            detection_dist = 2 * detection_dist
            D6 = detection_dist

            for x in range(0, 23):
                Ea += Log10Div10 * target_spectrum(x) + prop_loss_cum(x) + A_weight_levels(x)


    while abs(D6 - D5) < precision_fraction * detection_dist:
        detection_dist = (D5 + D6) / 2
        for x in range(0, 23):
            Ea += Log10Div10 * target_spectrum(x) + prop_loss_cum(x) + A_weight_levels(x)

        prop_loss = TenDivLog10 *  log(Ea, 2)

        if prop_loss > dBa:
            D5 = detection_dist
        else:
            D6 = detection_dist


# Main Function Declaration and Call
def main():
    # Possibly implement a switch case to consider user input
    # of which Macro to call
    Detection()

