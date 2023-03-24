#Put all helper functions here.
import data_dict
import math
#import ADM

#Runs parameters through set formulas and changes GEF values
def ComplexDiv(A, B, C, D):
    G = 1 / (C * C + D * D)
    E = (A * C + B * D) * G
    F = (B * C - A * D) * G
    return E, F

#Runs parameters through set formulas and changes EF values
def ComplexMul(A, B, C, D):
    E = A * C - B * D
    F = A * D + B * C
    return E, F

def Polar(X, Y):
    R = math.sqrt(X * X + Y * Y)
    Th = math.atan2(X, Y)
    return R, Th

def inverse_distance(detection_distance, prop_loss_cum, prop_loss_indiv):
    m_measure_distance = data_dict.given_cons['m_measure_distance']

    Log10Div10 = 0.230258509
    TenDivLog10 = 1 / Log10Div10
    #D4 should be mic distance source m. changes depending on what macro is used. Can also be detect distance
    p_inv = 2 * TenDivLog10 * math.log(m_measure_distance / detection_distance) #Inverse distance loss from D3 to D4. D3 is the mic distance from the target from excel.

    for i in range(24):
        if (p_inv == 0):
            prop_loss_indiv[i] = -0.001
        prop_loss_cum[i] = prop_loss_indiv[i]
    
    return prop_loss_cum, prop_loss_indiv

def Ingard(freq):
    Hs = data_dict.det_cons['source_height_det'] 
    Hr = data_dict.det_cons['listener_height_det'] 
    R = data_dict.given_cons['m_measure_distance']
    celsius_degrees_T1 = data_dict.det_cons['celsius_degrees_det']
    kelvin_Too = 273.15
    Co = 331.32
    Rhoa = 0.747
    Rhod = 0.747
    L = 1.1
    Pi = 4 * math.atan(1)
    Em2 = data_dict.det_cons['em2_det']
    Sigma = data_dict.meas_cons['sigma_meas']
    Log10Div10 = 0.230258509
    TenDivLog10 = 1 / Log10Div10
    ground_effect = []
  
    HspHr = Hs + Hr
    HsmHr = Hs - Hr
    R1 = math.sqrt(R * R + HsmHr * HsmHr)
    R2 = math.sqrt(R * R + HspHr * HspHr)
    R12 = R1 / R2
    Dr = R2 - R1
    SinTheta1 = HspHr / R2
    CosTheta1 = R / R2
    TanTheta1 = HspHr / R
    Tk = kelvin_Too + celsius_degrees_T1
    Cs = Co * math.sqrt(Tk / kelvin_Too)
    Bef = math.exp(math.log(2) / 6)
    Mu = (Bef - 1 / Bef) / 2
    Eta = (Bef + 1 / Bef) / 2
    I1c = math.sqrt(Pi) * Em2 * R1 * L

    for i in range(24):
        K1 = 2 * Pi * freq[i] / Cs
        K1dr = K1 * Dr
        Del = R1 / (K1 * L * L)
        Omega = math.sqrt(1 + 1 / (Del * Del)) - 1
        Dom1 = Del * Omega
        Dom2 = Del * math.sqrt(2 * Omega)
        At = math.atan(Dom1 / (1 - Dom2)) - math.atan(Dom1 / (1 + Dom2))
        L1 = 0.5 * Dom1 * math.log((1 + Dom2) / abs(1 - Dom2)) + At
        I1 = I1c * K1 * K1
        I2 = 0.5 * I1 * L1 / ((Dom1 + Del) * Dom2)
        Sd2 = 0.5 * (I1 + I2)
        X = 0.5 * (I1 - I2)
        if(X > 1):
            Ea2 = 0.27 * math.exp(math.log(X) * 0.33)
        else:
            Ea2 = X / (1 + 11 * X / 4)

        #Functions we need for this portion: ComplexDiv, ComplexMul, Polar
        Lnfos = math.log(freq[i] / Sigma)
        Rrc = 1 + 9.08 * math.exp(-0.75 * Lnfos)
        Xrc = 11.9 * math.exp(-0.73 * Lnfos)
        Rzr, Rzi = ComplexDiv(1, 0, Rrc, Xrc)
        Rpr, Rpi = ComplexDiv(SinTheta1 - Rzr, -Rzi, SinTheta1 + Rzr, Rzi)
        Srzr, Srzi = ComplexMul(SinTheta1 + Rzr, Rzi, SinTheta1 + Rzr, Rzi)
        Czr, Czi = ComplexDiv(Srzr, Srzi, 1 + SinTheta1 * Rzr, SinTheta1 * Rzi)
        Wr, Wi = ComplexMul(0, 0.5 * K1 * R2, Czr, Czi)
        Wm, Wp = Polar(Wr, Wi) #Means Wm = Sqr((Wr * Wr) + (Wi * Wi)) // and Wp = math.atan2(Wr, Wi)
    
        if Wm < 6:
            Intg = 3
            Fact = 1
            Wsr = Wr
            Wsi = Wi
            W1r = Wr
            W1i = Wi

            for J in range(1, 12):
                Fact = Fact * J
                Cons = Fact * Intg
                W2r, W2i = ComplexMul(Wr, Wi, W1r, W1i)
                Wsr = Wsr + W2r / Cons
                Wsi = Wsi + W2i / Cons
                W1r = W2r
                W1i = W2i
                Intg = Intg + 2

            Ewr = math.exp(-Wr)
            Ewrr = Ewr * math.cos(Wi)
            Ewri = -Ewr * math.sin(Wi)
            Whm = math.sqrt(Pi * Wm)
            Whp = Wp / 2
            Whr = Whm * math.cos(Whp + Pi / 2)
            Whi = Whm * math.sin(Whp + Pi / 2)
            Ws1r, Ws1i = ComplexMul(Ewrr, Ewri, Whr - 2 * Wsr, Whi - 2 * Wsi)
            Fr = 1 + Ws1r
            Fi = Ws1i
        else:
            W1r, W1i = ComplexDiv(1, 0, 2 * Wr, 2 * Wi)
            W2r, W2i = ComplexMul(W1r, W1i, W1r, W1i)
            Wsr = W1r + 3 * W2r
            Wsi = W1i + 3 * W2i
            W3r, W3i = ComplexMul(W1r, W1i, W2r, W2i)
            Wsr = Wsr + 15 * W3r
            Wsi = Wsi + 15 * W3i
            Fr = -Wsr
            Fi = -Wsi
        
        Q1r, Q1i = ComplexMul(Fr, Fi, 1 - Rpr, -Rpi)
        Fm, Fp = Polar(Fr, Fi)
        Qr = Rpr + Q1r
        Qi = Rpi + Q1i
        Qm, Qp = Polar(Qr, Qi)
        Mkdr = Mu * K1dr
        Sinc = math.sin(Mkdr) / Mkdr
        Cost = math.cos(Eta * K1dr + Qp)
        Qmr = Qm * R12
        I3 = (1 + Ea2) * (1 + Qmr * Qmr) + 2 * Qmr * (1 + Ea2 * Rhoa) * Cost * Sinc * math.exp(-Sd2 * (1 - Rhod))
        ground_effect.append(TenDivLog10 * math.log(I3))
    
    return ground_effect #is this the right placement for this return?

        
def ground_effect(freq, prop_loss_cum, prop_loss_indiv, ground_effect_ref):
    windspeed = data_dict.det_cons['wind_speed'] #this looks a bit different than the original, not sure if thats ok
    Iwthr1 = 0
    
    
    if (Iwthr1 == 0):
        ground_effect = Ingard(freq)
        if (windspeed >= 0):
            for i in range(24):
                prop_loss_indiv[i] = ground_effect[i] - ground_effect_ref[i]
                if (prop_loss_indiv[i] == 0):
                    prop_loss_indiv[i] = -0.01
                prop_loss_cum[i] = prop_loss_cum[i] + prop_loss_indiv[i]
        else:
            for i in range(24):
                if (prop_loss_indiv[i] == 0):
                    prop_loss_indiv[i] = -0.001
                prop_loss_indiv[i] = ground_effect[i] - ground_effect_ref[i]
        
    return prop_loss_cum, prop_loss_indiv


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

def Barrier(freq, prop_loss_cum, prop_loss_indiv, detection_dist):
    barr_atten = []
    barr_num = data_dict.det_cons['barrier_on']
    distance_from_source = data_dict.det_cons['barrier_dist']
    barrier_height_det = data_dict.det_cons['barrier_dist']
    source_height_det = data_dict.det_cons['source_height_det']
    listener_height_det = data_dict.det_cons['listener_height_det']
    celsius_degrees_T1 = data_dict.det_cons['celsius_degrees_det']
    Pi = 4 * math.atan(1)
    Log10Div10 = 0.230258509
    TenDivLog10 = 1 / Log10Div10
    
    for i in range(24):
        barr_atten.append(-0.001)
    if barr_num > 0:
        if detection_dist >= distance_from_source:
            Hba = barrier_height_det - (source_height_det + distance_from_source * (listener_height_det) / detection_dist)
            X1 = math.sqrt(distance_from_source ** 2 + ((barrier_height_det - source_height_det) **2))
            X1 = X1 + math.sqrt(((detection_dist - distance_from_source) ** 2) + ((barrier_height_det - listener_height_det) ** 2))
            X1 = X1 - math.sqrt((detection_dist **2) + ((listener_height_det - source_height_det)**2))
            Cs = 331.4 * math.sqrt(1 + celsius_degrees_T1 / 273.15)
            
            for i in range(24):
                if(Hba < 0):
                    X2 = 2 * X1 * freq[2] / Cs #not sure about the freq part
                    X3 = math.sqrt(2 * Pi * abs(X2))
                    X2 = -(X2)
                    if(X2 <= -0.1916):
                        barr_atten[i] = -0.01
                    else: 
                        barr_atten[i] = -(5 + 2 * TenDivLog10 * math.log(X3 / math.tan(X3))) - 0.01
                else:
                    X2 = 2 * X1 * freq[i] / Cs #not sure about the freq
                    X3 = math.sqrt(2 * Pi * abs(X2))
                    X2 = math.exp(2 * X3)
                    if(X2 == 1):
                        barr_atten[i] = -5 -0.01
                    else:
                        barr_atten[i] = -(5 + 2 * TenDivLog10 * math.log(X3 * (X2 + 1) / (X2 - 1))) - 0.01
                        if(barr_atten[i] < -20):
                            barr_atten[i] = -20
    
    for i in range(24):
        prop_loss_indiv[i] = barr_atten[i]
        prop_loss_cum[i] = prop_loss_cum[i] + prop_loss_indiv[i]

    return prop_loss_cum, prop_loss_indiv



def Foliage(freq, prop_loss_cum, prop_loss_indiv, detection_distance):
    Fo = [-0.001] * 24
    Fo_list = []
    foliated_zone_nums = data_dict.det_cons['foliage_on']
    distance_frmSource_to_fol = data_dict.det_cons['foliage_dist']
    foliage_depth = data_dict.det_cons['foliage_depth']
    leaf_width = data_dict.det_cons['leaf_width']
    leaf_areapervol = data_dict.det_cons['leaf_areapervol']
    celsius_degrees_det = data_dict.det_cons['celsius_degrees_det']
    Cs = 331.4 * math.sqrt(1 + celsius_degrees_det / 273.15)

    if foliated_zone_nums > 0:
        # If detection distance is greater than the distance from source to edge of foliage
        if detection_distance > distance_frmSource_to_fol:
            # Sets the distance difference to X2
            X2 = detection_distance - distance_frmSource_to_fol
            # If the new detection distance is greater than depth (extent) of foliage in meters
            if X2 > foliage_depth:
                # Set distance the depth (extent) of foliage in meters
                X2 = foliage_depth
            X2 = math.sqrt(X2)
            Cons = 2.647 / math.log(10)

            for I in range(10, 24):
                Ka = (2 * math.pi * freq[I] / Cs) * leaf_width / 100
                if Ka < 0.401:
                    Fo[I] = -0.01
                elif Ka < 5:
                    Fo[I] = -X2 * math.sqrt(leaf_areapervol) * (Cons * math.log(Ka) + 1.05)
                else:
                    Fo[I] = -X2 * math.sqrt(leaf_areapervol) * 2.9
                Fo_list.append(Fo[I])
    else:
        for I in range(0, 24):
            prop_loss_indiv[I] = Fo[I]
            prop_loss_cum[I] = prop_loss_cum[I] + prop_loss_indiv[I]
            Fo_list.append(prop_loss_cum[I])
    
    return prop_loss_cum, prop_loss_indiv


def normal_deviate(p):
    t = math.sqrt(-2 * math.log(p))
    return t - (2.30753 + t * 0.27061) / (1 + t * (0.99229 + t * 0.04481))


def Dprime(bkg, ht, third_obands, B3):
    """
    Functions called in Dprime -> NormalDeviate, mod_background_noise, background_noise_spec, B3,
    one_third_oct_band_weight, max_modBackground_thresh, hearing_threshold, mod_background_noise_2
    
    ***The rest in parameter list are variables
    """
    hit_prob = data_dict.lis_cons['hit_prob']
    false_alarm_rate = data_dict.lis_cons['false_alarm_rate']
    mod_background_noise = [0] * 24
    observer_efficiency = data_dict.lis_cons['observer_efficiency']
    max_modBackground_thresh = [0] * 24
    TenDivLog10 = 10 / math.log(10)
    Log10Div10 = 1/TenDivLog10
    mic_height = data_dict.meas_cons['mic_height_meas']

    if(hit_prob > 0.5):
        Z3 = normal_deviate(1 - hit_prob)
    else:
        Z3 = normal_deviate(hit_prob)
    
    if(false_alarm_rate > 0.5):
        Z = normal_deviate(1 - false_alarm_rate)
    else:
        Z = normal_deviate(false_alarm_rate)
        
    calculate_d = Z3 - Z
    if(hit_prob > 0.5):
        if(false_alarm_rate > 0.5):
            calculate_d = Z3 - Z
        else:
            calculate_d = Z3 + Z
    else:
        calculate_d = Z - Z3
        
    data_dict.lis_cons.update({'d_stat': calculate_d})

    for I in range(24):
        if I > 10:
            mod_background_noise[I] = bkg[I] + TenDivLog10 * math.log(calculate_d / observer_efficiency / B3[I])
        else:
            E4 = 0
            for J in range(5):
                I0 = I + J - 2
                W4 = third_obands[J][I]
                if((W4 > 0) and (I0 >= 0)):
                    E4 += W4 * math.exp(Log10Div10 * bkg[I0])
            mod_background_noise[I] = TenDivLog10 * math.log(E4 * calculate_d / observer_efficiency / B3[I])
            
        if ht[I] > mod_background_noise[I]:
            max_modBackground_thresh[I] = ht[I]
        else:
            max_modBackground_thresh[I] = mod_background_noise[I]

        return max_modBackground_thresh
        

def ansi_humidity(freq):
    rel_hum_H1 = data_dict.det_cons['relative_humid_percent_det']
    R = data_dict.given_cons['m_measure_distance']
    atmos_absorption = []
    kelvin_Too = 273.15
    Co = 331.32 
    celsius_degrees_T1 = data_dict.det_cons['celsius_degrees_det']
    t = kelvin_Too + celsius_degrees_T1
    PSTAR = 1
    acr1 = 0.000000016 * (1.377 * t / (t + 110.4)) / PSTAR
    thto = 2239.1 / t
    exthto = math.exp(-thto)
    rto = thto / (1 - exthto)
    mumaxo = 0.20948 * 4 * math.pi * rto * rto * exthto / 35
    thtn = 3352 / t #in the VBA code (line 57), 3352 has a "!" after it, not sure why
    exthtn = math.exp(-thtn)
    rtn = thtn / (1 - exthtn)
    mumaxn = 0.78084 * 4 * math.pi * rtn * rtn * exthtn / 35
    cs = Co * math.sqrt(t / kelvin_Too)
    h = rel_hum_H1 * math.exp(math.log(10) * (20.5318 - 2939 / t - 2.13759744 * math.log(t))) / PSTAR #H1 is RH% from excel
    fox = (24 + 44100 * h * (0.05 + h) / (0.391 + h)) * math.sqrt(293 / t)    # Is 293 used as 273 + 20degC?
    fni = (9 + 350 * h) * (293 / t)
    d1 = R * 0.01 #R can be a different value depending on what helper function is called

    for i in range(24):
        freq2 = freq[i] #freq is from the Freq. Hz collumn, need to maybe make a df for this?
        freq2 = freq2 * freq2
        acr = acr1 * freq2
        x2o = freq2 / (fox * fox)
        avibo = 868.5 * mumaxo * fox / cs * x2o / (1 + x2o)
        x2n = freq2 / (fni * fni)
        avibn = 868.5 * mumaxn * fni / cs * x2n / (1 + x2n)
        atmos_absorption.append(-d1 * (acr + avibo + avibn)) #AtmAbs, atmosphere absorption, a global variable this is where it is initialized
    
    return atmos_absorption



def atmosphere(freq, prop_loss_cum, prop_loss_indiv, atm_abs, atm_abs_ref):
    ansi_humidity(freq)

    for i in range(24):
        prop_loss_indiv[i] = atm_abs[i] - atm_abs_ref[i] #only place I see AtmAbsRef initialized is in Reference Calc(), but is just iniliatized to AtmAbs[i]

        if (prop_loss_indiv[i] == 0):
            prop_loss_indiv[i] = -0.001
        prop_loss_cum[i] = prop_loss_cum[i] + prop_loss_indiv[i]

    return prop_loss_cum, prop_loss_indiv

        
def signal_noise(target_spectrum, max_modBkg, third_oband, prop_loss_cum):
    detection_band = [0] * 24
    listener_spectrum = [0] * 24
    Log10Div10 = math.log(10) / 10
    TenDivLog10 = 10 / math.log(10)

    for i in range(24):
        if i > 10:
            detection_band[i] = target_spectrum[i] + prop_loss_cum[i] - max_modBkg[i]
        else:
            E3 = 0
            for j in range(5):
                I0 = i + j - 2
                W4 = third_oband[j][i]
                if (W4 > 0 and I0 >= 0):
                    E3 = E3 + W4 * math.exp(Log10Div10 * (target_spectrum[I0] + prop_loss_cum[I0]))
            detection_band[i] = TenDivLog10 * math.log(E3) - max_modBkg[i]
        listener_spectrum[i] = detection_band[i] + max_modBkg[i]
    M2 = detection_band[0]
    B1 = 0
    for i in range(24):
        if (detection_band[i] >= M2):
            M2 = detection_band[i]
            B1 = i
    
    return M2, B1, prop_loss_cum

def binary_search(freq, ground_effect_ref, atm_abs, atm_abs_ref, target, max_modBkg, third_oband):
    Z9 = -1
    m_measure_distance = data_dict.given_cons['m_measure_distance']
    detection_dist = m_measure_distance * 25
    precision_fraction = 0.001
    prop_loss_cum = [0] * 24
    prop_loss_indiv = [0] * 24
    Z9 = Z9 + 1
    M2, B1, prop_loss_cum, prop_loss_indiv = Propagate(freq, detection_dist, ground_effect_ref, atm_abs, atm_abs_ref, target, max_modBkg, third_oband, prop_loss_cum, prop_loss_indiv)

    while M2 >= 0:
        Z9 = Z9 + 1
        D5 = detection_dist
        detection_dist = 2 * detection_dist
        D6 = detection_dist
        M2, B1, prop_loss_cum, prop_loss_indiv = Propagate(freq, detection_dist, ground_effect_ref, atm_abs, atm_abs_ref, target, max_modBkg, third_oband, prop_loss_cum, prop_loss_indiv)
        print(M2)

    if Z9 == 0:
        while M2 <= 0:
            D6 = detection_dist
            detection_dist = detection_dist / 2
            D5 = detection_dist
            M2, B1, prop_loss_cum, prop_loss_indiv = Propagate(freq, detection_dist, ground_effect_ref, atm_abs, atm_abs_ref, target, max_modBkg, third_oband, prop_loss_cum, prop_loss_indiv)

    while abs(D6 - D5) >= precision_fraction * detection_dist:
        detection_dist = (D5 + D6) / 2
        M2, B1, prop_loss_cum, prop_loss_indiv = Propagate(freq, detection_dist, ground_effect_ref, atm_abs, atm_abs_ref, target, max_modBkg, third_oband, prop_loss_cum, prop_loss_indiv)
        if M2 > 0:
            D5 = detection_dist
        else:
            D6 = detection_dist

    #return 5 items - detection distance most important
    return detection_dist


"""
ground_effect_reference(26) is 'Reference Ground Effect during measurement' from .vbs file
"""
def reference_calc(freq):
    ground_effect = Ingard(freq)
    atmos_absorption = ansi_humidity(freq)
    AtmAbsRef = []
    ground_effect_ref = []
        
    for I in range(0,24):
        ground_effect_ref.append(ground_effect[I])
        AtmAbsRef.append(atmos_absorption[I])
    
    return ground_effect, atmos_absorption, ground_effect_ref, AtmAbsRef



def Reverse():
    Propagate()


    
def targetdBA():
    log_10_div_10 = 0.230258509
    ten_divided_by_log_10 = 1 / log_10_div_10

    E = float()
    E = 0.0
    for x in range(23):
        E = E + math.exp(log_10_div_10 * (S1[x] + S10[x]))

    targetdBA_result = ten_divided_by_log_10 * math.log(E)

    return targetdBA_result


def BkgdBA():
     E = float()
     E = 0.0
     
     for I in range(0,23):
         E = E + math.exp(Log10Div10 * (S6[I] + S10[I]))
         
     bkgdba_return = TenDivLog10 * math.log(E)
     
     return bkgdba_return    
 
 
def ListenerdBA():
    E = float()
    E = 0.0
    for I in range(0, 23):
        E = E + math.exp(Log10Div10 * (S1[I] + S3[I] + S10[I]))
    listener_dba_return = TenDivLog10 * math.log(E)
    return listener_dba_return


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

            for x in range(0, 24):
                Ea += Log10Div10 * target_spectrum(x) + prop_loss_cum(x) + A_weight_levels(x)


    while abs(D6 - D5) < precision_fraction * detection_dist:
        detection_dist = (D5 + D6) / 2
        for x in range(0, 24):
            Ea += Log10Div10 * target_spectrum(x) + prop_loss_cum(x) + A_weight_levels(x)

        prop_loss = TenDivLog10 *  log(Ea, 2)

        if prop_loss > dBa:
            D5 = detection_dist
        else:
            D6 = detection_dist


def Propagate(freq, detection_distance, ground_effect_ref, atm_abs, atm_abs_ref, target, max_modBkg, third_oband, prop_loss_cum, prop_loss_indiv):
    prop_loss_cum, prop_loss_indiv = inverse_distance(detection_distance, prop_loss_cum, prop_loss_indiv)
    prop_loss_cum, prop_loss_indiv = ground_effect(freq, prop_loss_cum, prop_loss_indiv, ground_effect_ref)
    prop_loss_cum, prop_loss_indiv = Barrier(freq, prop_loss_cum, prop_loss_indiv, detection_distance)
    prop_loss_cum, prop_loss_indiv = Foliage(freq, prop_loss_cum, prop_loss_indiv, detection_distance)
    prop_loss_cum, prop_loss_indiv = atmosphere(freq, prop_loss_cum, prop_loss_indiv, atm_abs, atm_abs_ref)
    M2, B1, prop_loss_cum = signal_noise(target, max_modBkg, third_oband, prop_loss_cum)

    return M2, B1, prop_loss_cum, prop_loss_indiv
    