import data_dict
import math
import ADM



#Runs parameters through set formulas and changes GEF values
def ComplexDiv(A, B, C, D, E, F):
    G = 1 / (C * C + D * D)
    E = (A * C + B * D) * G
    F = (B * C - A * D) * G
    return 1

#Runs parameters through set formulas and changes EF values
def ComplexMul(A, B, C, D, E, F):
    E = A * C - B * D
    F = A * D + B * C
    return 1

def Polar(X, Y, R, Th):
    R = math.sqrt(X * X + Y * Y)
    Th = math.atan2(X, Y)
    Polar = 1


def InverseDistance():
    prop_loss_indiv = [] #not sure how to use these prop loss
    prop_loss_cum = [] #not sure how to use these prop loss
    freq = [] #not sure how to use this
    m_measure_distance = data_dict.given_cons['m_measure_distance']
    detection_dist = [] #not sure where to get this one from, comes from calculated macros
    
    Log10Div10 = 0.230258509
    TenDivLog10 = 1 / Log10Div10
    p_inv = 2 * TenDivLog10 * math.log(m_measure_distance / detection_dist) #Inverse distance loss from D3 to D4. D3 is the mic distance from the target from excel.

    for i in range(24):
        prop_loss_indiv[i] = p_inv
        if (prop_loss_indiv[i] == 0):
            prop_loss_indiv[i] = -0.001
        prop_loss_cum[i] = prop_loss_indiv[i]
              
        
def ground_effect():
    source_height_det = data_dict.det_cons['source_height_det']
    listener_height_det = data_dict.det_cons['listener_height_det']
    SigmaDet = data_dict.det_cons['sigma_det']
    Em2Det = data_dict.det_cons['em2_det']
    detection_dist = []
    R = detection_dist
    Sigma = SigmaDet
    Em2 = Em2Det
    TrgHgt = source_height_det
    DetHgt = listener_height_det
    windspeed = data_dict.det_cons['wind_speed'] #this looks a bit different than the original, not sure if thats ok
    Iwthr1 = 0
    prop_loss_indiv = [] #not sure how to use these prop loss
    prop_loss_cum = [] #not sure how to use these prop loss
    ground_effect_reference = ground_effect #check this one, using line 468 of original macros
    
    
    if (Iwthr1 == 0):
        Ingard()
        if (windspeed >= 0):
            for i in range(24):
                prop_loss_indiv[i] = ground_effect[i] - ground_effect_reference[i]
                if (prop_loss_indiv[i] == 0):
                    prop_loss_indiv[i] = -0.01
                prop_loss_cum[i] = prop_loss_cum[i] + prop_loss_indiv[i]
        else:
            for i in range(24):
                if (prop_loss_indiv[i] == 0):
                    prop_loss_indiv[i] = -0.001
                prop_loss_indiv[i] = ground_effect[i] - ground_effect_reference[i]


def Barrier():
    barr_atten = []
    barr_num = data_dict.det_cons['barrier_on']
    detection_dist = [] #not sure about this one, since it is calculated inside a macro
    distance_from_source = data_dict.det_cons['barrier_dist']
    barrier_height_det = data_dict.det_cons['barrier_dist']
    source_height_det = data_dict.det_cons['source_height_det']
    listener_height_det = data_dict.det_cons['listener_height_det']
    celsius_degrees_T1 = data_dict.det_cons['celsius_degrees_det']
    Pi = 4 * math.atan(1)
    Log10Div10 = 0.230258509
    TenDivLog10 = 1 / Log10Div10
    prop_loss_indiv = [] #not sure how to use these prop loss
    prop_loss_cum = [] #not sure how to use these prop loss
    freq = [] #not sure how to use this
    
    for i in range(24):
        barr_atten[i] = -0.001
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



def Foliage():
    #-------->>>> needs detection distance and freq fixed up
    fol_atten = [-0.001] * 24
    detection_dist = [] #still not sure where this one comes from
    Fo_list = [] #looks like a local array that Bruck created to hold the prop_loss_cum, not sure if we need it
    freq = [] #not sure about this one (know its just the freq band in col A in Data)
    foliated_zone_nums = data_dict.det_cons['foliage_on']
    distance_frmSource_to_fol = data_dict.det_cons['foliage_dist']
    foliage_depth = data_dict.det_cons['foliage_depth']
    leaf_width = data_dict.det_cons['leaf_width']
    leaf_areapervol = data_dict.det_cons['leaf_areapervol']
    celsius_degrees_det = data_dict.det_cons['celsius_degrees_det']
    prop_loss_indiv = []
    prop_loss_cum = []
    Cs = 331.4 * math.sqrt(1 + celsius_degrees_det / 273.15) #<-------------- I believe this is brought in from Barrier (in the macros, Barrier is called before Foliage, then this Cs variable is used in Foliage without being initialized), but also given a different value in ingard
    
    if foliated_zone_nums > 0:
        # If detection distance is greater than the distance from source to edge of foliage
        if detection_dist > distance_frmSource_to_fol:
            # Sets the distance difference to X2
            X2 = detection_dist  - distance_frmSource_to_fol
            # If the new detection distance is greater than depth (extent) of foliage in meters
            if X2 > foliage_depth:
                # Set distance the depth (extent) of foliage in meters
                X2 = foliage_depth
            X2 = math.sqrt(X2)
            Cons = 2.647 / math.log(10)

            for I in range(10, 24):
                Ka = (2 * math.pi * freq[I] / Cs) * leaf_width / 100
                if Ka < 0.401:
                    fol_atten[I] = -0.01
                elif Ka < 5:
                    fol_atten[I] = -X2 * math.sqrt(leaf_areapervol) * (Cons * math.log(Ka) + 1.05)
                else:
                    fol_atten[I] = -X2 * math.sqrt(leaf_areapervol) * 2.9
                Fo_list.append(fol_atten)
    else:
        for I in range(0, 24):
            prop_loss_indiv[I] = fol_atten[I]
            prop_loss_cum[I] = prop_loss_cum[I] + prop_loss_indiv[I]
            Fo_list.append(prop_loss_cum)
    return Fo_list


def NormalDeviate(p):
    t = math.sqrt(-2 * math.log(p))
    return t - (2.30753 + t * 0.27061) / (1 + t * (0.99229 + t * 0.04481))


def Dprime():
    #is this function supposed to return anything?
    
    hit_prob = data_dict.lis_cons['hit_prob']
    false_alarm_rate = data_dict.lis_cons['hit_prob']
    mod_background_noise = [] # Corresponds with 'Bkg Cb' in calculated Macros (for example, in Detection Distance, its cells M8-M31)
    background_noise_spec = [] #From Data sheet, background noise (AKA Urban environment)
    observer_efficiency = data_dict.lis_cons['observer_efficiency']
    B3 = [] #from Data sheet,  Col AI
    one_third_oct_band_weight = [] #arrays from 'Data' cols AD - AH
    hearing_threshold = data_dict.given_cons['hth_name']
    max_modBackground_thresh = [] #Corresponds with 'Max Bkg/Thr' in calculated macros (ex: cells N8 - N31 in Detection Distance)
    mod_background_noise_2 =[] #Just holds the modified value for max_modBackground_thresh given in line 250
    
    
    
    Log10Div10 = 0.230258509
    TenDivLog10 = 1 / Log10Div10
    prop_loss_indiv = [] #not sure how to use these prop loss
    prop_loss_cum = []

    if(hit_prob > 0.5):
        Z3 = NormalDeviate(1 - hit_prob)
    else:
        Z3 = NormalDeviate(hit_prob)
    
    if(false_alarm_rate > 0.5):
        Z = NormalDeviate(1 - false_alarm_rate)
    else:
        Z = NormalDeviate(false_alarm_rate)
        
    calculate_d = Z3 - Z
    if(hit_prob > 0.5):
        if(false_alarm_rate > 0.5):
            calculate_d = Z3 - Z
        else:
            calculate_d = Z3 + Z
    else:
        calculate_d = Z - Z3
        
    Dprime = calculate_d #Not sure why this is a variable here <-----------
    for i in range(24):
        if(i > 10):
            mod_background_noise[i] = background_noise_spec[i] + TenDivLog10 * math.log(calculate_d / observer_efficiency / B3[i])
        else:
            E4 = 0
            for j in range(5):
                I0 = i + j - 2
                W4 = one_third_oct_band_weight[i, j]
                if((W4 > 0) and (I0 >= 0)):
                    E4 += W4 * math.exp(Log10Div10 * background_noise_spec[I0])
            mod_background_noise[i] = TenDivLog10 * math.log(E4 * calculate_d / observer_efficiency / B3[i])
            
        if(hearing_threshold[i] > mod_background_noise[i]):
            max_modBackground_thresh[i] = hearing_threshold[i]
        else:
            max_modBackground_thresh[i] = mod_background_noise[i]
            
        mod_background_noise_2[i] = max_modBackground_thresh[i]
        
        

def ansi_humidity():
    rel_hum_H1 = data_dict.det_cons['relative_humid_percent_det']
    R = data_dict.given_cons['m_measure_distance']
    atmos_absorption = []
    kelvin_Too = 273.15
    Co = 331.32 
    celsius_degrees_T1 = data_dict.meas_cons['celsius_degrees_det']
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
    h = rel_hum_H1 * math.exp(math.log10 * (20.5318 - 2939 / t - 2.13759744 * math.log(t))) / PSTAR #H1 is RH% from excel
    fox = (24 + 44100 * h * (0.05 + h) / (0.391 + h)) * math.sqrt(293 / t)    # Is 293 used as 273 + 20degC?
    fni = (9 + 350 * h) * (293 / t)
    d1 = R * 0.01 #R can be a different value depending on what helper function is called

    for i in range(24):
        freq2 = ADM.freq[i] #freq is from the Freq. Hz collumn, need to maybe make a df for this?
        freq2 = freq2 * freq2
        acr = acr1 * freq2
        x2o = freq2 / (fox * fox)
        avibo = 868.5 * mumaxo * fox / cs * x2o / (1 + x2o)
        x2n = freq2 / (fni * fni)
        avibn = 868.5 * mumaxn * fni / cs * x2n / (1 + x2n)
        atmos_absorption[i] = -d1 * (acr + avibo + avibn) #AtmAbs, atmosphere absorption, a global variable this is where it is initialized
    
    return atmos_absorption



def Atmosphere(s3, s8, atm_abs, atm_abs_ref):
    ansi_humidity()
    
    prop_loss_indiv = []
    atmos_absorption = []
    atm_abs_ref = []
    prop_loss_cum = []
    

    for i in range(24):
        atm_abs_ref[i] = atmos_absorption[i]
        prop_loss_indiv[i] = atmos_absorption[i] - atm_abs_ref[i] #only place I see AtmAbsRef initialized is in Reference Calc(), but is just iniliatized to AtmAbs[i]

        if (prop_loss_indiv[i] == 0):
            prop_loss_indiv[i] = -0.001
        prop_loss_cum[i] = prop_loss_cum[i] + prop_loss_indiv[i]

        
def SignalNoise():
    #does this function return anything? M2?
    prop_loss_cum = []
    detection_band = [] #Local variable to this function, just holds the calculations, and gives to M2 (which is used in Binary Search)
    target_spectrum = data_dict.given_cons['trg_name']
    max_modBackground_thresh = [] #Corresponds with 'Max Bkg/Thr' in calculated macros (ex: cells N8 - N31 in Detection Distance)
    one_third_oct_band_weight = [] #arrays from 'Data' cols AD - AH
    listener_spectrum = [] #Corresponds with green-colored Listener col K in Macros calculation (ex: K8-K31)
    Log10Div10 = 0.230258509
    TenDivLog10 = 1 / Log10Div10
    
    for i in range(24):
        if (i > 10):
            detection_band[i] = target_spectrum[i] + prop_loss_cum[i] - max_modBackground_thresh[i]
        else:
            E3 = 0
            for j in range(5):
                I0 = i + j - 2
                W4 = one_third_oct_band_weight(i, j)
                if (W4 > 0 and I0 >= 0):
                    E3 = E3 + W4 * math.exp(Log10Div10 * (target_spectrum[I0] + prop_loss_cum[I0]))
            detection_band[i] = TenDivLog10 * math.log[E3] - max_modBackground_thresh[i]
        listener_spectrum[i] = detection_band[i] + max_modBackground_thresh[i]
    M2 = detection_band[0]
    B1 = 0
    for i in range(1, 24):
        if (detection_band[i] >= M2):
            M2 = detection_band[i]
            B1 = i



def Ingard():
  Hs = data_dict.det_cons['source_height_det'] 
  Hr = data_dict.det_cons['listener_height_det'] 
  R = data_dict.given_cons['m_measure_distance']
  celsius_degrees_T1 = data_dict.meas_cons['celsius_degrees_det']
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

  for I in range(24):
    freq = ADM.freq[I] #not sure where this one comes from, will come back to it
    K1 = 2 * Pi * freq / Cs
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
    Lnfos = math.log(freq / Sigma)
    Rrc = 1 + 9.08 * math.exp(-0.75 * Lnfos)
    Xrc = 11.9 * math.exp(-0.73 * Lnfos)
    Rzr, Rzi = 0
    Done = ComplexDiv(1, 0, Rrc, Xrc, Rzr, Rzi)
    Rpr, Rpi = 0
    Done = ComplexDiv(SinTheta1 - Rzr, -Rzi, SinTheta1 + Rzr, Rzi, Rpr, Rpi)
    Srzr, Srzi = 0
    Done = ComplexMul(SinTheta1 + Rzr, Rzi, SinTheta1 + Rzr, Rzi, Srzr, Srzi)
    Czr, Czi = 0
    Done = ComplexDiv(Srzr, Srzi, 1 + SinTheta1 * Rzr, SinTheta1 * Rzi, Czr, Czi)
    Wr, Wi = 0
    Done = ComplexMul(0, 0.5 * K1 * R2, Czr, Czi, Wr, Wi)
    Wm, Wp = 0
    Done = Polar(Wr, Wi, Wm, Wp) #Means Wm = Sqr((Wr * Wr) + (Wi * Wi)) // and Wp = math.atan2(Wr, Wi)
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
            W2r, W2i = 0
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
        Ws1r, Ws1i = 0
        Done = ComplexMul(Ewrr, Ewri, Whr - 2 * Wsr, Whi - 2 * Wsi, Ws1r, Ws1i)
        Fr = 1 + Ws1r
        Fi = Ws1i
    else:
        Done = ComplexDiv(1, 0, 2 * Wr, 2 * Wi, W1r, W1i)
        Done = ComplexMul(W1r, W1i, W1r, W1i, W2r, W2i)
        Wsr = W1r + 3 * W2r
        Wsi = W1i + 3 * W2i
        W3r, W3i = 0
        Done = ComplexMul(W1r, W1i, W2r, W2i, W3r, W3i)
        Wsr = Wsr + 15 * W3r
        Wsi = Wsi + 15 * W3i
        Fr = -Wsr
        Fi = -Wsi
        Q1r, Q1i = 0
    Done = ComplexMul(Fr, Fi, 1 - Rpr, -Rpi, Q1r, Q1i)
    'Done = Polar(Fr, Fi, Fm, Fp)'
    Qr = Rpr + Q1r
    Qi = Rpi + Q1i
    Qm, Qp = 0
    Done = Polar(Qr, Qi, Qm, Qp)
    Mkdr = Mu * K1dr
    Sinc = math.sin(Mkdr) / Mkdr
    Cost = math.cos(Eta * K1dr + Qp)
    Qmr = Qm * R12
    I3 = (1 + Ea2) * (1 + Qmr * Qmr) + 2 * Qmr * (1 + Ea2 * Rhoa) * Cost * Sinc * math.exp(-Sd2 * (1 - Rhod))
    ground_effect[I] = TenDivLog10 * math.log(I3)
    
    return ground_effect #is this the right placement for this return?
    




def binary_search():
    #does this return anything?
    m_meas_distance = data_dict.given_cons['m_measure_distance']
    detection_dist = [] #not sure about this one, since it is calculated inside a macro
    precision_fraction = 0.001
    M2 = SignalNoise() #This value is calculated in Signal Noise function
    
    
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
    ground_effect = []
    atmos_absorption = []
    AtmAbsRef = []
        
    for I in range(0,23):
        ground_effect[I] = Ingard()
        atmos_absorption[I] = ansi_humidity()
        ground_effect_ref[I] = ground_effect[I]
        AtmAbsRef[I] = atmos_absorption[I]
    
    return ground_effect, atmos_absorption, ground_effect_ref, AtmAbsRef



def Reverse():
    Propagate()


    
def targetdBA():
    
    Log10Div10 = 0.230258509
    TenDivLog10 = 1 / Log10Div10
    target_spectrum = [] #Col B from Data
    A_weight_levels = [] # Col AJ from Data
    E = 0.0
    
    for i in range(24):
        E = E + math.exp(Log10Div10 * (target_spectrum[i] + A_weight_levels[i]))

    targetdBA_result = TenDivLog10 * math.log(E)

    return targetdBA_result


def BkgdBA():
    Log10Div10 = 0.230258509
    TenDivLog10 = 1 / Log10Div10
    E = 0.0
    background_noise_spec = [] #background noise from Data sheet (Col M)
    A_weight_levels = [] #A weight levels from Data sheet (Col AJ)
     
    for i in range(24):
         E = E + math.exp(Log10Div10 * (background_noise_spec[i] + A_weight_levels[i]))
         
    bkgdba_return = TenDivLog10 * math.log(E)
     
    return bkgdba_return    
 
 
def ListenerdBA():
    Log10Div10 = 0.230258509
    TenDivLog10 = 1 / Log10Div10
    E = 0.0
    target_spectrum = [] #Data Col B
    prop_loss_cum = [] #
    A_weight_levels = [] # Col AJ
    
    for i in range(24):
        E = E + math.exp(Log10Div10 * (target_spectrum[i] + prop_loss_cum[i] + A_weight_levels[i]))
    listener_dba_return = TenDivLog10 * math.log(E)
    return listener_dba_return


def calculate_measure_dist(detection_dist: float):
    return detection_dist * 25.0


def binary_search_A():
    #can you double check this function, and does it need to return anything?
    prop_loss_cum = []
    m_meas_distance = data_dict.given_cons['m_measure_distance']
    target_spectrum = data_dict.given_cons['trg_name']
    dBa = data_dict.given_cons['dba_given']
    A_weight_levels = [] # Col AJ
    precision_fraction = 0.001
    Log10Div10 = 0.230258509
    TenDivLog10 = 1 / Log10Div10
    
    
    
    z9 = -1
    detection_dist = m_meas_distance * 25
    prop_loss = prop_loss_cum(23) # <------------ Not sure about this line, not sure if it follows the logic
    Ea = 0

    while prop_loss < dBa:
        D5 = detection_dist
        detection_dist = 2 * detection_dist
        D6 = detection_dist

        for i in range(24):
            Ea += Log10Div10 * target_spectrum[i] + prop_loss_cum[i] + A_weight_levels[i]

    if z9 == 0:
        while prop_loss < dBa:
            D5 = detection_dist
            detection_dist = 2 * detection_dist
            D6 = detection_dist

            for i in range(24):
                Ea += Log10Div10 * target_spectrum[i] + prop_loss_cum[i] + A_weight_levels[i]


    while abs(D6 - D5) < precision_fraction * detection_dist:
        detection_dist = (D5 + D6) / 2
        for i in range(24):
            Ea += Log10Div10 * target_spectrum[i] + prop_loss_cum[i] + A_weight_levels[i]

        prop_loss = TenDivLog10 *  math.log(Ea, 2)

        if prop_loss > dBa:
            D5 = detection_dist
        else:
            D6 = detection_dist


def Propagate():
    InverseDistance()
    ground_effect()
    Barrier()
    Foliage()
    Atmosphere()
    SignalNoise()
    