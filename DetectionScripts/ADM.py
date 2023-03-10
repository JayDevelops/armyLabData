import pandas as pd  # Install pandas in python or use Anaconda environment
import data  # python module including helper functions (In our case, translated Macros from ADM by Jeol)
import math

# global variables
prop_loss_cum = 0

def Detection():
    # data_df = data.read_data()
    # print(data_df.iloc[:,29:34])


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
    return 2

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

        
def normal_deviate(p):
    t = math.sqrt(-2 * math.log(p))
    return t - (2.30753 + t * 0.27061) / (1 + t * (0.99229 + t * 0.04481))

def ground_effect():
    TrgHgt = source_height_det
    DetHgt = listener_height_det
    R = detection_dist
    Sigma = SigmaDet
    Em2 = Em2Det
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
                
                
#NEEDS WORK <------------------------------------------------
def binary_searchA(m_meas_distance, D5, D6, M2, precision_fraction):
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

                
#NEEDS WORK <------------------------------------------------                
def ingard():
    print("Needs Work")

#NEEDS WORK <------------------------------------------------     
def dprime():
    print("Needs Work")
    
    
"""
ground_effect_reference(26) is 'Reference Ground Effect during measurement' from .vbs file
"""
def reference_calc(ground_effect_ref):
    for I in range(0,23):
        ground_effect_ref(I) = ground_effect(I)
        AtmAbsRef(I) = atmos_absorption(I)


def Foliage(N1, D4, W1, W2, Fl, Cs, Al, S4, S8, S3):
    excel_file = "ADM - from Joel - Sept-2013.xls"
    modelSheet_df = pd.read_excel("ADM - from Joel - Sept-2013.xls", sheet_name="Model", usecols='A:V')
    dataSheet_df = pd.read_excel("ADM - from Joel - Sept-2013.xls", sheet_name="Data")
    
    N1 = pd.read_excel(excel_file, sheet_name='Model', usecols='K', nrows=5)
    D4 = pd.read_excel(excel_file, sheet_name='Model', usecols='C', nrows=5)
    W1 = pd.read_excel(excel_file, sheet_name='Model', usecols='L', nrows=5)
    W2 = pd.read_excel(excel_file, sheet_name='Model', usecols='M', nrows=5)
    Fl = pd.read_excel(excel_file, sheet_name='Model', usecols='N', nrows=5)
    Al = pd.read_excel(excel_file, sheet_name='Model', usecols='O', nrows=5)
    N1 = N1.iloc[-2, 0]
    D4 = D4.iloc[-3, 0]
    W1 = W1.iloc[-2, 0]
    W2 = W2.iloc[-2, 0]
    Fl = Fl.iloc[-2, 0]
    Al = Al.iloc[-2, 0]
    S4 = dataSheet_df.iloc[:-3,0].values.tolist()
    Cs = 340.29
    S8 = [0] * 24
    S3 = [0] * 24
    Fo = [-0.001] * 24
    Fo_list = []

    if N1 > 0:
        # If detection distance is greater than the distance from source to edge of foliage
        if D4 > W1:
            # Sets the distance difference to X2
            X2 = D4 - W1
            # If the new detection distance is greater than depth (extent) of foliage in meters
            if X2 > W2:
                # Set distance the depth (extent) of foliage in meters
                X2 = W2
            X2 = math.sqrt(X2)
            Cons = 2.647 / math.log(10)

            for I in range(10, 24):
                Ka = (2 * math.pi * S4[I] / Cs) * Al / 100
                if Ka < 0.401:
                    Fo[I] = -0.01
                elif Ka < 5:
                    Fo[I] = -X2 * math.sqrt(Fl) * (Cons * math.log(Ka) + 1.05)
                else:
                    Fo[I] = -X2 * math.sqrt(Fl) * 2.9
                Fo_list.append(Fo)
    else:
        for I in range(0, 24):
            S8[I] = Fo[I]
            S3[I] = S3[I] + S8[I]
            Fo_list.append(S3)
    return Fo_list

def initMacros():
    excel_file = "ADM - from Joel - Sept-2013.xls"
    modelSheet_df = pd.read_excel("ADM - from Joel - Sept-2013.xls", sheet_name="Model", usecols='A:V')
    dataSheet_df = pd.read_excel("ADM - from Joel - Sept-2013.xls", sheet_name="Data")
    Trg = pd.read_excel(excel_file,sheet_name="Model",usecols='B',nrows=5)
    Bkg = pd.read_excel(excel_file,sheet_name="Model",usecols='C',nrows=5)
    Hth = pd.read_excel(excel_file,sheet_name="Model",usecols='D',nrows=5)

    # Extract the frequency values from the first column of the DataFrame
    S4 = dataSheet_df.iloc[:-3,0].values.tolist()
    
    # Loading Frequency values from data to model sheet
    for i in range(len(S4)):
        # Get the value from the first column of the current row in "Data" sheet
        freqValues = S4[i]
        
        # Assign the value to S4 and the corresponding cell in "Model" sheet
        modelSheet_df.loc[7+i, 'A'] = freqValues

    if Trg.iloc[-1,0] != 0:
        # Target Spectrum
        S1 = dataSheet_df.iloc[:-3,1].values.tolist()
        
        for i in range(len(S1)):
            targetValues = S1[i]
            modelSheet_df.iloc[i+7, 1] = targetValues

        # Load D3 Value from Data sheet to Model sheet
        D3 = dataSheet_df.iloc[24,Trg.iloc[-1].astype(int)].values
        modelSheet_df.iloc[2,2] = D3

    # else:
    #     # Check if a distance value exists in Cell C3 for Microphone from Source
    #     Rg1 = modelSheet_df.iloc[2, 2]
    #
    #     if isinstance(Rg1, (int, float)) and len(str(Rg1)) > 0:
    #         D3 = Rg1
    #
    #     else:
    #         # MsgBox("Enter a Distance in meters in Cell C3 for Microphone from Source")
    #         return

    if Bkg.iloc[-1,0] != 0:
        S6 = dataSheet_df.iloc[:-3,Bkg.iloc[-1,0]].values.tolist()
        for i in range(len(S6)):
            bkgNoiseValues = S6[i]
            modelSheet_df.iloc[i+7, 2] = bkgNoiseValues

    if Hth.iloc[-1,0] != 0:
        s5 = dataSheet_df.iloc[:-3,Hth.iloc[-1,0]].values.tolist()
        for i in range(len(s5)):
            hearingValues = s5[i]
            modelSheet_df.iloc[i+7, 3] = hearingValues

    # -----------------------------------------------------------------------------------------------
    # Loading Ai weights
    s5 = dataSheet_df.iloc[:-3,35].values.tolist()
    AiWt = dataSheet_df.iloc[:-3,37].values.tolist()

    Log10Div10 = 0.230258509
    TenDivLog10 = 1 / Log10Div10
    Log10 = 2.302585093
    Co = 331.32
    Too = 273.15
    A1 = 0.001

    # m Source Height meas
    h2ref = pd.read_excel(excel_file, sheet_name='Model', usecols='A', nrows=3)
    h2ref = h2ref.iloc[-2,0]

    # m Mic Height meas
    h3ref = pd.read_excel(excel_file, sheet_name='Model', usecols='B', nrows=3)
    h3ref = h3ref.iloc[-2,0]

    # Deg C meas
    tRef = pd.read_excel(excel_file, sheet_name='Model', usecols='D', nrows=3)
    tRef = tRef.iloc[-2,0]

    # % r.h. meas
    hRef = pd.read_excel(excel_file, sheet_name='Model', usecols='E', nrows=3)
    hRef = hRef.iloc[-2,0]

    # m Flow resistivity meas
    sigmaRef = pd.read_excel(excel_file, sheet_name='Model', usecols='F', nrows=3)
    sigmaRef = sigmaRef.iloc[-2,0]

    # Em2 turbulence factor meas
    em2Ref = pd.read_excel(excel_file, sheet_name='Model', usecols='G', nrows=3)
    em2Ref = em2Ref.iloc[-2,0]

    t1 = tRef
    h1 = hRef

    # Reference Calc call
    reference_calc()

    # m Source Height det
    h2 = pd.read_excel(excel_file, sheet_name='Model', usecols='H', nrows=3)
    h2 = h2.iloc[-2,0]

    # m Listener Height det
    h3 = pd.read_excel(excel_file, sheet_name='Model', usecols='I', nrows=3)
    h3 = h3.iloc[-2,0]

    # Deg C det
    t1 = pd.read_excel(excel_file, sheet_name='Model', usecols='J', nrows=3)
    t1 = t1.iloc[-2,0]

    # % r.h. det
    h1 = pd.read_excel(excel_file, sheet_name='Model', usecols='K', nrows=3)
    h1 = h1.iloc[-2,0]

    # m Flow resistivity det
    sigmaDelt = pd.read_excel(excel_file, sheet_name='Model', usecols='L', nrows=3)
    sigmaDelt = sigmaDelt.iloc[-2,0]

    # Em2 turbulence factor det
    em2Det = pd.read_excel(excel_file, sheet_name='Model', usecols='M', nrows=3)
    em2Det = em2Det.iloc[-2,0]

    # Wind speed det
    windSpeed = pd.read_excel(excel_file, sheet_name='Model', usecols='N', nrows=3)
    windSpeed = windSpeed.iloc[-2,0]

    # Observer efficiency
    e1 = pd.read_excel(excel_file, sheet_name='Model', usecols='R', nrows=3)
    e1 = e1.iloc[-2,0]

    # Hit prob
    p1 = pd.read_excel(excel_file, sheet_name='Model', usecols='S', nrows=3)
    p1 = p1.iloc[-2,0]

    # False alarm prop
    p2 = pd.read_excel(excel_file, sheet_name='Model', usecols='T', nrows=3)
    p2 = p2.iloc[-2,0]

    # Calculate d' statistic
    D1 = pd.read_excel(excel_file, sheet_name='Model', usecols='U', nrows=3)
    D1 = D1.iloc[-2, 0]

    windFlag = 0
    windDir = "Upwind"

    # barrier? 0 or 1
    b9 = pd.read_excel(excel_file, sheet_name='Model', usecols='H', nrows=5)

    # distance from source m
    b7 = pd.read_excel(excel_file, sheet_name='Model', usecols='I', nrows=5)

    # height m
    b8 = pd.read_excel(excel_file, sheet_name='Model', usecols='J', nrows=5)

    # foliage? 0 or 1
    N1 = pd.read_excel(excel_file, sheet_name='Model', usecols='K', nrows=5)
    N1 = N1.iloc[-2, 0]

    # distance in meters from source to near edge of foliage
    W1 = pd.read_excel(excel_file, sheet_name='Model', usecols='L', nrows=5)
    W1 = W1.iloc[-2, 0]

    # depth (extent) of foliage in meters
    W2 = pd.read_excel(excel_file, sheet_name='Model', usecols='M', nrows=5)
    W2 = W2.iloc[-2, 0]

    # leaf area per unit vol dense hardwood brush in m^-1
    Fl = pd.read_excel(excel_file, sheet_name='Model', usecols='N', nrows=5)
    Fl = Fl.iloc[-2, 0]

    # average leaf width in cm
    Al = pd.read_excel(excel_file, sheet_name='Model', usecols='O', nrows=5)
    Al = Al.iloc[-2, 0]

    # Type of surface
    F7 = 1
    Iwthr1 = 0
    Surface = "Grass"

    # Default background number 12=G.C. 11=Low EPA
    Bnumber = 11

    # Typical vehicle
    Tnumber = 3

    # ISO Hearing Threshold for Pure tones
    Hnumber = 2

def targetdBA():
    log_10_div_10 = 0.230258509
    ten_divided_by_log_10 = 1 / log_10_div_10

    E = float()
    E = 0.0
    for x in range(23):
        E = E + math.exp(log_10_div_10 * (S1[x] + S10[x]))

    targetdBA_result = ten_divided_by_log_10 * math.log(E)

    return targetdBA_result

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

# adds the previous E value in the function with the exp function of Log10Div10 * 3 separate array values
def ListenerdBA():
    E = float()
    E = 0.0
    for I in range(0, 23):
        E = E + math.exp(Log10Div10 * (S1[I] + S3[I] + S10[I]))
    listener_dba_return = TenDivLog10 * math.log(E)
    return listener_dba_return

def calculate_measure_dist(detection_dist: float):
    return detection_dist * 25.0

def targetdBA():
    log_10_div_10 = 0.230258509
    ten_divided_by_log_10 = 1 / log_10_div_10

    E = float()
    E = 0.0
    for x in range(23):
        E = E + math.exp(log_10_div_10 * (S1[x] + S10[x])) 
    
    targetdBA_result = ten_divided_by_log_10 * math.log(E)

    return targetdBA_result


# Main Function Declaration and Call
def main():
    # Possibly implement a switch case to consider user input
    # of which Macro to call
    Detection()


if __name__ == '__main__':
    main()
