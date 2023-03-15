import pandas as pd  # Install pandas in python or use Anaconda environment
import data  # python module including helper functions (In our case, translated Macros from ADM by Jeol)
import math

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



# Main Function Declaration and Call
def main():
    # Possibly implement a switch case to consider user input
    # of which Macro to call
    Detection()


