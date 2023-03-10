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
        ground_effect_ref(I) = ground_effect(I)
        AtmAbsRef(I) = atmos_absorption(I)




def calculate_measure_dist(detection_dist: float):
    return detection_dist * 25.0



# Main Function Declaration and Call
def main():
    # Possibly implement a switch case to consider user input
    # of which Macro to call
    Detection()


if __name__ == '__main__':
    main()