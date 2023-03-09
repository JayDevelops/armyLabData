import math

import pandas as pd


# Reads "Data" sheet from ADM by Joel for reference values and returns it
def read_data():
    adm_data_df = pd.read_excel('ADM - from Joel - Sept-2013.xls', sheet_name='Data')
    # adm_data_df = pd.read_excel('adm.xls', sheet_name='Data')
    return adm_data_df


# Extract relevant reference values to be used in calcuating detectability
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


def Ingard():
    return 2

def AnsiHumidity():
    return 7

def inverse_distance(measure_dist: float, detection_dist: float) -> list:
    # 10 divided by log 10 static, and log of measure distance divided by detection distance log in base 10
    ten_divided_by_log_10, log_m_dist_by_d_dist = (10 / math.log(10, 10)), math.log(measure_dist / detection_dist, 10)

    inverse_distances = []

    for x in range(0, 23):
        inv_dist = 2 * ten_divided_by_log_10 * log_m_dist_by_d_dist
        if inv_dist == 0:
            inverse_distances.append(inv_dist)
        else:
            inverse_distances.append(-0.001)

    return inverse_distances

"""
ground_effect_reference(26) is 'Reference Ground Effect during measurement' from .vbs file
"""
def reference_calc(sigma: float, source_height_condition: float,
                   measure_dist: float, mic_height_coord: float,
                   em2_ref: float, ground_effect_ref: float) -> list:
    arr = []
    for x in range(0, 23):
        new_num = Ingard() + AnsiHumidity()
        arr.append(Ingard() + AnsiHumidity())

    return arr



def calculate_measure_dist(detection_dist: float):
    return detection_dist * 25.0

def targetdBA(E: double):
    ten_divided_by_log_10, log_10_div_10 = (10 / math.log(10, 10)), (math.log(10,10) / 10)
    E = 0
    for x in range(0, 23):
        E = E + math.exp(log_10_div_10 * s1[x], s10[x]) 
    
    targetdBA = ten_divided_by_log_10 * math.log(E, 10)

    return targetdBA