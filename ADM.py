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


def calculate_measure_dist(detection_dist: float):
    return detection_dist * 25.0


def test_inv_dist():
    measure_dist = 30.0
    detection_dist = calculate_measure_dist(30.0)

    inv_distances = inverse_distance(measure_dist, detection_dist)

    for x in inv_distances:
        print(x)


test_inv_dist()

import pandas as pd

modelSheet_df = pd.read_excel('admDataSet.xls', sheet_name="Model", usecols='A:V')
dataSheet_df = pd.read_excel('admDataSet.xls', sheet_name="Data")

excel_file = 'admDataSet.xls'

def InitMacros():
    # Rg1 = Range()
    # Delete all charts
    # DeleteAllCharts()
    S1 = []
    S4 = []
    S5 = []
    S6 = []
    W3 = [[], []]
    B3 = []

    # Read target, background, and health indicators from Model sheet
    # Trg = modelSheet_df.Range("B6").Value
    Trg = pd.read_excel(excel_file, sheet_name="Model", usecols='B', nrows=5)
    Bkg = pd.read_excel(excel_file, sheet_name="Model", usecols='C', nrows=5)
    Hth = pd.read_excel(excel_file, sheet_name="Model", usecols='D', nrows=5)

    # Load S4 values from Data sheet to Model sheet
    for i in range(24):
        S4[i] = dataSheet_df.Range("A2").Offset(i, 0).Value
        modelSheet_df.Range("A8").Offset(i, 0).Value = S4[i]

    if Trg != 0:
        # Load S1 values from Data sheet to Model sheet
        for i in range(24):
            S1[i] = dataSheet_df.Range("A2").Offset(i, Trg.iloc[-1.0]).Value
            modelSheet_df.Range("A8").Offset(i, 1).Value = S1[i]

        # Load D3 value from Data sheet to Model sheet
        D3 = dataSheet_df.Range("A2").Offset(24, Trg).Value
        modelSheet_df.Range("C3").Value = D3
    else:
        # Check if a distance value exists in Cell C3 for Microphone from Source
        Rg1 = modelSheet_df.Range("C3")
        if isinstance(Rg1.value, (int, float)) and len(str(Rg1.value)) > 0:
            D3 = Rg1.value
        else:
            # MsgBox("Enter a Distance in meters in Cell C3 for Microphone from Source")
            return

    if Bkg != 0:
        # Load S6 values from Data sheet to Model sheet
        for i in range(24):
            S6[i] = dataSheet_df.Range("A2").Offset(i, Bkg).Value
            modelSheet_df.Range("A8").Offset(i, 2).Value = S6[i]

    if Hth != 0:
        # Load S5 values from Data sheet to Model sheet
        for i in range(24):
            S5[i] = dataSheet_df.Range("A2").Offset(i, Hth).Value
            modelSheet_df.Range("A8").Offset(i, 3).Value = S5[i]

    # Load S10 and AiWt values from Data sheet
    for i in range(0, 24):
        S10 = [dataSheet_df.Range("A2").Offset(i, 35).Value for i in range(24)]
    for i in range(0, 24):
        AiWt = [dataSheet_df.Range("A2").Offset(i, 37).Value for i in range(24)]

    ColOffset = 29
    for i in range(0, 24):
        # Load W3 values from Data sheet
        for j in range(5):
            W3[i][j] = dataSheet_df.Range("A2").Offset(i, j + ColOffset).Value
        # Load B3 values from Data sheet
        j = 5
        B3[i] = dataSheet_df.Range("A2").Offset(i, j + ColOffset).Value

    # Set some constants
    Pi = 4 * math.atan(1)
    Log10Div10 = 0.230258509
    TenDivLog10 = 1 / Log10Div10
    Log10 = 2.302585093
    Co = 331.32
    Too = 273.15
    A1 = 0.001
    H2ref = modelSheet_df.Range("A3").Value  # m Source Height meas
    H3ref = modelSheet_df.Range("B3").Value  # m Mic Height meas
    # D3 m in "C3"  Mic distance from source
    Tref = modelSheet_df.Range("D3").Value  # Deg C meas
    Href = modelSheet_df.Range("E3").Value  # % r.h. meas
    SigmaRef = modelSheet_df.Range("F3").Value  # m Flow resistivity meas
    Em2Ref = modelSheet_df.Range("G3").Value  # Em2 turbulence factor meas
    T1 = Tref
    H1 = Href

    # ReferenceCalc call
    H2 = modelSheet_df.Range("H3").Value                  # m Source Height det
    H3 = modelSheet_df.Range("I3").Value                  # m Listener Height det
    T1 = modelSheet_df.Range("J3").Value                  # Deg C det
    H1 = modelSheet_df.Range("K3").Value                  # % r.h. det
    SigmaDet = modelSheet_df.Range("L3").Value
    Em2Det = modelSheet_df.Range("M3").Value
    WindSpeed = modelSheet_df.Range("N3").Value
    E1 = modelSheet_df.Range("r3").Value
    P1 = modelSheet_df.Range("s3").Value
    P2 = modelSheet_df.Range("t3").Value

    # Assigning D1 to the function Dprime
    # D1 = Dprime
    # modelSheet_df.Range("u3").Value = D1

    WindFlag = 0
    WindDir = " Upwind"
    b9 = modelSheet_df.Range("H5").Value
    b7 = modelSheet_df.Range("I5").Value
    b8 = modelSheet_df.Range("J5").Value
    N1 = modelSheet_df.Range("K5").Value
    W1 = modelSheet_df.Range("L5").Value
    W2 = modelSheet_df.Range("M5").Value
    Fl = modelSheet_df.Range("N5").Value
    Al = modelSheet_df.Range("O5").Value
    F7 = 1
    Iwthr1 = 0
    Surface = "Grass"
    Bnumber = 11
    Tnumber = 3
    Hnumber = 2


import math
import pandas as pd

N1 = pd.read_excel('admDataSet.xls', sheet_name='Model', usecols='K', nrows=5)
D4 = pd.read_excel('admDataSet.xls', sheet_name='Model', usecols='C', nrows=5)
W1 = pd.read_excel('admDataSet.xls', sheet_name='Model', usecols='L', nrows=5)
W2 = pd.read_excel('admDataSet.xls', sheet_name='Model', usecols='M', nrows=5)
Fl = pd.read_excel('admDataSet.xls', sheet_name='Model', usecols='N', nrows=5)
Al = pd.read_excel('admDataSet.xls', sheet_name='Model', usecols='O', nrows=5)
N1 = N1.iloc[-2, 0]
D4 = D4.iloc[-3, 0]
W1 = W1.iloc[-2, 0]
W2 = W2.iloc[-2, 0]
Fl = Fl.iloc[-2, 0]
Al = Al.iloc[-2, 0]

S4 = [50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000,
      6300, 8000, 10000]

Cs = 340.29
S8 = [0] * 24
S3 = [0] * 24


def Foliage(N1, D4, W1, W2, Fl, Cs, Al, S4, S8, S3):
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


Fo_values = Foliage(N1, D4, W1, W2, Fl, Cs, Al, S4, S8, S3)
print(Fo_values)