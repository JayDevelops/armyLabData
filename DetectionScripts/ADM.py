import pandas as pd  # Install pandas in python or use Anaconda environment
import math
import data_dict
import ADM_Functions
import json

def Detection():
    data_df = pd.read_excel('./DetectionScripts/ADM - from Joel - Sept-2013.xls', sheet_name='Data')
    # # Instantiating specific columns to be pulled from the data sheet in ADM by Joel
    # Trg_name = 'Drone'  # Drone db column
    # Bkg_name = 'Urban'  # Ambient setting background noise
    # Hth_name = 'ISO Std'  # Hearing Threshold based on ISO Standard

    # Pull from data dictionary
    Trg_name = data_dict.given_cons['trg_name']
    Bkg_name = data_dict.given_cons['bkg_name'] # Ambient setting background noise
    Hth_name = data_dict.given_cons['hth_name']  # Hearing Threshold based on ISO Standard

    freq = data_df['Freq Hz'].tolist()
    target = data_df[Trg_name].tolist()
    bkg = data_df[Bkg_name].tolist()
    ht = data_df[Hth_name].tolist()
    B3 = data_df.iloc[:24, 34].tolist()
    awt_weights = data_df['Awt weights'].tolist()
    ai_weights = data_df['A.I. weights'].tolist()
    third_obands_df = data_df.iloc[:, 29:34].copy()
    third_obands = []
    for i in range(5):
        third_obands.append(third_obands_df.iloc[:24,i].tolist())

    # print(third_obands)
    # print(len(third_obands))
    # print(freq)
    # print(target)
    # print(bkg)
    # print(ht)
    # print(awt_weights)
    # print(ai_weights)
    # print(third_octave_bands)

    #Reference Calc Call
    ground_effect, atmos_absorption, ground_effect_ref, AtmAbsRef = ADM_Functions.reference_calc(freq)
    #Call Dprime
    max_modBkg_thresh = ADM_Functions.Dprime(bkg, ht, third_obands, B3)
    #Binary Search Call
    detection_distance = ADM_Functions.binary_search(freq, ground_effect_ref, atmos_absorption, AtmAbsRef, target, max_modBkg_thresh, third_obands)
    
    return detection_distance

# Main Function Declaration and Call
def main():
    # Possibly implement a switch case to consider user input
    # of which Macro to call
    detection_distance = Detection()
    detection_data = {}
    detection_data['detection_distance'] = detection_distance
    with open('detection_data.json', 'w') as f:
        json.dump(detection_data, f)

    #print(detection_distance)


if __name__ == '__main__':
    main()
