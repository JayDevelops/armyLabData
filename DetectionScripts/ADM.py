import pandas as pd  # Install pandas in python or use Anaconda environment
import math
import data_dict

def Detection():
    # data_df = data.read_data()
    # print(data_df.iloc[:,29:34])

    data_df = pd.read_excel('ADM - from Joel - Sept-2013.xls', sheet_name='Data')

    # # Instantiating specific columns to be pulled from the data sheet in ADM by Joel
    # Trg_name = 'Drone'  # Drone db column
    # Bkg_name = 'Urban'  # Ambient setting background noise
    # Hth_name = 'ISO Std'  # Hearing Threshold based on ISO Standard

    # Pull from data dictionary
    Trg_name = data_dict.given_cons['trg_name']
    Bkg_name = data_dict.given_cons['bkg_name'] # Ambient setting background noise
    Hth_name = data_dict.given_cons['hth_name']  # Hearing Threshold based on ISO Standard
    

    # freq = data_df['Freq Hz'].tolist()
    # target = data_df[Trg_name].tolist()
    # bkg = data_df[Bkg_name].tolist()
    # ht = data_df[Hth_name].tolist()
    # awt_weights = data_df['Awt weights'].tolist()
    # ai_weights = data_df['A.I. weights'].tolist()
    # third_octave_bands_df = data_df.iloc[:, 29:34].copy()

    freq_df = data_df['Freq Hz']
    target_df = data_df[Trg_name]
    bkg_df = data_df[Bkg_name]
    ht_df = data_df[Hth_name]
    awt_weights_df = data_df['Awt weights']
    ai_weights_df = data_df['A.I. weights']
    third_octave_bands_df = data_df.iloc[:, 29:34].copy()

    # print(freq)
    # print(target)
    # print(bkg)
    # print(ht)
    # print(awt_weights)
    # print(ai_weights)

# Main Function Declaration and Call
def main():
    # Possibly implement a switch case to consider user input
    # of which Macro to call
    Detection()


if __name__ == '__main__':
    main()
