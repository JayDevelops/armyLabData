import pandas as pd #Install pandas in python or use Anaconda environment
import data #python module including helper functions (In our case, translated Macros from ADM by Jeol)

def Detection():
    # data_df = data.read_data()
    # print(data_df.iloc[:,29:34])
    
    data_df = pd.read_excel('ADM - from Joel - Sept-2013.xls', sheet_name='Data')
    
    #Instantiating specific columns to be pulled from the data sheet in ADM by Joel
    Trg_name = 'Drone' #Drone db column 
    Bkg_name = 'Urban' #Ambient setting background noise
    Hth_name = 'ISO Std' #Hearing Threshold based on ISO Standard

    freq_df = data_df['Freq Hz']
    target_df = data_df[Trg_name]
    bkg_df = data_df[Bkg_name]
    ht_df = data_df[Hth_name]
    awt_weights_df = data_df['Awt weights']
    ai_weights_df = data_df['A.I. weights']
    third_octave_bands_df = data_df.iloc[:,29:34].copy()

    # #Find and set all necessary dataframes for further processing
    # freq_df, target_df, bkg_noise_df, \
    #     hear_thresh_df, awt_weights, ai_weights_df, \
    #     third_obands_df = data.get_data(Trg_name, Bkg_name, Hth_name, data_df)
    
    # print(third_obands_df)
    

#Main Function Declaration and Call
def main():
    # Possibly implement a switch case to consider user input 
    # of which Macro to call
    Detection()

if __name__=='__main__':
    main()

