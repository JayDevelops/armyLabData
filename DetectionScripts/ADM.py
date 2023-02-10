import pandas as pd #Install pandas in python or use Anaconda environment
import data #python module including helper functions (In our case, translated Macros from ADM by Jeol)

def Detection():
    data_df = data.read_data()
    print(data_df.iloc[:,29:34])

    #Instantiating specific columns to be pulled from the data sheet in ADM by Joel
    Trg_name = 'Drone' #Drone db column 
    Bkg_name = 'Urban' #Ambient setting background noise
    Hth_name = 'ISO Std' #Hearing Threshold based on ISO Standard

    data.get_data(Trg_name, Bkg_name, Hth_name, data_df)

#Main Function Declaration and Invocation
def main():
    Detection()

if __name__=='__main__':
    main()

