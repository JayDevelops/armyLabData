import pandas as pd  # Install pandas in python or use Anaconda environment
import math
from DetectionScripts import data_dict
from DetectionScripts import ADM_Functions

def atmosphere(freq, prop_loss_cum, prop_loss_indiv, atm_abs, atm_abs_ref):
    ansi_humidity(freq)

    for i in range(24):
        prop_loss_indiv[i] = atm_abs[i] - atm_abs_ref[i] #only place I see AtmAbsRef initialized is in Reference Calc(), but is just iniliatized to AtmAbs[i]

        if (prop_loss_indiv[i] == 0):
            prop_loss_indiv[i] = -0.001
        prop_loss_cum[i] = prop_loss_cum[i] + prop_loss_indiv[i]
    print('atmosphere', prop_loss_cum)
    return prop_loss_cum, prop_loss_indiv