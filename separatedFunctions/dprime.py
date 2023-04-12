import pandas as pd  # Install pandas in python or use Anaconda environment
import math
from DetectionScripts import data_dict
from DetectionScripts import ADM_Functions

def Dprime(bkg, ht, third_obands, B3):
    """
    Functions called in Dprime -> NormalDeviate, mod_background_noise, background_noise_spec, B3,
    one_third_oct_band_weight, max_modBackground_thresh, hearing_threshold, mod_background_noise_2
    
    ***The rest in parameter list are variables
    """
    hit_prob = data_dict.lis_cons['hit_prob']
    false_alarm_rate = data_dict.lis_cons['false_alarm_rate']
    mod_background_noise = [0] * 24
    observer_efficiency = data_dict.lis_cons['observer_efficiency']
    max_modBackground_thresh = [0] * 24
    TenDivLog10 = 10 / math.log(10)
    Log10Div10 = 1/TenDivLog10
    mic_height = data_dict.meas_cons['mic_height_meas']

    if(hit_prob > 0.5):
        Z3 = normal_deviate(1 - hit_prob)
    else:
        Z3 = normal_deviate(hit_prob)
    
    if(false_alarm_rate > 0.5):
        Z = normal_deviate(1 - false_alarm_rate)
    else:
        Z = normal_deviate(false_alarm_rate)
        
    calculate_d = Z3 - Z
    if(hit_prob > 0.5):
        if(false_alarm_rate > 0.5):
            calculate_d = Z3 - Z
        else:
            calculate_d = Z3 + Z
    else:
        calculate_d = Z - Z3
        
    data_dict.lis_cons.update({'d_stat': calculate_d})

    for I in range(24):
        if I > 10:
            mod_background_noise[I] = bkg[I] + TenDivLog10 * math.log(calculate_d / observer_efficiency / B3[I])
        else:
            E4 = 0
            for J in range(5):
                I0 = I + J - 2
                W4 = third_obands[J][I]
                if((W4 > 0) and (I0 >= 0)):
                    E4 += W4 * math.exp(Log10Div10 * bkg[I0])
            mod_background_noise[I] = TenDivLog10 * math.log(E4 * calculate_d / observer_efficiency / B3[I])
            
        if ht[I] > mod_background_noise[I]:
            max_modBackground_thresh[I] = ht[I]
        else:
            max_modBackground_thresh[I] = mod_background_noise[I]

    return max_modBackground_thresh