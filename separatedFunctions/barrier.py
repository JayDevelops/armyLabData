import pandas as pd  # Install pandas in python or use Anaconda environment
import math
from DetectionScripts import data_dict
from DetectionScripts import ADM_Functions

def Barrier(freq, prop_loss_cum, prop_loss_indiv, detection_dist):
    barr_atten = []
    barr_num = data_dict.det_cons['barrier_on']
    distance_from_source = data_dict.det_cons['barrier_dist']
    barrier_height_det = data_dict.det_cons['barrier_dist']
    source_height_det = data_dict.det_cons['source_height_det']
    listener_height_det = data_dict.det_cons['listener_height_det']
    celsius_degrees_T1 = data_dict.det_cons['celsius_degrees_det']
    Pi = 4 * math.atan(1)
    Log10Div10 = 0.230258509
    TenDivLog10 = 1 / Log10Div10
    
    for i in range(24):
        barr_atten.append(-0.001)
    if barr_num > 0:
        if detection_dist >= distance_from_source:
            Hba = barrier_height_det - (source_height_det + distance_from_source * (listener_height_det) / detection_dist)
            X1 = math.sqrt(distance_from_source ** 2 + ((barrier_height_det - source_height_det) **2))
            X1 = X1 + math.sqrt(((detection_dist - distance_from_source) ** 2) + ((barrier_height_det - listener_height_det) ** 2))
            X1 = X1 - math.sqrt((detection_dist **2) + ((listener_height_det - source_height_det)**2))
            Cs = 331.4 * math.sqrt(1 + celsius_degrees_T1 / 273.15)
            
            for i in range(24):
                if(Hba < 0):
                    X2 = 2 * X1 * freq[2] / Cs #not sure about the freq part
                    X3 = math.sqrt(2 * Pi * abs(X2))
                    X2 = -(X2)
                    if(X2 <= -0.1916):
                        barr_atten[i] = -0.01
                    else: 
                        barr_atten[i] = -(5 + 2 * TenDivLog10 * math.log(X3 / math.tan(X3))) - 0.01
                else:
                    X2 = 2 * X1 * freq[i] / Cs #not sure about the freq
                    X3 = math.sqrt(2 * Pi * abs(X2))
                    X2 = math.exp(2 * X3)
                    if(X2 == 1):
                        barr_atten[i] = -5 -0.01
                    else:
                        barr_atten[i] = -(5 + 2 * TenDivLog10 * math.log(X3 * (X2 + 1) / (X2 - 1))) - 0.01
                        if(barr_atten[i] < -20):
                            barr_atten[i] = -20
    
    for i in range(24):
        prop_loss_indiv[i] = barr_atten[i]
        prop_loss_cum[i] = prop_loss_cum[i] + prop_loss_indiv[i]
    print('barrier', prop_loss_cum)
    return prop_loss_cum, prop_loss_indiv
