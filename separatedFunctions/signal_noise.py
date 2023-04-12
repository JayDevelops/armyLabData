import pandas as pd  # Install pandas in python or use Anaconda environment
import math
from DetectionScripts import data_dict
from DetectionScripts import ADM_Functions

def signal_noise(target_spectrum, max_modBkg, third_oband, prop_loss_cum):
    detection_band = [0] * 24
    listener_spectrum = [0] * 24
    Log10Div10 = math.log(10) / 10
    TenDivLog10 = 10 / math.log(10)
    print('detection_band', detection_band)

    for i in range(24):
        if i > 10:
            detection_band[i] = target_spectrum[i] + prop_loss_cum[i] - max_modBkg[i]
            #print(target_spectrum)
            #print(prop_loss_cum)
            #print(max_modBkg)
            #print('detection_band', detection_band)
        else:
            E3 = 0
            for j in range(5):
                I0 = i + j - 2
                W4 = third_oband[j][i]
                if (W4 > 0 and I0 >= 0):
                    E3 = E3 + W4 * math.exp(Log10Div10 * (target_spectrum[I0] + prop_loss_cum[I0]))
                    #print(target_spectrum)
                    #print(prop_loss_cum)
                    #print(E3)
            detection_band[i] = TenDivLog10 * math.log(E3) - max_modBkg[i]
            #print('detection_band', detection_band)
        listener_spectrum[i] = detection_band[i] + max_modBkg[i]
    M2 = detection_band[0]
    print('detection_band', detection_band[0])

    B1 = 0
    for i in range(24):
        if (detection_band[i] >= M2):
            M2 = detection_band[i]
            B1 = i
    print('signal_noise', prop_loss_cum)
    return M2, B1, prop_loss_cum