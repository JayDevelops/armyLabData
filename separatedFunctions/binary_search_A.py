import pandas as pd  # Install pandas in python or use Anaconda environment
import math
from DetectionScripts import data_dict
from DetectionScripts import ADM_Functions

def binary_search_A(dBa, m_meas_distance, precision_fraction, D6, D5):
    z9 = -1
    detection_dist = m_meas_distance * 25
    prop_loss = prop_loss_cum(23)
    Ea = 0

    while prop_loss < dBa:
        D5 = detection_dist
        detection_dist = 2 * detection_dist
        D6 = detection_dist

        for x in range(0, 23):
            Ea += Log10Div10 * target_spectrum(x) + prop_loss_cum(x) + A_weight_levels(x)

    if z9 == 0:
        while prop_loss < dBa:
            D5 = detection_dist
            detection_dist = 2 * detection_dist
            D6 = detection_dist

            for x in range(0, 24):
                Ea += Log10Div10 * target_spectrum(x) + prop_loss_cum(x) + A_weight_levels(x)


    while abs(D6 - D5) < precision_fraction * detection_dist:
        detection_dist = (D5 + D6) / 2
        for x in range(0, 24):
            Ea += Log10Div10 * target_spectrum(x) + prop_loss_cum(x) + A_weight_levels(x)

        prop_loss = TenDivLog10 *  log(Ea, 2)

        if prop_loss > dBa:
            D5 = detection_dist
        else:
            D6 = detection_dist
