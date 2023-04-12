import pandas as pd  # Install pandas in python or use Anaconda environment
import math
from DetectionScripts import data_dict
from DetectionScripts import ADM_Functions

def binary_search(freq, ground_effect_ref, atm_abs, atm_abs_ref, target, max_modBkg, third_oband):
    Z9 = -1
    m_measure_distance = data_dict.given_cons['m_measure_distance']
    detection_dist = m_measure_distance * 25
    precision_fraction = 0.001
    prop_loss_cum = [0] * 24
    prop_loss_indiv = [0] * 24
    Z9 = Z9 + 1
    M2, B1, prop_loss_cum, prop_loss_indiv = Propagate(freq, detection_dist, ground_effect_ref, atm_abs, atm_abs_ref, target, max_modBkg, third_oband, prop_loss_cum, prop_loss_indiv)

    while M2 >= 0:
        Z9 = Z9 + 1
        D5 = detection_dist
        detection_dist = 2 * detection_dist
        D6 = detection_dist
        M2, B1, prop_loss_cum, prop_loss_indiv = Propagate(freq, detection_dist, ground_effect_ref, atm_abs, atm_abs_ref, target, max_modBkg, third_oband, prop_loss_cum, prop_loss_indiv)
        #print(M2)

    if Z9 == 0:
        while M2 <= 0:
            D6 = detection_dist
            detection_dist = detection_dist / 2
            D5 = detection_dist
            M2, B1, prop_loss_cum, prop_loss_indiv = Propagate(freq, detection_dist, ground_effect_ref, atm_abs, atm_abs_ref, target, max_modBkg, third_oband, prop_loss_cum, prop_loss_indiv)

    while abs(D6 - D5) >= precision_fraction * detection_dist:
        detection_dist = (D5 + D6) / 2
        M2, B1, prop_loss_cum, prop_loss_indiv = Propagate(freq, detection_dist, ground_effect_ref, atm_abs, atm_abs_ref, target, max_modBkg, third_oband, prop_loss_cum, prop_loss_indiv)
        if M2 > 0:
            D5 = detection_dist
        else:
            D6 = detection_dist

    #return 5 items - detection distance most important
    return detection_dist