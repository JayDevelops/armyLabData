import pandas as pd  # Install pandas in python or use Anaconda environment
import math
from DetectionScripts import data_dict
from DetectionScripts import ADM_Functions

def Propagate(freq, detection_distance, ground_effect_ref, atm_abs, atm_abs_ref, target, max_modBkg, third_oband, prop_loss_cum, prop_loss_indiv):
    prop_loss_cum, prop_loss_indiv = inverse_distance(detection_distance, prop_loss_cum, prop_loss_indiv)
    prop_loss_cum, prop_loss_indiv = ground_effect(freq, prop_loss_cum, prop_loss_indiv, ground_effect_ref)
    prop_loss_cum, prop_loss_indiv = Barrier(freq, prop_loss_cum, prop_loss_indiv, detection_distance)
    prop_loss_cum, prop_loss_indiv = Foliage(freq, prop_loss_cum, prop_loss_indiv, detection_distance)
    prop_loss_cum, prop_loss_indiv = atmosphere(freq, prop_loss_cum, prop_loss_indiv, atm_abs, atm_abs_ref)
    M2, B1, prop_loss_cum = signal_noise(target, max_modBkg, third_oband, prop_loss_cum)

    return M2, B1, prop_loss_cum, prop_loss_indiv