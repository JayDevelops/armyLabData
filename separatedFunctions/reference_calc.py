import pandas as pd  # Install pandas in python or use Anaconda environment
import math
from DetectionScripts import data_dict
from DetectionScripts import ADM_Functions

def reference_calc(freq):
    ground_effect = Ingard(freq)
    atmos_absorption = ansi_humidity(freq)
    AtmAbsRef = []
    ground_effect_ref = []
        
    for I in range(0,24):
        ground_effect_ref.append(ground_effect[I])
        AtmAbsRef.append(atmos_absorption[I])
    
    return ground_effect, atmos_absorption, ground_effect_ref, AtmAbsRef