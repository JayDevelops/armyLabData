import pandas as pd  # Install pandas in python or use Anaconda environment
import math
from DetectionScripts import data_dict
from DetectionScripts import ADM_Functions

def ListenerdBA():
    E = float()
    E = 0.0
    for I in range(0, 23):
        E = E + math.exp(Log10Div10 * (target_spectrum[I] + prop_loss_cum[I] + A_weight_levels[I]))
    listener_dba_return = TenDivLog10 * math.log(E)
    return listener_dba_return