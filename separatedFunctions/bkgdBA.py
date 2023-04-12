import pandas as pd  # Install pandas in python or use Anaconda environment
import math
from DetectionScripts import data_dict
from DetectionScripts import ADM_Functions

def BkgdBA():
     E = float()
     E = 0.0
     
     for I in range(0,23):
         E = E + math.exp(Log10Div10 * (prop_loss_cum[I] + A_weight_levels[I]))
         
     bkgdba_return = TenDivLog10 * math.log(E)
     
     return bkgdba_return    