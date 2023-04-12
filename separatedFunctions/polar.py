import pandas as pd  # Install pandas in python or use Anaconda environment
import math
from DetectionScripts import data_dict
from DetectionScripts import ADM_Functions

def Polar(X, Y):
    R = math.sqrt(X * X + Y * Y)
    Th = math.atan2(X, Y)
    return R, Th
