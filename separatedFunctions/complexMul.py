import pandas as pd  # Install pandas in python or use Anaconda environment
import math
from DetectionScripts import data_dict
from DetectionScripts import ADM_Functions

def ComplexMul(A, B, C, D):
    E = A * C - B * D
    F = A * D + B * C
    return E, F