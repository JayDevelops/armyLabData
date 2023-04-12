import pandas as pd  # Install pandas in python or use Anaconda environment
import math
from DetectionScripts import data_dict
from DetectionScripts import ADM_Functions


def calculate_measure_dist(detection_dist: float):
    return detection_dist * 25.0