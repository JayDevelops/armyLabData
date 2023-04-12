import pandas as pd  # Install pandas in python or use Anaconda environment
import math
from DetectionScripts import data_dict
from DetectionScripts import ADM_Functions


def ansi_humidity(freq):
    rel_hum_H1 = data_dict.det_cons['relative_humid_percent_det']
    R = data_dict.given_cons['m_measure_distance']
    atmos_absorption = []
    kelvin_Too = 273.15
    Co = 331.32 
    celsius_degrees_T1 = data_dict.det_cons['celsius_degrees_det']
    t = kelvin_Too + celsius_degrees_T1
    PSTAR = 1
    acr1 = 0.000000016 * (1.377 * t / (t + 110.4)) / PSTAR
    thto = 2239.1 / t
    exthto = math.exp(-thto)
    rto = thto / (1 - exthto)
    mumaxo = 0.20948 * 4 * math.pi * rto * rto * exthto / 35
    thtn = 3352 / t #in the VBA code (line 57), 3352 has a "!" after it, not sure why
    exthtn = math.exp(-thtn)
    rtn = thtn / (1 - exthtn)
    mumaxn = 0.78084 * 4 * math.pi * rtn * rtn * exthtn / 35
    cs = Co * math.sqrt(t / kelvin_Too)
    h = rel_hum_H1 * math.exp(math.log(10) * (20.5318 - 2939 / t - 2.13759744 * math.log(t))) / PSTAR #H1 is RH% from excel
    fox = (24 + 44100 * h * (0.05 + h) / (0.391 + h)) * math.sqrt(293 / t)    # Is 293 used as 273 + 20degC?
    fni = (9 + 350 * h) * (293 / t)
    d1 = R * 0.01 #R can be a different value depending on what helper function is called

    for i in range(24):
        freq2 = freq[i] #freq is from the Freq. Hz collumn, need to maybe make a df for this?
        freq2 = freq2 * freq2
        acr = acr1 * freq2
        x2o = freq2 / (fox * fox)
        avibo = 868.5 * mumaxo * fox / cs * x2o / (1 + x2o)
        x2n = freq2 / (fni * fni)
        avibn = 868.5 * mumaxn * fni / cs * x2n / (1 + x2n)
        atmos_absorption.append(-d1 * (acr + avibo + avibn)) #AtmAbs, atmosphere absorption, a global variable this is where it is initialized
    
    return atmos_absorption