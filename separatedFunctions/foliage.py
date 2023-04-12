import pandas as pd  # Install pandas in python or use Anaconda environment
import math
from DetectionScripts import data_dict
from DetectionScripts import ADM_Functions

def Foliage(freq, prop_loss_cum, prop_loss_indiv, detection_distance):
    Fo = [-0.001] * 24
    Fo_list = []
    foliated_zone_nums = data_dict.det_cons['foliage_on']
    distance_frmSource_to_fol = data_dict.det_cons['foliage_dist']
    foliage_depth = data_dict.det_cons['foliage_depth']
    leaf_width = data_dict.det_cons['leaf_width']
    leaf_areapervol = data_dict.det_cons['leaf_areapervol']
    celsius_degrees_det = data_dict.det_cons['celsius_degrees_det']
    Cs = 331.4 * math.sqrt(1 + celsius_degrees_det / 273.15)

    if foliated_zone_nums > 0:
        # If detection distance is greater than the distance from source to edge of foliage
        if detection_distance > distance_frmSource_to_fol:
            # Sets the distance difference to X2
            X2 = detection_distance - distance_frmSource_to_fol
            # If the new detection distance is greater than depth (extent) of foliage in meters
            if X2 > foliage_depth:
                # Set distance the depth (extent) of foliage in meters
                X2 = foliage_depth
            X2 = math.sqrt(X2)
            Cons = 2.647 / math.log(10)

            for I in range(10, 24):
                Ka = (2 * math.pi * freq[I] / Cs) * leaf_width / 100
                if Ka < 0.401:
                    Fo[I] = -0.01
                elif Ka < 5:
                    Fo[I] = -X2 * math.sqrt(leaf_areapervol) * (Cons * math.log(Ka) + 1.05)
                else:
                    Fo[I] = -X2 * math.sqrt(leaf_areapervol) * 2.9
                Fo_list.append(Fo[I])
    else:
        for I in range(24):
            prop_loss_indiv[I] = Fo[I]
            prop_loss_cum[I] = prop_loss_cum[I] + prop_loss_indiv[I]
            Fo_list.append(prop_loss_cum[I])
    print('Foliage', prop_loss_cum)
    return prop_loss_cum, prop_loss_indiv