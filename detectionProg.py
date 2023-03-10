""" Based on the meeting we had on Monday, where we decided it is best
to work on making a small, simple program using values that we know, as
opposed to trying to decipher all the macros, this is the program I think 
we can build on. These values were obtained from the ADM from Joel model
as output when the given parameters were run through it. """

#This is assuming a mic distance of 30m, which we will need to address later
#will first use the measurement conditions on the sheet so we can get an accurate model
#Then we can try with the drone data

import math

"""The variables that are used in this function are described as:

barrier_attenuation -> Global array variable (declared in line 35 of VBA Macros). Inititalized for first time inside Barrier()
detection_dist -> Global 'double' variable (declared in line 13 of VBA Macros). It is also used in GroundEffect() on line 186
distance_from_source -> initialized from excel cell "I5" (line 579)
barrier_height_det -> initialized from excel cell "J5" (line 580)
source_height_det -> initialized from excel cell "H3" (line 564)
listener_height_det -> initialized from excel cell "I3" (line 565)
celsius_degrees_det -> initialized from excel cell "J3" (line 566)
TenDivLog10 -> Static number (described in lines 546 and 547)
freq -> Global integer array declared in line 21. Looks like initialized in InitMacros as 'A' column values of Data Sheet (aka freq Hz)
prop_loss_indiv -> Global 'double' array declared in line 25. Is used in many functions, and changes based on the actions of those functions
prop_loss_cum -> Global 'double' array declared in line 20. Is used to hold value after all changes to prop_loss_indiv throughout the program.  
"""

def Barrier():
    for i in range(24):
        barrier_attenuation(i) = -0.001
    if barrier_number > 0:
        if detection_dist >= distance_from_source:
            Hba = barrier_height_det - (source_height_det + distance_from_source * (listener_height_det) / detection_dist)
            X1 = math.sqrt(distance_from_source ** 2 + ((barrier_height_det - source_height_det) **2))
            X1 = X1 + math.sqrt(((detection_dist - distance_from_source) ** 2) + ((barrier_height_det - listener_height_det) ** 2))
            X1 = X1 - math.sqrt((detection_dist **2) + ((listener_height_det - source_height_det)**2))
            Cs = 331.4 * math.sqrt(1 + celsius_degrees_det / 273.15)
            
            for i in range(24):
                if(Hba < 0):
                    X2 = 2 * X1 * freq(2) / Cs
                    X3 = math.sqrt(2 * Pi * abs(X2))
                    X2 = -(X2)
                    if(X2 <= -0.1916):
                        barrier_attenuation(i) = -0.01
                    else: 
                        barrier_attenuation(i) = -(5 + 2 * TenDivLog10 * math.log(X3 / math.tan(X3))) - 0.01
                else:
                    X2 = 2 * X1 * freq(i) / Cs
                    X3 = math.sqrt(2 * Pi * abs(X2))
                    X2 = math.exp(2 * X3)
                    if(X2 == 1):
                        barrier_attenuation(i) = -5 -0.01
                    else:
                        barrier_attenuation(i) = -(5 + 2 * TenDivLog10 * math.log(X3 * (X2 + 1) / (X2 - 1))) - 0.01
                        if(barrier_attenuation(i) < -20):
                            barrier_attenuation(i) = -20
    
    for i in range(24):
        prop_loss_indiv(i) = barrier_attenuation(i)
        prop_loss_cum(i) = prop_loss_cum(i) + prop_loss_indiv(i)
                        
                        
        
    
