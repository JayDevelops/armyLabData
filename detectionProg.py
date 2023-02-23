""" Based on the meeting we had on Monday, where we decided it is best
to work on making a small, simple program using values that we know, as
opposed to trying to decipher all the macros, this is the program I think 
we can build on. These values were obtained from the ADM from Joel model
as output when the given parameters were run through it. """

#This is assuming a mic distance of 30m, which we will need to address later
#will first use the measurement conditions on the sheet so we can get an accurate model
#Then we can try with the drone data

import math

freq, bkg, targetNoise = 0
hearing_thresh = 0

