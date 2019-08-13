import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.special as scp
import scipy.integrate as int
import mpmath as mp
dracoBackground = 2.32e5 #Predicted background from Lines and Boxes.
fiveSigmaSignals = 2.41e3 #num req. signals from line and box paper.
diffuseDracoContribution = 2.74e-3 #Constant from Paper.
omega = 1.6e-3 # Solid Angle
effectiveArea=500 # Just as it says, the effective area in centimeters squared.
duration = 31536000 # Duration of observation time.

# This function accounts for the stochastic variation in the background data by throwing the number of expected events
# into a Poisson distribution.
# This function takes in the following:
# binStart - this is the starting energy of the bin whose background count is to be generated.
# binEnd - this is the ending energy of the bin whose background count is to be generated.
def poissonSelection(binStart,binEnd):
    if binStart == 0 or binEnd == 0:
        return 0
    else:
        numExpectedEvents = findExpectedEventsDraco(binStart,binEnd)
        numObserved = np.random.poisson(numExpectedEvents)
        return numObserved

# This function returns the expected number of background events for a given bin.
# This function takes in the following:
# binStart - This is the starting energy of the bin whose background count is to be generated.
# binEnd - This is the ending energy of the bin whose background count is to be generated.
def findExpectedEventsDraco(binStart,binEnd):
    if binStart == 0 or binEnd == 0:
        return 0
    else:
        return omega*effectiveArea*duration*diffuseDracoContribution*((binStart)**-1 - (binEnd)**-1)