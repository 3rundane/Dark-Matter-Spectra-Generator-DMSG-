from DMSG import DMSG_Generation,DMSG_Bin
import numpy as np
import math

#This function is meant to simply compute the standard deviation of a poisson distribution whose mean is given as
# as a list of arguments (signalCount).
# This returns a list of standard deviations computed by iterating through signalCount and using its elements as inputs for
# poisson distributions to draw from. It is then square-rooted to obtain the sample standard deviation. Note I may have
# done the statistics wrong here.
def computeSigma(signalCount):
    lengthOfList = len(signalCount)
    yourSigma = [0]*lengthOfList
    for i in range(lengthOfList):
        if(signalCount[i]==0):
            yourSigma[i]=1
        else:
            yourSigma[i] = math.sqrt(np.random.poisson(signalCount[i]))
            # yourSigma[i] = math.sqrt(signalCount[i])
    return yourSigma