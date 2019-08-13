import DracoBackground as dar
import Bin
import random as rnd
import SortingAlgorithms as sort
import numpy as np
import matplotlib.pyplot as plt
import math

#This function takes in:
# an empty list with which to store bins (binList).
# the starting energy of the very first bin to be created. (startEnergy)
# the width for all the bins. (In here I assume the bin width to be the same for all bins.) (binWidth)
# the number of bins desired (numBins)
# the function returns a non-empty (in principle) list full of bin objects which are ordered in terms of increasing starting energy. Which is the same thing
# as ordering them by increasing ending energy. The convention I've taken for the code however is to use the starting energies for pretty much all functionality.
# Another thing to note about this function is that it also implements the background model the user has "selected". Right now the model that is "plugged in" is derived from
# the fit to COMPTEL and EGRET data that was used in the Lines and Boxes paper by Boddy and al. I believe that it is equation 4.2. Essentially I integrated that differential flux fit
# (4.2) with respect to bin energy and the solid angle for both the solid angle given in the paper and two arbitrary starting and ending energies. That general function is then
# multiplied by the effective area used in that same paper as well as the duration of the observation time in seconds. This function is stored in DracoBackground and is used to
# generate for each bin the number of background counts. Note that the analytic solution to the integration mentioned above is important as the limits of integration for the bin
# energy depends on the energy range of the bin in question.
def binCreatorBG(binList, startEnergy, binWidth, numBins):
        for i in range(numBins):
            binList[i] = Bin.Bin(binWidth, startEnergy + i*(binWidth),dar.poissonSelection(startEnergy+i*binWidth,startEnergy+(i+1)*binWidth))
        return binList

#This function's entire purpose in life is to help findNumSig compute the total photon flux via the form I derive in my notes.
def computeDenom(denom, totalNum, delta, deltaSSquare, sqrtSNaught, eta):
    for n in range(1,totalNum+1):
        denom = denom + (1 + ((n) ** delta) * (deltaSSquare / sqrtSNaught)) ** eta
    return denom

# This function "places" the generated signals via the energies in the proper bins from which those bins then update their counts then returns the list of updated bins to the user.
# Note that I use Loop5 in my code because it is note very useful to return the binList at the point where I use these methods. I am only keeping this one here in case someone needs
# it. The function takes a total number of signals determined from the background signal amount, a list of energy values ord
def placeEnergiesMk1(totalSignals, signalEnergyStorage, binList):
        signalToPlace = 0
        for i in range(len(signalEnergyStorage)):
            signalToPlace = signalEnergyStorage[i]
            for j in range(len(binList)):
                if binList[j].acceptEnergy(signalToPlace):
                    binList[j].addSignalCounts()
        return binList

# In the Lines and Boxes paper by Boddy and Al. They consider the number of signals a detector would need to detect from Draco in order to claim a 5-sigma discovery. Furthermore
# they assert that for each dark matter particle chi in the ensemble of states the ratio of the number of signals form chi to the total signals should be proportional to the ratio
# of that particle's flux to the total ensemble flux. This gives us N sub s / N tot. is proportional to capital Phi sub n / captial Phi tot. Substituting in equations
# 2.7 and 3.7 from their paper gives us a formula to directly calculate the number of signals from each dark particle we need to generate as a function of delta, the number of particles in
# the ensemble, the mass splitting and the lightest state's center of mass energy value. This function is implemented below.
def findNumSig(n,totalSignals,sqrtSNaught,delta,deltaSSquare,eta,totalNum):
    currentSN = sqrtSNaught + ((n) ** delta) * deltaSSquare
    numer = (1 + (n ** delta) * (deltaSSquare / sqrtSNaught)) ** eta
    denom = 0
    denom = computeDenom(denom, totalNum, delta, deltaSSquare, sqrtSNaught, eta)
    numsig = int(round(totalSignals * (numer / denom))) #compute the fraction of the total number of signals required to claim a 5-sigma discovery given Draco's total background signals.
    return numsig

# This function creates separate lists of different length to be passed into fillOutArrays().
# The problem that createStorageDictionary was designed to solve is caused by there being multiple particles with different
# amounts of signals coming from them combined with the condition where roughly 1/3 of those signals for each type of particle in the
# ensemble must have an energy according to equation 3.11 with the rest being uniformly distributed between the two energies given by
# equation 3.15. Since there is, in principle, 1 list for each particle type with all of these lists being potentially different lengths I created
# a dictionary who holds the lengths of the nth particle's list. This allowed me to mizimize the number of for loops used to create the signal energies.
# This function takes in:
# totalSignals - this is the number of signals needed for discovery
# sqrtSNaught - this is the center-of-mass energy of the lightest particle in the ensemble.
# delta - this is the shape parameter in equation 2.7
# deltaSSquare - this is the energy splitting
# eta - this one of the shape parameters in Lines and Boxes.
# totalNum - this is the total number of particles in the ensemble.
def createStorageDictionary(totalSignals,sqrtSNaught,delta,deltaSSquare,eta,totalNum):
    n=1
    x= np.arange(1,totalNum+1).tolist()
    dct = {}
    for j in x:
        dct['lst_%s' %j] = [0]*findNumSig(n,totalSignals,sqrtSNaught,delta,deltaSSquare,eta,totalNum)
        n = n + 1
    return dct

# This function fills out an unordered list of energies and signal uncertainties then returns them to the user.
# The function takes in the following:
# dct - a dictionary containing the lengths of each list for each particle in the ensemble.
# nrange - This is a list containing values from 1 to N where N is the number of particles and each element is greater than the previous by 1.
# newList - This is an empty list which will be used to contain the ordered energies.
# sqrtSNaught - This is the center-of-mass energy of the lightest particle in the ensemble.
# delta - This is the shape parameter in equation 2.7
# deltaSSquare - This is the energy splitting.
def fillOutArrays(dct,nrange,newList,sqrtSNaught,delta,deltaSSquare):
    currentSN = 0
    eBoxMin = 0
    eBoxMax = 0
    eLine = 0
    signalSigmas=[1]*len(newList)
    for n in nrange:
        currentSN = sqrtSNaught + ((n) ** delta) * deltaSSquare
        eBoxMin = (135 ** 2) / (2 * currentSN)
        eBoxMax = currentSN / 2 ##still need to test these
        eLine = ((currentSN ** 2 - 135 ** 2) / (2 * currentSN))

        nameOfList = ('lst_%s' % n)
        currentList = []
        currentSigmas=[]
        numberOfSignals = 0
        currentList = dct.get(nameOfList)
        numberOfSignals = len(currentList)
        for i in range(len(dct.get(nameOfList))):
            currentSigmas = [1]*numberOfSignals

            if i <=numberOfSignals/3:
                # currentSigmas[i] = eLine
                # # currentList[i] = np.random.poisson(eLine)
                currentList[i] = eLine

            else:
                currentSigmas[i] = ((eBoxMax-eBoxMin)**2)/12
                currentList[i] = rnd.uniform(eBoxMin, eBoxMax)

        newList = newList + currentList
        signalSigmas = signalSigmas + currentSigmas
    return newList,signalSigmas

# This functions computes a sorted list of signal energies to be placed into some bin list. It also returns a
# list of uncertainties associated with the signals. The following arguments that are passed in are:
# signalEnergyStorage - This is an empty list that will store the generated energies in increasing order.
# totalNum - This is the total number of particles in the ensemble.
# sqrtSNaught - This is the center-of-mass energy of the lightest particle
# totalSignals - This si the total number of signals needed to claim discovery for a given background.
# delta - This is the shape parameter in equation 2.7
# deltaSSquare - this is the energy splitting
# eta - This is the shape parameter used throughout Boddy and Al.

def createSignalEnergies(totalNum, sqrtSNaught, totalSignals, delta, deltaSSquare, eta): #signalEnergyStorage used to be the first argument.
    dct = createStorageDictionary(totalSignals,sqrtSNaught,delta,deltaSSquare,eta, totalNum)
    nrange = np.arange(1,totalNum+1).tolist()
    newList = []
    newList,signalSigmas= fillOutArrays(dct,nrange,newList,sqrtSNaught,1,2)

    return sorted(newList), sorted(signalSigmas)

# This function "places" the generated signals via the energies in the proper bins from which those bins thens update their counts then returns the list of updated bins to the user.
# The function takes in the following arguments:
# totalSignals - The total signals needed to claim discovery for a given background.
# signalEnergyStorage - This is a NONEMPTY list of the ordered energies.
def placeEnergiesMk2(totalSignals, signalEnergyStorage, binList):
    signalToPlace = 0
    for i in range(len(signalEnergyStorage)):
        signalToPlace = signalEnergyStorage[i]
        for j in range(len(binList)):

            # the algorithm is designed to add signals to the first bin it fits inside of. This could cause some design flaws so I'll double check this later.
            if binList[j].acceptEnergy(signalToPlace):
                binList[j].addSignalCounts()

