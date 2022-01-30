
import numpy as np
from DMSG import DMSG_shmearStuff as smear
from scipy import signal
from DMSG import DMSG_DDMGenerationLoops as dg
# import DMSG_DDMGenerationLoops as dg

# dracoBackground is the total number of expected background events observed in a year as outlined by the Line and Box paper.
# totalSignals is the total number of signals from dark matter needed to claim a five-sigma discovery as outlined by the Line and Box paper.
# diffuseDracoContribution is a constant put out front in eq. 4.2 in the Line and Box paper.
# numBins is the same number of bins passed into Generation's constructor. It's purpose here is to save the value
# that was passed in just in case it is needed.
# binWidth is the same binWidth passed into this class.
# backgroundCountList is a list of the number of background signal counts of the bin objects
# signalCountList is a list of the signalCounts for each bin object that a generation object is holding.
# energyList is a list of the starting energies of each bin, ordered from smallest energy to largest.
# signalEnergyStorage is a list that holds the  sorted energies of each signal event generated. It is empty at first and is made to
# be thrown into a function from whatever dark matter generation model that both generates the events/energies and sorts them.
class Generation:
    # __dracoBackground = 2.32e5
    # __totalSignals = 2.41e3  # num req. signals from line and box paper
    # __diffuseDracoContribution = 2.74e-3  # inverse squared events per MeV) I think. It's eq. 4.2 in Line and Box Paper.
    __numBins = 0
    __binWidth = 2
    __backgroundCountList = []
    __signalCountList = []

    __energyList = []
    __signalEnergyStorage = []

    # numBins is the total number of bins that the user wants to create. Typically I keep this near to the number of particles.
    # startEnergy is the starting energy of the first bin desired. Typically this is used to set the "origin" of the energy axis in spectra plots.
    # binWidth is the width of each energy bin. Note that the width of each bin is assumed to be the same for all bins in this code.
    # totalNum is the total number of particles in the ensemble.
    # totalSignals is the total number of signals that we need to generate in order to claim a 5 sigma discovery for a given background.
    # sqrtSNaught is the center of mass energy of the lightest particle in teh ensemble.
    # delta is a shape parameter used in equation 2.7 of the line and box paper
    # deltaSSquare is the energy splitting between states, here it is assumed to be a constant.
    # eta is the shape parameter used throughout the line and box paper. I'd list an equation but it's found in so many of them that it would be redundant.
    def __init__(self,numBins,startEnergy,binWidth,totalNum,totalSignals,sqrtSNaught,delta,deltaSSquare,eta):
        ########
        #initialize any necessary constants for below.



        #########
        #calculate relevant quantities
        self.__binWidth = binWidth
        self.__endEnergyList = [0]*numBins
        self.__totalSignals = totalSignals
        self.__signalCountList = [0]*numBins
        self.__binList = [0]*numBins
        self.__binList = dg.binCreatorBG(self.__binList,startEnergy,binWidth,numBins)
        self.__energyList = [0]*numBins
        self.__backgroundCountList = [0] * numBins
        self.__signalEnergyStorage,self.__signalSigmas = dg.createSignalEnergies(totalNum,sqrtSNaught,totalSignals,delta,deltaSSquare,eta)

        # the code below deals with implementing the effects of telescope energy resolution on the spectra. However the functions
        # I have written to smear the energies are not currently working which is why I have commented out this feature.
        # x = np.arange(-100, 100)
        # y = np.array([smear.gaussianSmearingFunction2(t, 0.5) for t in x])
        # filtered2 = signal.convolve(self.__signalEnergyStorage, y, mode='same') / sum(y)
        # self.__signalEnergyStorage = filtered2

        # operates on your bins to update their counts. After the line below the bins will contain all of count data.
        dg.placeEnergiesMk2(self.__totalSignals,self.__signalEnergyStorage,self.__binList)



    #Returns a list of the bin objects stored in this Generation instance.
    def getBinList(self):
        return self.__binList

    #Returns a list containing the number of background events for each energy bin.
    def getBackgroundBinCount(self):
        for i in range(0,len(self.__binList)):
            self.__backgroundCountList[i] = ((self.__binList[i]).getBackgroundCount()) / self.__binWidth #TEST MODULE
        return self.__backgroundCountList

    #Returns a list containing the number of signal events for each energy bin.
    def getSignalBinCount(self):
        for i in range(0,len(self.__binList)):
            #print(self.__binList[i].getSignalCount())
            self.__signalCountList[i] = ((self.__binList[i]).getSignalCount()) / self.__binWidth #TEST MODULE
        return self.__signalCountList

    #Returns a list containing the ordered starting energies for each bin. Typically I use this list for plotting.
    def getBinEnergy(self):
        for i in range(0, len(self.__binList)):
            self.__energyList[i] = (self.__binList[i]).getStartingEnergy() #TEST MODULE
        return self.__energyList
        # inherited classes will take care of initializing the number number of expected and observed counts

    #Returns a list containing the ordered ending energies for each bin.
    def getEndingEnergy(self):
        for i in range(0,len(self.__binList)):
            self.__endEnergyList[i] = (self.__binList[i]).getEndingEnergy()
        return self.__endEnergyList

    #Gives the uncertainty list. This might be redundant...
    def getUncertaintyList(self):
        #this might be incorrect...
        return self.__signalSigmas
    #Gives you the original generated energies for the photons.
    def getSignalEnergyStorage(self):
        return self.__signalEnergyStorage
        # new processing functionality for background.

