import math
import numpy as np

class Bin:
    def __init__(self, width, start, backgroundCounts=0, signalCount=0, totalCount=0):

        self.__binWidth=width
        self.__binStart=start
        self.__binEnd = start + width
        self.__backgroundCounts=math.floor(backgroundCounts)
        self.__signalCount =int(signalCount)
        self.__totalCount = int(totalCount)

        if self.__signalCount >= self.__backgroundCounts:
            self.__totalCount = signalCount - backgroundCounts

    def addBackgroundEvent(self):
        self.__backgroundCounts = self.__backgroundCounts + 1
        self.updateTotalCount()
    def removeBackgroundEvent(self):
        if self.__backgroundCounts > 0:
            self.__backgroundCounts = self.__backgroundCounts - 1
        else:
            print("Count is 0. You cannot have negative counts")
        self.updateTotalCount()

    ##Operate on background Counts
    def setBackgroundCount(self,numBackgroundCounts):
        if numBackgroundCounts > 0:
            self.__backgroundCounts = math.floor(numBackgroundCounts)
        self.updateTotalCount()
        
    def getBackgroundCount(self):
        return self.__backgroundCounts

    ## Operate on signalCounts
    def setSignalCount(self,numSignalCounts):
        if numSignalCounts > 0:
            self.__signalCount = math.floor(numSignalCounts)

        self.updateTotalCount()
    def getSignalCount(self):
        # return self.__signalCount
        return np.random.poisson(self.__signalCount)

    def addSignalCounts(self):
        self.__signalCount = self.__signalCount + 1
        self.updateTotalCount()

    ## Operate on TotalCount

    def setTotalCount(self,totalCount):
        if totalCount > 0:
            self.__totalCount = math.floor(totalCount)

    def getTotalCount(self):
        return self.__totalCount

    def addTotalCount(self):
        self.__totalCount = self.__totalCount + 1

    def updateTotalCount(self):
        if self.__backgroundCounts <= self.__signalCount:
            self.__totalCount = self.__signalCount - self.__backgroundCounts

    def acceptEnergy(self,energyToAdd):
        smearedEnergy = self.smearSignal(energyToAdd,0.01)
        if (smearedEnergy >=self.__binStart and self.__binEnd >smearedEnergy):
            return True
        else:
            return False

    def getStartingEnergy(self):
        return self.__binStart

    def getEndingEnergy(self):
        return self.__binEnd

    def getBinWidth(self):
        return self.__binWidth

    def smearSignal(self,acceptedEnergy, epsilon):
        smearedSignal = np.random.normal(acceptedEnergy, epsilon * acceptedEnergy)
        return smearedSignal

