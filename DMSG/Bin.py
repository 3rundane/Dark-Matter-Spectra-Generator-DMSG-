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
       # else:
            #print("Error: backgroundCounts greater than signalCounts")



    def addBackgroundEvent(self): #unittest
        self.__backgroundCounts = self.__backgroundCounts + 1
        self.updateTotalCount()
    def removeBackgroundEvent(self): #unittest
        if self.__backgroundCounts > 0:
            self.__backgroundCounts = self.__backgroundCounts - 1
        else:
            print("Count is 0. You cannot have negative counts")
        self.updateTotalCount()
    ##Operate on background Counts
    def setBackgroundCount(self,numBackgroundCounts): #unittest
        if numBackgroundCounts > 0:
            self.__backgroundCounts = math.floor(numBackgroundCounts)
        self.updateTotalCount()
    def getBackgroundCount(self): #unittest
        return self.__backgroundCounts

    ## Operate on signalCounts
    def setSignalCount(self,numSignalCounts): #unittest
        if numSignalCounts > 0:
            self.__signalCount = math.floor(numSignalCounts)

        self.updateTotalCount()
    def getSignalCount(self): #unittest
        # return self.__signalCount
        return np.random.poisson(self.__signalCount)
    def addSignalCounts(self): #unittest
        self.__signalCount = self.__signalCount + 1
        self.updateTotalCount()

    ## Operate on TotalCount
    # need to update background and signal count each time one of the methods below are run
    #you need to code this further Dane.
    def setTotalCount(self,totalCount): #unittest
        if totalCount > 0:
            self.__totalCount = math.floor(totalCount)
    def getTotalCount(self): #unittest
        return self.__totalCount
    def addTotalCount(self): #unittest
        self.__totalCount = self.__totalCount + 1
    def updateTotalCount(self): #unittest
        if self.__backgroundCounts <= self.__signalCount:
            self.__totalCount = self.__signalCount - self.__backgroundCounts
        #else:
            #print('Error: Background count is greater than signalCount')




    def acceptEnergy(self,energyToAdd):
        smearedEnergy = self.smearSignal(energyToAdd,0.01)
        if (smearedEnergy >=self.__binStart and self.__binEnd >smearedEnergy):
            return True
        else:
            return False


    def getStartingEnergy(self): #unittest
        return self.__binStart

    def getEndingEnergy(self): #unittest
        return self.__binEnd

    def getBinWidth(self): #unittest
        return self.__binWidth

    def smearSignal(self,acceptedEnergy, epsilon):
        ##does this even do what we want it to do? does it even have the right units?
        smearedSignal = np.random.normal(acceptedEnergy, epsilon * acceptedEnergy)
        return smearedSignal

    # def signalList(self):
    #     return self.signalList #parantheses self.signalList() or not self.signalList??

