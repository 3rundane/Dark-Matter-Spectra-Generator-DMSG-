import math
import numpy as np
class Bin:
    # Width is the energy range of particles you want your bin to accept.
    # start is the smallest energy value that I bin can accept.
    # Using start and width the bin object then determines the last energy it can accept by adding Width to Start.
    # backgroundCounts is the number of background events you want the bin to start out its life with.
    # signalCount is the number of signal events you want the bin to start out its life with.
    # totalCount is the total number of event counts (background + signal) that you want your bin to start out with.
    # By default backgroundCounts, signalCounts, and totalCount start out at 0.

    #The bin object has several class members that it stores. These are:
    # binWidth, binStart, binEnd, backgroundCounts, signalCount, and totalCount
    # binWidth is the energy width of the bin.
    # binStart is the starting energy of the bin.
    # binEnd is the ending energy of the bin.
    # backgroundCounts is the current number of background events stored in the bin.
    # signalCount is the current number of signal events stored in the bin.
    # totalCount is the current number of total events in stored in the bin.

    # The bin object has several class methods available for use. These are, in top-down order:
    # addBackgroundEvent, removeBackgroundEvent, setBackgroundCount,getBackgroundCount,setSignalCount,getSignalCount,getSignalCountDirectly,
    # addSignalCounts,setTotalCount,getTotalCount,addTotalCount,updateTotalCount,acceptEnergy,getStartingEnergy,EndingEnergy,getBinWidth, smearSignal
    # I will describe the methods and their functions above the actual function definition.
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


    # increments the number of background events contained in this bin. I typically used this method in some for loop.
    def addBackgroundEvent(self): #unittest
        self.__backgroundCounts = self.__backgroundCounts + 1
        self.updateTotalCount()
    # decrements the number of background events contained in this bin. It screams at you if you try to make it negative.
    def removeBackgroundEvent(self): #unittest
        if self.__backgroundCounts > 0:
            self.__backgroundCounts = self.__backgroundCounts - 1
        else:
            print("Count is 0. You cannot have negative counts")
        self.updateTotalCount()

    #This allows the user to directly set the number of background counts at any time. I didn't use it often but if you want it it's there.
    def setBackgroundCount(self,numBackgroundCounts): #unittest
        if numBackgroundCounts > 0:
            self.__backgroundCounts = math.floor(numBackgroundCounts)
        self.updateTotalCount()
    #This returns the number of background counts "contained" in the bin at the present time. I mainly used this for plotting purposes.
    def getBackgroundCount(self): #unittest
        return self.__backgroundCounts

    # This does the same thing as setBackgroundCount above but with signal counts. It doesn't scream at you however if you try to make counts negative.
    # it also calls updateTotalCount so that the bin object is always up to date with how many total counts there are present.
    def setSignalCount(self,numSignalCounts): #unittest
        if numSignalCounts > 0:
            self.__signalCount = math.floor(numSignalCounts)

        self.updateTotalCount()

    # This method is meant to return the observed signal count to the user. Namely the signal count subject to random variation.
    # This is accomplished by taking signalCount as the mean of a poisson distribution returning a draw from that distribution.
    # The reason we draw from a poisson distribution is because the code tries to mimic a telescope's storage process,
    # and fundamentally,a telescope is really a counting experiment with a low probability of success and high number of trials.
    def getSignalCount(self):
        # return self.__signalCount
        return np.random.poisson(self.__signalCount)

    # This method does the same thing as above but gives the actual signal count stored, i.e it returns an exact value and not one
    # subject to random variation.
    def getSignalCountDirectly(self):
        return self.__signalCount

    #This increments the number of signal counts and updates total counts.
    def addSignalCounts(self):
        self.__signalCount = self.__signalCount + 1
        self.updateTotalCount()

    #sets the total number of counts in the bin object.
    def setTotalCount(self,totalCount):
        if totalCount > 0:
            self.__totalCount = math.floor(totalCount)

    #returns the total number of counts in the bin object.
    def getTotalCount(self):
        return self.__totalCount

    # Increments the total number of counts in the bin object.
    def addTotalCount(self):
        self.__totalCount = self.__totalCount + 1

    # method is used to update total counts whenever a user operates on either signal counts or backgroundCounts
    def updateTotalCount(self):
        if self.__backgroundCounts <= self.__signalCount:
            self.__totalCount = self.__signalCount - self.__backgroundCounts
        #else:
            #print('Error: Background count is greater than signalCount')



    #This method checks whether or not an energy should be added based upon if it is contained within the bin range (binStart and binEnd).
    #The method returns True or False.
    def acceptEnergy(self,energyToAdd):
        smearedEnergy = self.smearSignal(energyToAdd,0.01) # this is where I tried to implement the energy resolution smearing.
        # smearedEnergy=energyToAdd
        if (smearedEnergy >=self.__binStart and self.__binEnd >smearedEnergy):
            return True
        else:
            return False

    #Returns the starting energy of the bin.
    def getStartingEnergy(self): #unittest
        return self.__binStart

    #Returns the ending energy of the bin.
    def getEndingEnergy(self): #unittest
        return self.__binEnd

    #Returns the energy width of the bin.
    def getBinWidth(self): #unittest
        return self.__binWidth

    # This is my attempt at a method that can smear signals for energy resolution considerations. You pass in an energy, then it returns a smeared signal to check against
    # the bin's energy range. If the smeared energy lies within the bin's range then we finally add a count.
    # acceptedEnergy is meant to be the energy that wants to be added to the bin.
    # epsilon is contributes to the variance of the gaussian function.
    def smearSignal(self,acceptedEnergy, epsilon):
        ##does this even do what we want it to do? does it even have the right units?
        smearedSignal = np.random.normal(acceptedEnergy, epsilon * acceptedEnergy)
        return smearedSignal

    # def signalList(self):
    #     return self.signalList #parantheses self.signalList() or not self.signalList??

