import numpy as np
import matplotlib.pyplot as plt
import math as mt
import scipy.special as scp
import scipy.integrate as int
from IndirectMain import Bin

class Generator:
    __dracoBackground = 2.32e5
    __fiveSigmaSignals = 2.41e3 #num req. signals from line and box paper
    __diffuseDracoContribution = 2.74e-3 #inverse squared events per MeV) I think. It's eq. 4.2 in Line and Box Paper.
    __numBins = 0
    __binWidth = 1
    __plottedEnergyRange = (0.3,120) #measured in MeV
    def __init__(self):
        #generate proper bin stuff based upon energy range and bin width
        self.__binList = [Bin(i,self.__binWidth,self.__plottedEnergyRange(0)+(i*self.__binWidth)) for i in range(0,120,3)]
        #inherited classes will take care of initializing the number number of expected and observed counts







