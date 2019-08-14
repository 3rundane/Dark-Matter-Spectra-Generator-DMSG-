import numpy as np
import matplotlib.pyplot as plt
import math
from pylab import *
import matplotlib.pyplot as plt
import scipy.optimize as optimize
from DMSG import DMSG_Generation as Generation
from DMSG import DMSG_Bin as Bin

numBins = 40
startEnergy = 10
binWidth = 2.75
numParticles = 23
sqrtSNaught = 135
draco = Generation.Generation(numBins, startEnergy, binWidth, numParticles, 2410, sqrtSNaught, 1, 2, 1)
binEnergy = draco.getBinEnergy()
binSignalCount = draco.getSignalBinCount()
binBackgroundCount = draco.getBackgroundBinCount()
print(binEnergy)
print(binSignalCount)
print(binBackgroundCount)
