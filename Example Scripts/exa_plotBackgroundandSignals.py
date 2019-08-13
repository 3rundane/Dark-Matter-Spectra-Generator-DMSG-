import numpy as np
import matplotlib.pyplot as plt
import math
from pylab import *
import matplotlib.pyplot as plt
import scipy.optimize as optimize
from DMSG import DMSG_Generation as Generation
from DMSG import DMSG_Start as Start
from DMSG import DMSG_Bin as Bin
# This file is a script used to compute plots of generated background spectra, generated mock signal events
# and background spectra + generated mock signal events.
# This file is organized as follows: first we use the Start "wrapper" to create a new Generation object for plotting uses.
# Second, we plot 3 different figures. Figure(1) plots in loglog scale the background events created by DracoBackground.py
# Figure(2) plots the DDM signal events generated with xscale in log form.
# Figure(3) plots the bin by bin sum of DDM signal events with the background events with xscale in log form.
draco=Start.start(48,10,120,48,135)
lengthArray = len(draco.getSignalBinCount())

lengthBG = len(draco.getBackgroundBinCount())

energyLength = len(draco.getBinEnergy())

totalcountlist = [0]*len(draco.getSignalBinCount())
for i in range(lengthArray):
    totalcountlist[i] = draco.getBackgroundBinCount()[i] + draco.getSignalBinCount()[i] #This has not been unittested yet.

plt.figure(1)
#plt.plot(draco.getBinEnergy(),draco.getBackgroundBinCount(),'g^')
plt.loglog(draco.getBinEnergy(),draco.getBackgroundBinCount(),'g.')
#plt.loglog(draco.getBinEnergy(),draco.getSignalBinCount(),'g^')
# plt.plot(draco.getBinEnergy(),draco.getSignalBinCount(),'r^')
plt.title('Draco Background Events')
plt.ylabel('total Events per bin width (#events/MeV)')
plt.xlabel('Bin energy(MeV)')




figure(2)
plt.plot(draco.getBinEnergy(),draco.getSignalBinCount(),'b.')
plt.title('DDM Signal Events')
plt.ylabel('total Events per bin width (#events/MeV)')
plt.xlabel('Bin energy(MeV)')
plt.xscale('log')

figure(3)
plt.plot(draco.getBinEnergy(),totalcountlist,'k.')
plt.title('Total Events')
plt.ylabel('total Events per bin width (#events/MeV)')
plt.xlabel('Bin energy(MeV)')
plt.xscale('log')
plt.show()