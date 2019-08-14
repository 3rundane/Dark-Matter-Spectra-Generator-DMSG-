import numpy as np
import matplotlib.pyplot as plt
import math
from pylab import *
import matplotlib.pyplot as plt
import scipy.optimize as optimize
from DMSG import DMSG_Generation as Generation
from DMSG import DMSG_Start as Start
from DMSG import DMSG_Bin as Bin

draco = Generation.Generation(48,10,2.29,33,2410,135,1,2,1)
figure()
plt.plot(draco.getBinEnergy(),draco.getSignalBinCount(),'b.')
plt.title('DDM Signal Events')
plt.ylabel('total Events per bin width (#events/MeV)')
plt.xlabel('Bin energy(MeV)')
plt.xscale('log')
plt.show()