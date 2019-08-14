import math
import numpy as np
from DMSG import DMSG_Bin as Bin

# a bin with 0 counts for all counts
bin1 = Bin.Bin(2,10)
# a bin with some initial counts
bin2 = Bin.Bin(2,10,1,1,1)
# a bin that gets set with some counts after initializing
bin3 = Bin.Bin(2,10)
print(bin3.getSignalCount())
bin3.setSignalCount(4)
print(bin3.getSignalCount())
#get the bin's signal count without any random variation.
print(bin3.getSignalCountDirectly())