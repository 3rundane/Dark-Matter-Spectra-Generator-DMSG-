
import Generation

# Start is a generic method to create a single generation object for plotting any other purposes. It's purpose is to streamline
# the process of creating generation objects so that you do not have to think about calculating binwidth. It may actually be
# simpler to alter the generation object to include the functionality of start...
# Anyways, in order to call start 5 parameters are required:
# numBins is the number of bins the user desires
# startEnergy is the starting energy value of the first bin object
# endingEnergy is the ending energy of the last bin object
# numParticles is the number of particles in the DDM ensemble denoted captial N in the Line and Box Paper
# sqrtSNaught is the center of mass energy of the lightest particle in the DDM ensemble as defined by the Line and Box Paper in the table on
# page 8.
# Returns: a generation object who has initialized and created a list of bin objects with both a pre-defined background event count list and a DDM count list
# Author: Dane Gollero
# Date: 7/17/2019
def start(numBins,startEnergy,endEnergy,numParticles=33,sqrtSNaught=135):
    binWidth = (endEnergy - startEnergy)/numBins
    draco = Generation.Generation(numBins, startEnergy, binWidth, numParticles, 2410, sqrtSNaught, 1, 2, 1) #TEST MODULE
    return draco

