import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.special as scp
import scipy.integrate as int
import mpmath as mp
import time
# This file demonstrates the use of the basic template functions derived by Boddy and Al. The home of the template functions
# are in this file, however there are other copies of the same ones in pieceFileMk2 and slightly modified copies to be used in a curve
# fitting analysis which are in modelFittingMark1.py.
#Author: Dane Gollero


# hypergeometric analytic continuation of the incomplete euler beta function. For use with the below differential flux software.
# This function takes in the follwing:
# a - This is some float negative or positive. I believe this function from mpmath also works with complex numbers but I would double check
# their documentation to be sure.
# b - This is again some floating point number.
# x - Same story as above.
# Note that the specific roles a, b, and x have in the context of the Lines and Boxes paper's equation 3.19 is specific and I've set up
# the argument order of this function to reflect that. (or I at least tried to lol).
def incompleteBetaFunctionExtension(a,b,z):
    return mp.betainc(a, b, 0, z, regularized=False)


# Computes primary differential flux contribution for a given energy point(energyGamma). This function returns equation 3.13 for plotting purposes.
# This function takes in:
# energyGamma - This is the incoming photon energy. Note that this function only works for energyGamma as a python scalar
# NOT as a list. I think this problem can be fixed by replacing math.sqrt() expressions with np.sqrt() but I ran out of time to
# try it. The vectorized version is contained in modelFittingMark1.
# phiNaught - This is the flux of the lightest ensemble particle. Or simply the normalization.
# sNaught - This is the square of the center of mass energy of the lightest ensemble particle.
# I accidentally switched conventions between this code and the data generation
# code by sometimes using the center of mass energy (sqrtSNaught) and here (other code uses where I use sNaught.
# deltaSSquare - This is the energy splitting first introduced in equation 2.7 of the Lines and Boxes paper.
# eta - The greek letter xi was used in the paper as the shape parameter which I call eta. Wherever you see eta I mean xi from the paper.
# pionMass - The mass of the neutral pion
# sN - The squared center of mass energy for the heaviest particle in the ensemble.
def differentialPhotonFluxP(energyGamma,phiNaught, sNaught, deltaSSquare,eta, pionMass,sN):
    # I am computing everything here by breaking equation 3.13 down into smaller pieces. I am doing this for the sake
    # of code readability.

    sqrtStar = math.sqrt(energyGamma**2 + pionMass**2) + energyGamma
    pieceOne = ((phiNaught/(3*deltaSSquare)))*((sqrtStar/math.sqrt(sNaught))**eta)
    pieceTwo = (2*(sqrtStar**2))/(sqrtStar**2 + pionMass**2)
    pieceThree = np.heaviside((sqrtStar - math.sqrt(sNaught)),0.5)*np.heaviside((math.sqrt(sN)-sqrtStar),0.5)

    return pieceOne*pieceTwo*pieceThree

#Computes secondary differential flux contribution for a given energy point. i.e it returns equation 3.19 for plotting purposes.
# The arguments here are the same as those above.
def differentialPhotonFluxS(energyGamma,phiNaught,sNaught,deltaSSquare,eta,pionMass,sN):
    # I am breaking equation 3.19 down into smaller, more manageable pieces.
    sqrtMin = min(math.sqrt(sN),max(math.sqrt(sNaught),2*energyGamma,(pionMass**2)/(2*energyGamma)))
    zOne = (pionMass**2)/(sqrtMin**2)
    zTwo = (pionMass**2)/sN
    betaOne = incompleteBetaFunctionExtension(-1*eta/2,0,zOne)
    betaTwo = incompleteBetaFunctionExtension(-1*eta/2, 0, zTwo)
    yourSecondDifferentialSir = (((2*phiNaught)/(3*deltaSSquare))*((pionMass/math.sqrt(sNaught))**eta)*(betaOne - betaTwo))
    return yourSecondDifferentialSir

# This function plots the derivative of the photon flux with respect to energy and plots it against energy. It is
# setup to do this for the benchmarks A-D given in the Lines and Boxes paper. It is very slow.
def computeBenchmarkFlux():
    #Initialize all relevant constants/ start code timer.
    start = time.time()
    energyGammaFull = np.arange(0.1,120,0.01)
    dim = len(energyGammaFull)
    pionMass = 135
    phiNaught = 1
    sNaught = np.array([18225, 18225, 26896, 26896])
    deltaSSquare = 2
    eta = 1
    sN = np.array([32761, 53361, 32400, 52900])# These are the squared center of mass energy values used in the benchmarks

    #Initialize some empty lists
    deltaFluxTA = np.zeros(dim,dtype=float)
    deltaFluxTB = np.zeros(dim, dtype=float)
    deltaFluxTC = np.zeros(dim, dtype=float)
    deltaFluxTD = np.zeros(dim, dtype=float)

    #Use template functions to fill out differential flux arrays.
    for i in range(0,dim):
        deltaFluxTA[i] = differentialPhotonFluxP(energyGammaFull[i],phiNaught,sNaught[0],deltaSSquare,eta,pionMass,sN[0]) +  differentialPhotonFluxS(energyGammaFull[i],phiNaught,sNaught[0],deltaSSquare,eta,pionMass,sN[0])
        deltaFluxTB[i] =differentialPhotonFluxP(energyGammaFull[i],phiNaught,sNaught[1],deltaSSquare,eta,pionMass,sN[1]) +  differentialPhotonFluxS(energyGammaFull[i],phiNaught,sNaught[1],deltaSSquare,eta,pionMass,sN[1])
        deltaFluxTC[i] =differentialPhotonFluxP(energyGammaFull[i],phiNaught,sNaught[2],deltaSSquare,eta,pionMass,sN[2]) +  differentialPhotonFluxS(energyGammaFull[i],phiNaught,sNaught[2],deltaSSquare,eta,pionMass,sN[2])
        deltaFluxTD[i] =differentialPhotonFluxP(energyGammaFull[i],phiNaught,sNaught[3],deltaSSquare,eta,pionMass,sN[3]) +  differentialPhotonFluxS(energyGammaFull[i],phiNaught,sNaught[3],deltaSSquare,eta,pionMass,sN[3])

    #Compute normalization constants
    n1 = int.trapz(deltaFluxTA,dx=0.01)
    n2 = int.trapz(deltaFluxTB,dx=0.01)
    n3 = int.simps(deltaFluxTC,energyGammaFull)
    n4 = int.simps(deltaFluxTD,energyGammaFull)

    #Plot everything
    plt.figure(1)
    plt.plot(energyGammaFull,deltaFluxTA/n1)
    plt.xscale('log')
    plt.ylim(0,0.075)
    plt.xlabel('Energy (MeV)')
    plt.ylabel('Normalized Differential Flux (1/MeV)')
    plt.title('Benchmark A')
    #
    plt.figure(2)
    plt.plot(energyGammaFull, deltaFluxTB/n2)
    plt.xscale('log')
    plt.ylim(0,0.05)
    plt.xlabel('Energy (MeV)')
    plt.ylabel('Normalized Differential Flux (1/MeV)')
    plt.title('Benchmark B')
    # # #
    plt.figure(3)
    plt.plot(energyGammaFull,deltaFluxTC)
    plt.title('Benchmark C')
    plt.xlabel('Energy (MeV)')
    plt.ylabel('Normalized Differential Flux (1/MeV)')
    plt.xscale('log')
    plt.ylim(0, 0.4)
    #
    plt.figure(4)
    plt.plot(energyGammaFull, deltaFluxTD)
    plt.xlabel('Energy (MeV)')
    plt.ylabel('Normalized Differential Flux (1/MeV)')
    plt.title('Benchmark D')
    plt.xscale('log')
    end = time.time()

    plt.show()


computeBenchmarkFlux()