import matplotlib.pyplot as plt
import numpy as np
from lmfit import Model,Parameters,Minimizer,minimize,report_fit
from DMSG import DMSG_Generation,DMSG_Bin,DMSG_DracoBackground,DMSG_DDMGenerationLoops
from DMSG import DMSG_computePoissonUncertainty as sig
import math
from TemplateFunctions import temp_templateFunctions as pfile
import csv

import scipy.optimize as opt


from collections.abc import Iterable

from DMSG import DMSG_Start as Start
# This file contains duplicates of the functions in vectorizedTemplates.py as well as an example script that fits
# DMSG-generated data via lmfit.

# def diffP(energy, sNaught, eta, deltaSSquare, sN, psi):
#     sqrtStar = math.sqrt(energy ** 2 + 135 ** 2) + energy
#     pieceOne = ((psi / 3) * ((sqrtStar) ** eta))
#     pieceTwo = (2 * (sqrtStar ** 2)) / (sqrtStar ** 2 + 135 ** 2)
#     pieceThree = np.heaviside((sqrtStar - math.sqrt(sNaught)), 0.5) * np.heaviside((math.sqrt(sN) - sqrtStar), 0.5)
#     # pieceThree = np.heaviside((math.sqrt(sNaught)), 0.5) * np.heaviside((sqrtStar), 0.5)
#     # print(pieceOne*pieceTwo*pieceThree)
#     if isinstance(energy, Iterable):
#         return np.array([diff1(energy_i,sNaught,eta,deltaSSquare,sN,psi) + diff2(energy_i,sNaught,eta,deltaSSquare,sN,psi) for energy_i in energy])
#         #return [diff1(energy_i,sNaught,eta,deltaSSquare,sN) for energy_i in energy]
#
#     else:
#         return np.array(diff1(energy,sNaught,eta,deltaSSquare,sN,phiNaught) + diff2(energy,sNaught,eta,deltaSSquare,sN,phiNaught)
#     return pieceOne * pieceTwo * pieceThree


def diffS(energy, sNaught, eta, deltaSSquare, sN, psi):
    sqrtMin = min(math.sqrt(sN), max(math.sqrt(sNaught), 2 * energy, (135 ** 2) / (2 * energy)))
    zOne = (135 ** 2) / (sqrtMin ** 2)
    zTwo = (135 ** 2) / sN
    betaOne = float(str(pfile.incompleteBetaFunctionExtension(-1 * eta / 2, 0, zOne).real))
    betaTwo = float(str(pfile.incompleteBetaFunctionExtension(-1 * eta / 2, 0, zTwo).real))
    piece1 = ((2 * psi) / (3))
    piece2 = (135) ** eta
    yourSecondDifferentialSir = piece1 * piece2 * (betaOne - betaTwo)
    return yourSecondDifferentialSir
def diff(energy,sNaught,eta,deltaSSquare,sN,phiNaught):
    def diff1(energy,sNaught,eta,deltaSSquare,sN,phiNaught):
        sqrtStar = math.sqrt(energy ** 2 + 135 ** 2) + energy
        pieceOne = ((phiNaught/ (3 * deltaSSquare))) * ((sqrtStar / math.sqrt(sNaught)) ** eta)
        pieceTwo = (2 * (sqrtStar ** 2)) / (sqrtStar ** 2 + 135 ** 2)
        pieceThree = np.heaviside((sqrtStar - math.sqrt(sNaught)), 0.5) * np.heaviside((math.sqrt(sN) - sqrtStar), 0.5)
        # pieceThree = np.heaviside((math.sqrt(sNaught)), 0.5) * np.heaviside((sqrtStar), 0.5)
        #print(pieceOne*pieceTwo*pieceThree)
        return pieceOne*pieceTwo*pieceThree
    def diff2(energy,sNaught,eta,deltaSSquare,sN,phiNaught):
        # print(type(min(math.sqrt(sN), max(math.sqrt(sNaught), 2 * energy, (135 ** 2) / (2 * energy)))))
        sqrtMin = min(math.sqrt(sN), max(math.sqrt(sNaught), 2 * energy, (135 ** 2) / (2 * energy)))
        zOne = (135 ** 2) / (sqrtMin ** 2)
        zTwo = (135 ** 2) / sN
        betaOne = float(str(pfile.incompleteBetaFunctionExtension(-1 * eta / 2, 0, zOne).real))
        betaTwo = float(str(pfile.incompleteBetaFunctionExtension(-1 * eta / 2, 0, zTwo).real))
        # print(type(betaOne))
        # print(betaOne)
        # print(type(betaTwo))
        piece1 = ((2 * phiNaught) / (3 * deltaSSquare))
        piece2 = (135 / math.sqrt(sNaught)) ** eta
        yourSecondDifferentialSir = piece1 * piece2*(betaOne-betaTwo)
        return yourSecondDifferentialSir
    if isinstance(energy, Iterable):
        return np.array([diff1(energy_i,sNaught,eta,deltaSSquare,sN,phiNaught) + diff2(energy_i,sNaught,eta,deltaSSquare,sN,phiNaught) for energy_i in energy])
        #return [diff1(energy_i,sNaught,eta,deltaSSquare,sN) for energy_i in energy]

    else:
        return np.array(diff1(energy,sNaught,eta,deltaSSquare,sN,phiNaught) + diff2(energy,sNaught,eta,deltaSSquare,sN,phiNaught))
        #return diff1(energy,sNaught,eta,deltaSSquare,sN)

def diffmk2(energy,sNaught,eta,deltaSSquare,sN,psi):
    # print(energy)
    # print(sNaught)
    # print(eta)
    # print(deltaSSquare)
    # print(sN)
    # print(phiNaught)
    def diff1(energy,sNaught,eta,deltaSSquare,sN,psi):
        sqrtStar = math.sqrt(energy ** 2 + 135 ** 2) + energy
        pieceOne = ((psi/3) * ((sqrtStar) ** eta))
        pieceTwo = (2 * (sqrtStar ** 2)) / (sqrtStar ** 2 + 135 ** 2)
        pieceThree = np.heaviside((sqrtStar - math.sqrt(sNaught)), 0.5) * np.heaviside((math.sqrt(sN) - sqrtStar), 0.5)
        # pieceThree = np.heaviside((math.sqrt(sNaught)), 0.5) * np.heaviside((sqrtStar), 0.5)
        #print(pieceOne*pieceTwo*pieceThree)
        return pieceOne*pieceTwo*pieceThree
    def diff2(energy,sNaught,eta,deltaSSquare,sN,psi):
        # print(type(min(math.sqrt(sN), max(math.sqrt(sNaught), 2 * energy, (135 ** 2) / (2 * energy)))))
        sqrtMin = min(math.sqrt(sN), max(math.sqrt(sNaught), 2 * energy, (135 ** 2) / (2 * energy)))
        zOne = (135 ** 2) / (sqrtMin ** 2)
        zTwo = (135 ** 2) / sN
        betaOne = float(str(pfile.incompleteBetaFunctionExtension(-1 * eta / 2, 0, zOne).real))
        betaTwo = float(str(pfile.incompleteBetaFunctionExtension(-1 * eta / 2, 0, zTwo).real))
        # print(type(betaOne))
        # print(betaOne)
        # print(type(betaTwo))
        piece1 = ((2 * psi) / (3))
        piece2 = (135) ** eta
        yourSecondDifferentialSir = piece1 * piece2*(betaOne-betaTwo)
        return yourSecondDifferentialSir
    if isinstance(energy, Iterable):
        return np.array([diff1(energy_i,sNaught,eta,deltaSSquare,sN,psi) + diff2(energy_i,sNaught,eta,deltaSSquare,sN,psi) for energy_i in energy])
        #return [diff1(energy_i,sNaught,eta,deltaSSquare,sN) for energy_i in energy]

    else:
        return np.array(diff1(energy,sNaught,eta,deltaSSquare,sN,phiNaught) + diff2(energy,sNaught,eta,deltaSSquare,sN,phiNaught))
        #return diff1(energy,sNaught,eta,deltaSSquare,sN)
def diffResidual(params,energy,data,dataError):
    residual = np.array([])
    params = params.valuesdict()
    sNaught = params['sNaught']
    eta = params['eta']
    deltaSSquare = params['deltaSSquare']
    sN = params['sN']
    psi = params['psi']

    residual = np.divide(np.subtract(diffmk2(energy,sNaught,eta,deltaSSquare,sN,psi),data),dataError)
    return residual
def diffPResidualP(params,energy,data,dataError):
    residual = np.array([])
    params = params.valuesdict()
    sNaught = params['sNaught']
    eta = params['eta']
    deltaSSquare = params['deltaSSquare']
    sN = params['sN']
    psi = params['psi']

    residual = np.divide(np.subtract(diffP(energy,sNaught,eta,deltaSSquare,sN,psi),data),dataError)
    return residual
def diffSResidualS(params,energy,data,dataError):
    residual = np.array([])
    params = params.valuesdict()
    sNaught = params['sNaught']
    eta = params['eta']
    deltaSSquare = params['deltaSSquare']
    sN = params['sN']
    psi = params['psi']

    residual = np.divide(np.subtract(diffS(energy,sNaught,eta,deltaSSquare,sN,psi),data),dataError)
    return residual

# Begin model fitting script.

draco=Start.start(40, 10, 120,33,164)

# energy values to plot
energyValues = draco.getBinEnergy()
# DDM event counts to plot
signalCounts = draco.getSignalBinCount()
# background event counts to help obtain uncertainty
backgroundCounts = draco.getBackgroundBinCount()

#This is how I obtained the errorbars in the plots.
dataError = np.array(sig.computeSigma(np.add(signalCounts,backgroundCounts)))

# The commented out code below is an example of using the Model class to perform the fits.
# gmodel = Model(diff)
#
# # init_parvals = gmodel.make_params(sNaught=164**2,eta=1,deltaSSquare=2,sN=230**2,phiNaught=65)
#result = gmodel.fit(signalCounts, energy=energyValues, params=init_parvals, weights=dataError)
#
# print(result.fit_report())
# final = signalCounts + result.residual
# plt.plot(energyValues,result.init_fit,'k--')
# plt.plot(energyValues,result.best_fit,'r-')

# The code below is an example of using the minimize function to perform the fits.

# Create and initialize Parameters objects for benchmarks A-D.
init_parvalsA = Parameters()
init_parvalsA.add_many(('sNaught',135**2,True,135**2,136**2,),
                      ('eta',1,True,-1,1),
                      ('deltaSSquare',1.8,True,0.5,2.1),
                      ('sN',181**2,True,180**2,182**2),
                      ('psi',6,True,2,8)) #A
init_parvalsB = Parameters()
init_parvalsB.add_many(('sNaught',135**2,True,135**2,136**2,),
                      ('eta',1,True,-1,1),
                      ('deltaSSquare',2,False,0,2.1),
                      ('sN',231**2,True,230**2,232**2),
                      ('psi',10,True,2,20)) #B
init_parvalsC = Parameters()
init_parvalsC.add_many(('sNaught',164**2,True,163**2,165**2,),
                      ('eta',1,True,-1,1),
                      ('deltaSSquare',2,False,0,2.1),
                      ('sN',180**2,True,179**2,181**2),
                      ('psi',6,True,2,8)) #C
init_parvalsD = Parameters()
init_parvalsD.add_many(('sNaught',164**2,True,163**2,165**2,),
                      ('eta',1,True,-1,1),
                      ('deltaSSquare',2,False,0,2.1),
                      ('sN',230**2,True,229**2,232**2),
                      ('psi',10,True,2,15)) #D

# Perform the minimization of the parameters with respect to the residual functions defined above.
result = minimize(diffResidual,init_parvalsD,method='least_squares',args=(energyValues,signalCounts,dataError))

print(report_fit(result))

# Get the best fit parameter values
updatedParams = result.params

# Fill out a list according to the best fit parameter values.
final = diffmk2(energyValues,updatedParams['sNaught'],updatedParams['eta'],updatedParams['deltaSSquare'],updatedParams['sN'],
              updatedParams['psi'])

# plot everything and format everything.
plt.errorbar(energyValues, signalCounts, yerr=dataError, fmt='k.', label="Mock Data")
plt.plot(energyValues,final,'r-',label="Fitted Curve")
chisq = result.redchi
plt.figtext(0.15, 0.55, r'$\chi^2/dof$: %.2f' % chisq, fontweight="bold")
plt.figtext(0.15,0.65,r'$S_o^1/2$: %.2f' %np.sqrt(updatedParams['sNaught']),fontweight="bold")
plt.figtext(0.15,0.6,r'$S_n^1/2$: %.2f' %np.sqrt(updatedParams['sN']),fontweight="bold")
plt.figtext(0.15,0.7,r'$\xi$: %.2f' %updatedParams['eta'],fontweight="bold")
plt.legend()
plt.xlabel('Energy[MeV]')

plt.ylabel('Counts per Bin Width [MeV]^-1')
plt.title('My Photon Spectrum B')
plt.xscale('log')
plt.xticks([10,15,20,30,50,70,100],['10','15','20','30','50','70','100'])
plt.xlim(10,120)
plt.ylim(-1.5)

plt.show()