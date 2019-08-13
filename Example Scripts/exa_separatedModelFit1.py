import matplotlib.pyplot as plt
import numpy as np
from lmfit import Model,Parameters,Minimizer,minimize,report_fit
import Generation,Bin,DistributionsAndRules,DDMGenerationLoops
import computeSigma as sig
import math
import pieceFileMk2 as pfile
import csv

import scipy.optimize as opt


from collections.abc import Iterable

import Start

# This file contains duplicates of the functions in vectorizedTemplates.py, an example script that fits
# either the primary energy-derivative photon flux or the secondary energy-derivative photon flux (equations 3.13 and 3.19)
# to some DMSG-generated data via lmfit.

def diffP(energy, sNaught, eta, deltaSSquare, sN, psi):

    # pieceThree = np.heaviside((math.sqrt(sNaught)), 0.5) * np.heaviside((sqrtStar), 0.5)
    # print(pieceOne*pieceTwo*pieceThree)
    if isinstance(energy, Iterable):
        return np.array([diffP(energy_i,sNaught,eta,deltaSSquare,sN,psi) for energy_i in energy])
        #return [diff1(energy_i,sNaught,eta,deltaSSquare,sN) for energy_i in energy]

    else:
        sqrtStar = math.sqrt(energy ** 2 + 135 ** 2) + energy
        pieceOne = ((psi / 3) * ((sqrtStar) ** eta))
        pieceTwo = (2 * (sqrtStar ** 2)) / (sqrtStar ** 2 + 135 ** 2)
        pieceThree = np.heaviside((sqrtStar - math.sqrt(sNaught)), 0.5) * np.heaviside((math.sqrt(sN) - sqrtStar), 0.5)
        return pieceOne * pieceTwo * pieceThree
        #return diff1(energy,sNaught,eta,deltaSSquare,sN)




def diffS(energy, sNaught, eta, deltaSSquare, sN, psi):
    # print(type(min(math.sqrt(sN), max(math.sqrt(sNaught), 2 * energy, (135 ** 2) / (2 * energy)))))
    if isinstance(energy, Iterable):
        return np.array([diffS(energy_i, sNaught, eta, deltaSSquare, sN, psi) for energy_i in energy])
        # return [diff1(energy_i,sNaught,eta,deltaSSquare,sN) for energy_i in energy]

    else:
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
        yourSecondDifferentialSir = piece1 * piece2 * (betaOne - betaTwo)
        return yourSecondDifferentialSir
# def diff(energy,sNaught,eta,deltaSSquare,sN,phiNaught):
#     # print(energy)
#     # print(sNaught)
#     # print(eta)
#     # print(deltaSSquare)
#     # print(sN)
#     # print(phiNaught)
#     def diff1(energy,sNaught,eta,deltaSSquare,sN,phiNaught):
#         sqrtStar = math.sqrt(energy ** 2 + 135 ** 2) + energy
#         pieceOne = ((phiNaught/ (3 * deltaSSquare))) * ((sqrtStar / math.sqrt(sNaught)) ** eta)
#         pieceTwo = (2 * (sqrtStar ** 2)) / (sqrtStar ** 2 + 135 ** 2)
#         pieceThree = np.heaviside((sqrtStar - math.sqrt(sNaught)), 0.5) * np.heaviside((math.sqrt(sN) - sqrtStar), 0.5)
#         # pieceThree = np.heaviside((math.sqrt(sNaught)), 0.5) * np.heaviside((sqrtStar), 0.5)
#         #print(pieceOne*pieceTwo*pieceThree)
#         return pieceOne*pieceTwo*pieceThree
#     def diff2(energy,sNaught,eta,deltaSSquare,sN,phiNaught):
#         # print(type(min(math.sqrt(sN), max(math.sqrt(sNaught), 2 * energy, (135 ** 2) / (2 * energy)))))
#         sqrtMin = min(math.sqrt(sN), max(math.sqrt(sNaught), 2 * energy, (135 ** 2) / (2 * energy)))
#         zOne = (135 ** 2) / (sqrtMin ** 2)
#         zTwo = (135 ** 2) / sN
#         betaOne = float(str(pfile.incompleteBetaFunctionExtension(-1 * eta / 2, 0, zOne).real))
#         betaTwo = float(str(pfile.incompleteBetaFunctionExtension(-1 * eta / 2, 0, zTwo).real))
#         # print(type(betaOne))
#         # print(betaOne)
#         # print(type(betaTwo))
#         piece1 = ((2 * phiNaught) / (3 * deltaSSquare))
#         piece2 = (135 / math.sqrt(sNaught)) ** eta
#         yourSecondDifferentialSir = piece1 * piece2*(betaOne-betaTwo)
#         return yourSecondDifferentialSir
#     if isinstance(energy, Iterable):
#         return np.array([diff1(energy_i,sNaught,eta,deltaSSquare,sN,phiNaught) + diff2(energy_i,sNaught,eta,deltaSSquare,sN,phiNaught) for energy_i in energy])
#         #return [diff1(energy_i,sNaught,eta,deltaSSquare,sN) for energy_i in energy]
#
#     else:
#         return np.array(diff1(energy,sNaught,eta,deltaSSquare,sN,phiNaught) + diff2(energy,sNaught,eta,deltaSSquare,sN,phiNaught))
#         #return diff1(energy,sNaught,eta,deltaSSquare,sN)

def diffmk2(energy,sNaught,eta,deltaSSquare,sN,psi):

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
        return np.array(diff1(energy,sNaught,eta,deltaSSquare,sN,psi) + diff2(energy,sNaught,eta,deltaSSquare,sN,psi))
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
def diffResidualP(params,energy,data,dataError):
    residual = np.array([])
    params = params.valuesdict()
    sNaught = params['sNaught']
    eta = params['eta']
    deltaSSquare = params['deltaSSquare']
    sN = params['sN']
    psi = params['psi']

    residual = np.divide(np.subtract(diffP(energy,sNaught,eta,deltaSSquare,sN,psi),data),dataError)
    return residual
def diffResidualS(params,energy,data,dataError):
    residual = np.array([])
    params = params.valuesdict()
    sNaught = params['sNaught']
    eta = params['eta']
    deltaSSquare = params['deltaSSquare']
    sN = params['sN']
    psi = params['psi']

    residual = np.divide(np.subtract(diffS(energy,sNaught,eta,deltaSSquare,sN,psi),data),dataError)
    return residual


draco=Start.start(40, 10, 120,8,164)
energyValues = draco.getBinEnergy()
signalCounts = draco.getSignalBinCount()
backgroundCounts = draco.getBackgroundBinCount()
dataError = np.array(sig.computeSigma(np.add(signalCounts,backgroundCounts)))



#bound your parameters first
# gmodel = Model(diff)
#
# # init_parvals = gmodel.make_params(sNaught=164**2,eta=1,deltaSSquare=2,sN=230**2,phiNaught=65)
init_parvals = Parameters()
init_parvals.add_many(('sNaught',164**2,True,163**2,165**2,),
                      ('eta',1,True,-1,1),
                      ('deltaSSquare',2,False,0,2.1),
                      ('sN',181**2,True,180**2,182**2),
                      ('psi',20,True,5,30))

result = minimize(diffResidualS,init_parvals,method='least_squares',args=(energyValues,signalCounts,dataError))

print(report_fit(result))
#
# result = gmodel.fit(signalCounts,energy=energyValues,params=init_parvals,weights=dataError)
#
# print(result.fit_report())
# final = signalCounts + result.residual
updatedParams = result.params
print(updatedParams)
final = diffmk2(energyValues,updatedParams['sNaught'],updatedParams['eta'],updatedParams['deltaSSquare'],updatedParams['sN'],
              updatedParams['psi'])
plt.errorbar(energyValues, signalCounts, yerr=dataError, fmt='k.', label="data")

# plt.plot(energyValues,result.init_fit,'k--')
# plt.plot(energyValues,result.best_fit,'r-')
plt.plot(energyValues,final,'r-')
plt.xscale('log')
plt.ylim(-1.5)
plt.show()