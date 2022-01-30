import matplotlib.pyplot as plt
import numpy as np
from lmfit import Model,Parameters,Minimizer,minimize,report_fit
from DMSG import DMSG_Generation,DMSG_Bin,DMSG_DracoBackground,DMSG_DDMGenerationLoops
from DMSG import DMSG_computePoissonUncertainty as sig
import math
import temp_templateFunctions as pfile
import csv

import scipy.optimize as opt


from collections.abc import Iterable

import Start

# These are all of the relative functions that I created for plotting/fitting the data in modelFittingMark1.
# They are all used to either compute predictions or to be thrown into something from the lmfit package.

# Below, commented out, is a half-finished primary flux energy derivative function which attempts to implement a change in variables
# for the normalization constant (phi to psi).
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

# This is the "vectorized" version of the sum of differentialPhotonFluxP and differentialPhotonFluxS from pieceFileMk2.
# In other words it can accept energy as either a list or scalar. One thing to note is that it is not meant to be used
# directly in mizimization methods from the lmfit package. I threw this function into lmfit's Model class.
# This function takes in the following:
# energy - This is either a list or scalar describing the energies of equations 3.13 and 3.19
# sNaught - This is the squared center of mass energy of the lightest ensemble particle.
# eta - This is the shape parameter in equations 3.13 and 3.19. It's actually called xi but I didn't know that at the time.
# deltaSSquare - This is the energy splitting
# sN - This is the square of the heaviest ensemble particle's center of mass energy.
# phiNaught - This is the normalization coefficient.
def diff(energy,sNaught,eta,deltaSSquare,sN,phiNaught):
    def diff1(energy,sNaught,eta,deltaSSquare,sN,phiNaught):
        sqrtStar = math.sqrt(energy ** 2 + 135 ** 2) + energy
        pieceOne = ((phiNaught/ (3 * deltaSSquare))) * ((sqrtStar / math.sqrt(sNaught)) ** eta)
        pieceTwo = (2 * (sqrtStar ** 2)) / (sqrtStar ** 2 + 135 ** 2)
        pieceThree = np.heaviside((sqrtStar - math.sqrt(sNaught)), 0.5) * np.heaviside((math.sqrt(sN) - sqrtStar), 0.5)
        return pieceOne*pieceTwo*pieceThree
    def diff2(energy,sNaught,eta,deltaSSquare,sN,phiNaught):
        sqrtMin = min(math.sqrt(sN), max(math.sqrt(sNaught), 2 * energy, (135 ** 2) / (2 * energy)))
        zOne = (135 ** 2) / (sqrtMin ** 2)
        zTwo = (135 ** 2) / sN
        betaOne = float(str(pfile.incompleteBetaFunctionExtension(-1 * eta / 2, 0, zOne).real))
        betaTwo = float(str(pfile.incompleteBetaFunctionExtension(-1 * eta / 2, 0, zTwo).real))
        piece1 = ((2 * phiNaught) / (3 * deltaSSquare))
        piece2 = (135 / math.sqrt(sNaught)) ** eta
        yourSecondDifferentialSir = piece1 * piece2*(betaOne-betaTwo)
        return yourSecondDifferentialSir
    if isinstance(energy, Iterable):
        return np.array([diff1(energy_i,sNaught,eta,deltaSSquare,sN,phiNaught) + diff2(energy_i,sNaught,eta,deltaSSquare,sN,phiNaught) for energy_i in energy])

    else:
        return np.array(diff1(energy,sNaught,eta,deltaSSquare,sN,phiNaught) + diff2(energy,sNaught,eta,deltaSSquare,sN,phiNaught))
# This is the version of the above function which tries to implement some of the fitting ideas present in Lines and Boxes. Keep in mind
# it is also only for use with the Model class from lmfit.
# This function takes in the following:
# energy - This is either a list or scalar describing the energies of equations 3.13 and 3.19
# sNaught - This is the squared center of mass energy of the lightest ensemble particle.
# eta - This is the shape parameter in equations 3.13 and 3.19. It's actually called xi but I didn't know that at the time.
# deltaSSquare - This is the energy splitting
# sN - This is the square of the heaviest ensemble particle's center of mass energy.
# psi - This is an aggregate normalization coefficient defined by equation 5.5.
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
        sqrtMin = min(math.sqrt(sN), max(math.sqrt(sNaught), 2 * energy, (135 ** 2) / (2 * energy)))
        zOne = (135 ** 2) / (sqrtMin ** 2)
        zTwo = (135 ** 2) / sN
        betaOne = float(str(pfile.incompleteBetaFunctionExtension(-1 * eta / 2, 0, zOne).real))
        betaTwo = float(str(pfile.incompleteBetaFunctionExtension(-1 * eta / 2, 0, zTwo).real))

        piece1 = ((2 * psi) / (3))
        piece2 = (135) ** eta
        yourSecondDifferentialSir = piece1 * piece2*(betaOne-betaTwo)
        return yourSecondDifferentialSir
    if isinstance(energy, Iterable):
        return np.array([diff1(energy_i,sNaught,eta,deltaSSquare,sN,psi) + diff2(energy_i,sNaught,eta,deltaSSquare,sN,psi) for energy_i in energy])

    else:
        return np.array(diff1(energy,sNaught,eta,deltaSSquare,sN,phiNaught) + diff2(energy,sNaught,eta,deltaSSquare,sN,phiNaught))

# This function computes the total differential residual for plotting in lmfit minimization methods.
# This function takes in the following:
# params - This is an instantiated Parameters object from the lmfit package containing all model constraints to be optimized.
# energy - This is the list of the possible energy values you want to fit over. (has to be ordered)
# data - This is the mock spectral data the "DMSG" package generated.
# dataError - This is the error in the data that "DMSG" generated.
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
# This function computes the primary differential flux residual for plotting in lmfit minimization methods.
# This function takes in the following:
# params - This is an instantiated Parameters object from the lmfit package containing all model constraints to be optimized.
# energy - This is the list of the possible energy values you want to fit over. (has to be ordered)
# data - This is the mock spectral data the "DMSG" package generated.
# dataError - This is the error in the data that "DMSG" generated.
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
# This function computes the secondary differential flux residual for plotting in lmfit minimization methods.
# This function takes in the following:
# params - This is an instantiated Parameters object from the lmfit package containing all model constraints to be optimized.
# energy - This is the list of the possible energy values you want to fit over. (has to be ordered)
# data - This is the mock spectral data the "DMSG" package generated.
# dataError - This is the error in the data that "DMSG" generated.
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