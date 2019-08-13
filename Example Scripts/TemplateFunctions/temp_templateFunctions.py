import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.special as scp
import scipy.integrate as int
import mpmath as mp


#hypergeometric analytic continuation of the incomplete euler beta function. For use with the below differential flux software.
def incompleteBetaFunctionExtension(a,b,x):
    return mp.betainc(a, b, 0, x, regularized=False)


#Computes primary differential flux contribution
def differentialPhotonFluxP(energyGamma,phiNaught, sNaught, deltaSSquare,eta, pionMass,sN):
    # #This function returns equation 3.13 for plotting purposes.
    # sqrtStar = math.sqrt(energyGamma**2 + pionMass**2) + energyGamma
    #
    # pieceOne = ((phiNaught/(3*deltaSSquare)))*((sqrtStar/math.sqrt(sNaught))**eta)
    # pieceTwo = (2*(sqrtStar**2))/(sqrtStar**2 + pionMass**2)
    # pieceThree = np.heaviside((sqrtStar - math.sqrt(sNaught)),0.5)*np.heaviside((math.sqrt(sN)-sqrtStar),0.5)
    # return pieceOne*pieceTwo*pieceThree
    # This function returns equation 3.13 for plotting purposes.
    sqrtStar = math.sqrt(energyGamma ** 2 + pionMass ** 2) + energyGamma

    pieceOne = ((phiNaught / (3 * deltaSSquare))) * ((sqrtStar / math.sqrt(sNaught)) ** eta)
    pieceTwo = (2 * (sqrtStar ** 2)) / (sqrtStar ** 2 + pionMass ** 2)
    pieceThree = np.heaviside((sqrtStar - math.sqrt(sNaught)), 0.5) * np.heaviside((math.sqrt(sN) - sqrtStar), 0.5)
    return pieceOne * pieceTwo * pieceThree
#Computes secondary differential flux contribution
def differentialPhotonFluxS(energyGamma,phiNaught,sNaught,deltaSSquare,eta,pionMass,sN):
    sqrtMin = min(math.sqrt(sN),max(math.sqrt(sNaught),2*energyGamma,(pionMass**2)/(2*energyGamma)))
    zOne = (pionMass**2)/(sqrtMin**2)
    zTwo = (pionMass**2)/sN
    betaOne = incompleteBetaFunctionExtension(-1*eta/2,0,zOne)
    betaTwo = incompleteBetaFunctionExtension(-1*eta/2, 0, zTwo)
    piece1 = ((2*phiNaught)/(3*deltaSSquare))
    piece2 = (pionMass/math.sqrt(sNaught))**eta
    yourSecondDifferentialSir = piece1*piece2*(betaOne - betaTwo)
    return yourSecondDifferentialSir

print(differentialPhotonFluxP(2,1,135**2,2,1,135,181**2))
