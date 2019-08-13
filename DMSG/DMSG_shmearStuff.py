import numpy as np
import matplotlib.pyplot as plt
import math as mt
import scipy.special as scp
import scipy.integrate as int

# This file contains a list of smearing functions for fun (or scienctific use).

# This function returns a value from a gaussian defined by the arguments. Note that this was intended to be the gaussian smearing
# function in the Lines and Boxes paper (equation 3.22).
# This function takes in the following parameters:
# eGamma - This is the energy whose probability of registering is given by the function.
# eGammaPrime - This is the actual energy of the photon coming in.
# epsilon - This is epsilon in the equation.
def gaussianSmearingFunction(eGamma, eGammaPrime,epsilon):
    return ((1.0/mt.sqrt(2*mt.pi)*epsilon*eGammaPrime))*(mt.exp(-1*((eGamma-eGammaPrime)**2/(2*(epsilon*eGammaPrime)**2))))

# This function returns a value from a gaussian defined by the arguments. Note that this is not the gaussian smearing
# function of the Lines and Boxes paper. The reason I did this was because I could not utilize the above function correctly.
# This function takes in the following parameters:
# eGamma - This is the energy whose probability of registering is given by the function.
# epsilon - This is epsilon in equation 3.22 of Lines and Boxes.
def gaussianSmearingFunction2(eGamma,epsilon):
    return ((1.0/mt.sqrt(2*mt.pi)*epsilon))*(mt.exp(-1*((eGamma)**2/(2*(epsilon)**2))))

