# MixtureProperties

# Class to get properties of a mixture given temperature (K) and pressure (Pa) conditions.

# Example:
# Mixture = MixtureProperties([['nitrogen',0.81],['oxygen',0.18],['CO2',0.01]],components,298.15,101325,'Mixture')

# Content

# MixtureProperties
# Property setters
# set.xi(array)
# getMWs
# computeMW
# computemassfracs
# computecp
# computedensity
# computelambda
# computeeta
# computePr
# computebeta
# computeallproperties

from GasProperties import GasProperties
import numpy as np
import math

class MixtureProperties(GasProperties): # MixtureProperties is a subclass of class GasProperties
    T = None # temperature K
    P = None # pressure Pa
    xi = None # array of doubles with molar fractions
    nComps = None  # integer of the number of components
    gases = None # list of GasProperties objects
    components = None  # list of component names

# MixtureProperties
# Construct an instance of MixtureProperties
#
# INPUTS
# components: cell array of nComponents x 2, where the first item is a char array of the gas and the second one the molar frac.Example: 'nitrogen', 0.81; 'oxygen', 1 - 0.82; 'CO2', 0.01
# temperature: In K, float.Example: 298.15
# pressure: In Pascal , float.Example: 101325
# name: A different name can be given, original name is mixture

    def __init__(self, components, temperature, pressure, name="mixture"):
        args = [0]*4
        if isinstance(components, list):
            list_length = len(components)
            for dimension in components:
                if len(dimension) == 2:
                    args[0] = components
                else:
                    raise ValueError('molar Fractions needs to have 2 dimensions')
        else:
            raise TypeError('molarFractions is not a list')

        args[1] = temperature
        args[2] = pressure
        args[3] = name

        self.T = args[1]
        self.P = args[2]
        self.name = args[3]

        ncomp = len(args[0])  # num rows in the list
        gasesarray = [GasProperties(args[0][i][0]) for i in range(ncomp)]  # create list of GasProperties objects
        comps = [args[0][i][0] for i in range(ncomp)]  # create list of components
        x_i = [args[0][i][1] for i in range(ncomp)]  # create list of x_i values
        #
        self.gases = gasesarray
        self.components = comps
        self.xi = x_i
        self.nComps = ncomp

# Property setters

    @property
    def T(self):
        return self._T

    @T.setter
    def T(self, value):
        temp = value
        if temp > 0:
            self._T = temp
        else:
            raise ValueError('Temperature value must be positive')

    @property
    def P(self):
        return self._P

    @P.setter
    def P(self, value):
        pressure = value
        if pressure > 0:
            self._P = pressure
        else:
            raise ValueError('Pressure value must be positive')

# set.xi(array)

    def set_xi(self, array):
        TOL = 1e-4
        if abs(sum(array) - 1) < TOL:
            self.xi = array
        else:
            raise ValueError("Molar fractions need to be equal to 1. Sum xi = {}".format(sum(array)))

# getMWs
# returns an ordered vector with the MW of the mixture components.
#
# OUTPUT
# MW: double. Molecular weight in kg/mol of the mixture.

    def getMWs(self):
        MW = [0] * self.nComps
        for iComp in range(self.nComps):
            MW[iComp] = self.gases[iComp].getMW()
        return MW

# computeMW
# computes the molecular weight in kg/mol of the mixture by weighting of molar fractions.
#
# OUTPUT
# M: array of float. Molecular weight in kg/mol of the mixture.

    def computeMW(self):
        M = sum(x * y for x, y in zip(self.xi, self.getMWs()))
        return M

# computemassfracs
# Computes the mass fractions of the components and returns it in order
#
# INPUTS
# None
#
# OUTPUT
# wi: Ordered array of floats. mass fractions of the components in kg/kg.

    def computemassfracs(self):
        MW = self.getMWs()
        M = self.computeMW()
        wi = [xi * mw / M for xi, mw in zip(self.xi, MW)]
        return wi

# computecp
# Computes heat capacity of the mixture at given conditions in [J/(kg*K)] and weighing by mass fraction of the components.
#
# INPUTS
# None
#
# OUTPUT
# cp: float. Heat capacity in units: [J/(kg*K)].

    def computecp(self):
        icp = [0] * self.nComps
        for iComp in range(self.nComps):
            icp[iComp] = self.gases[iComp].computecp(self.T, self.P)
        cp = sum(x * y for x, y in zip(self.computemassfracs(), icp))
        return cp

# computedensity
# computes mixture density.
#
#Inputs
#None
#
# OUTPUT
# rho: float. density in kg/m3

    def computedensity(self):
        irho = [0] * self.nComps
        for iComp in range(self.nComps):
            irho[iComp] = self.gases[iComp].computedensity(self.T, self.P)
        rho = sum(x * y for x, y in zip(self.computemassfracs(), irho))
        return rho

# computelambda
# Computes Heat conductivity of the mixture at given conditions
#
# INPUTS
# None
#
# OUTPUT
# lambda: float. Heat conductivity units: [W/(m*K)].

    def computelambda(self):

        ilambda = [0] * self.nComps
        for iComp in range(self.nComps):
            ilambda[iComp] = self.gases[iComp].computelambda(self.T, self.P)

        # Note
        # Both methods can be used to calculate lambda but in method 2 there are errors which is need to be rectified
        # Can look into Wassiljewa, Alexandra. "Heat conduction in gas mixtures." Physics. Z. 5 (1904): 737-742.
        # Originally adapted from get_values_look_up_mixture_properties_fcn
        # Method 1 is more prone to error but it is less than 1 %

        # Method 1
        lambda_  = sum(x * y for x, y in zip(self.computemassfracs(), ilambda))

        # Method 2
        # MW = np.array(self.getMWs()).reshape((1, 3))
        # ilambda = np.array(ilambda).reshape((1, 3))
        #
        # Corr1 = 1 / (2 * math.sqrt(2)) * (np.power((1 + (np.matmul(np.transpose(MW), np.power(MW, -1)))), -0.5)) \
        #         * (np.power((1 + np.power(np.matmul(np.transpose(ilambda), np.power(ilambda, -1)), 0.5)) * \
        #                     (np.power(np.matmul(np.transpose(MW), np.power(MW, -1)), 0.25)), 2))
        #
        # lambda_ = np.sum(self.xi * ilambda*  np.power( np.transpose(np.sum((self.xi * Corr1), axis=1)),-1))

        return lambda_

# computeeta
# Computes dynamic viscosity of the mixture at given conditions
#
# INPUTS
# None
#
# OUTPUT
# eta: float. dynamic viscosity units: [Pa·s].

    def computeeta(self):
        ieta = [0] * self.nComps
        for iComp in range(self.nComps):
            ieta[iComp] = self.gases[iComp].computeeta(self.T, self.P)

        # Note
        # Both methods can be used to calculate lambda but in method 2 there are errors which is need to be rectified
        # Can look into Wassiljewa, Alexandra. "Heat conduction in gas mixtures." Physics. Z. 5 (1904): 737-742.
        # Originally adapted from get_values_look_up_mixture_properties_fcn
        # Method 1 is more prone to error but it is less than 1 %

        # Method 1
        eta = sum(x * y for x, y in zip(self.computemassfracs(), ieta))

        # Method 2
        # MW = np.array(self.getMWs()).reshape((1, 3))
        #
        # ieta = np.array(ieta).reshape((1, 3))
        #
        # Corr = 1 / (2 * math.sqrt(2)) * (np.power((1 + (np.matmul(np.transpose(MW), np.power(MW, -1)))), -0.5)) \
        #         * (np.power((1 + np.power(np.matmul(np.transpose(ieta), np.power(ieta, -1)), 0.5)) * \
        #                     (np.power(np.matmul(np.transpose(MW), np.power(MW, -1)), 0.25)), 2))
        #
        # eta = np.sum(self.xi * ieta * np.power(np.transpose(np.sum((self.xi * Corr), axis=1)), -1))

        return eta


# computePr
# Computes Prandtl number at given conditions in [-].
#
# INPUT
# None
#
# OUTPUT
# Pr: float. Prandtl number. Units: Adimensional.

    def computePr(self):
        Pr = self.computeeta() * self.computecp() / self.computelambda()
        return Pr

# computebeta
# Computes beta.
#
# INPUT
# None
#
# OUTPUT
# beta: float. beta to 1/K

    def computebeta(self):
        beta = 1 / self.T
        return beta

# computeallproperties
# Computes selected properties at given conditions.
#
# INPUT
# None
#
# OUTPUT
# [rho, eta, lambda, cp, Pr]: An array containing by order: [ rho [kg/m³], eta [Pa s], lambda [W/(m K)], cp, Pr [-] ]

    def computeallproperties(self):
        rho = self.computedensity()  # kg/m³ -- Density
        eta = self.computeeta()  # Pas -- dyne viscosity
        lambda_ = self.computelambda()  # W/(m K) -- thermal conductivity
        cp = self.computecp()  # 1/K -- isobaric thermal expansion coefficient
        Pr = self.computePr()  # - -- Prandtl number
        return rho, eta, lambda_, cp, Pr


