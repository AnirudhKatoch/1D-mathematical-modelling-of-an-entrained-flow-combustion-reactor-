# GasProperties

# Note
# Derived from the MATLAB class GasProperties written by  authors Netter, Tobias (tobias.netter@tum.de) ,
# Ceruti, Amedeo (amedeo.ceruti@tum.de) within the project "Entrained Flow Gasification Reactor Modeling with MATLAB"

# Class to get properties of a gas given temperature (K) and pressure (Pa) conditions.

# Example:
# oxygen = GasProperties('oxygen')

# Contents
#
# GasProperties
# mutable properties
# Protected proeprties
# Public methods
# GasProperties
# getMW
# computedensity
# computecp
# computebeta
# computelambda
# computeeta
# computePr
# computeH
# computeallproperties
# Private methods
# getcpconstants

import os
import numpy as np
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary # Using Refrop by setting up an refrop environment

class GasProperties:

    R = 8.314 # ideal gas constant. Units: J/(mol K)
    stepT = 10  # temperature interval steps for property extrapolation
    deltaT = 100 # total temperature range steps for property extrapolation
    name = None # Contains the name of the gas

# Construct an instance of GasProperties and checks the inputs to construct an instance of the class.
#
# INPUTS
# gasString: char array of a given support char array list: supportedGases =
#  {'methane', 'hydrogen', 'nitrogen', 'oxygen', 'water', 'CO', 'CO2', 'ash'}. Example: 'methane'.

    def __init__(self, gasString):
        supportedGases = ['methane', 'hydrogen', 'nitrogen', 'oxygen', 'water', 'CO', 'CO2', 'ash']
        if isinstance(gasString, str):
            if gasString in supportedGases:
                self.name = gasString
            else:
                raise ValueError(f'{gasString} is not supported. Please check documentation.')
        else:
            raise TypeError('gasString must be a string.')

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, gas):
        if isinstance(gas, str):
            self._name = gas
        else:
            raise TypeError('name must be a string.')

# getMW
# returns molecular weight in kg/mol of the given self.name
#
# INPUTS
# None
#
# OUTPUT
# MW: molecular weight in a double in kg/mol..

    def getMW(self):
        gas_name = self.name
        if gas_name == 'methane':
            MW = 16.04e-3  # [kg/mol]
        elif gas_name == 'hydrogen':
            MW = 2.01588e-3  # [kg/mol]
        elif gas_name == 'CO2':
            MW = 44.01e-3  # [kg/mol]
        elif gas_name == 'CO':
            MW = 28.01e-3  # [kg/mol]
        elif gas_name == 'water':
            MW = 18.01528e-3  # [kg/mol]
        elif gas_name == 'nitrogen':
            MW = 28.01348e-3
        elif gas_name == 'oxygen':
            MW = 31.998e-3  # [kg/mol]
        else:
            raise ValueError(f'Molecular weight for {gas_name} is not defined.')
        return MW

# getfluidname
# Get the fluid names for self.name that is used as input in  Refrop
#
# INPUTS
# None
#
# OUTPUT
# name: Fluid name according to Refrop directory

    def getfluidname(self):
        gas_name = self.name
        if gas_name == 'methane':
            name = "METHANE.FLD"
        elif gas_name == 'hydrogen':
            name = "HYDROGEN.FLD"
        elif gas_name == 'CO2':
            name = "CO2.FLD"
        elif gas_name == 'CO':
            name = "CO.FLD"
        elif gas_name == 'water':
            name = "WATER.FLD"
        elif gas_name == 'nitrogen':
            name = "NITROGEN.FLD"
        elif gas_name == 'oxygen':
            name = "OXYGEN.FLD"
        else:
            raise ValueError(f'Refrop name for {gas_name} is not defined.')
        return name

# computedensity
# Computes density of the fuel based on composition at given conditions.
#
# INPUTS
# T: Temperature in K as float
# P: Pressure in Pa as float
#
# OUTPUT
# rho: density in a float in [kg/m³].

    def computedensity(self, T, P):
        RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
        RP.SETPATHdll(os.environ['RPPREFIX'])
        r = RP.SETUPdll(1, self.getfluidname(), "HMX.BNC", "DEF")
        assert (r.ierr == 0)
        P = P / 1000
        x = RP.TPFLSHdll(T, P, [1.0])
        D = ((x.D) * self.getMW()) / 0.001
        return D

# computecp
# Computes heat capacity c_p at given conditions in [J/(kg*K)] .
#
# INPUTS
# T: Temperature in K as float
# P: Pressure in Pa as float
#
# OUTPUT
# c_p: double. heat capacity as float in units: [J/(kg*K)].

    def computecp(self, T, P):
        RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
        RP.SETPATHdll(os.environ['RPPREFIX'])
        r = RP.SETUPdll(1, self.getfluidname(), "HMX.BNC", "DEF")
        assert (r.ierr == 0)
        P = P / 1000
        x = RP.TPFLSHdll(T, P, [1.0])
        Cp = x.Cp/self.getMW()
        return Cp

# computebeta
# Computes Beta at given conditions in 1/K.
#
# INPUTS
# T: Temperature in K as float
#
# OUTPUT
# beta: float. Units: [1/K].

    def computebeta(self, T):
        beta = 1 / T
        return beta

# computelambda
# Computes Heat conductivity at given conditions in [W/(m*K)].
#
# INPUTS
# T: Temperature in K as double
# P: Pressure in Pa as double
#
# OUTPUT
# lambda: float. Heat conductivity in units: [W/(m*K)].

    def computelambda(self, T, P):
        RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
        RP.SETPATHdll(os.environ['RPPREFIX'])
        r = RP.SETUPdll(1, self.getfluidname(), "HMX.BNC", "DEF")
        assert (r.ierr == 0)
        P = P / 1000
        x = RP.TPFLSHdll(T, P, [1.0])
        y = RP.TRNPRPdll(T, x.D, [1.0])
        tcx = y.tcx
        return tcx

# computeeta
# Computes dynamic viscosity at given conditions in [Pa s].
#
# INPUTS
# T: Temperature in K as float
# P: Pressure in Pa as float
#
# OUTPUT
# eta: dynamic viscosity as float. Units: [Pa s].

    def computeeta(self, T, P):
        RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
        RP.SETPATHdll(os.environ['RPPREFIX'])
        r = RP.SETUPdll(1, self.getfluidname(), "HMX.BNC", "DEF")
        assert (r.ierr == 0)
        P = P / 1000
        x = RP.TPFLSHdll(T, P, [1.0])
        y = RP.TRNPRPdll(T, x.D, [1.0])
        eta = y.eta/1e6
        return eta

# computePr
# Computes Prandtl-number at given conditions in [-].
#
# INPUTS
# T: Temperature in K as float
# P: Pressure in Pa as float
#
# OUTPUT
# Pr: Prandtl-number as float. Units: Adimensional.

    def computePr(self, T, P):
        Pr = ((self.computeeta(T,P) * self.computecp(T, P)) / self.computelambda(T, P))
        return Pr

# computeH
# Computes Enthalpy at given conditions in [J/kg].
#
# INPUTS
# T: Temperature in K as float
# P: Pressure in Pa as float
#
# OUTPUT
# Pr: isobaric thermal expansion coefficient as float. Units: [J/kg].

    def computeH(self, T, P):
        RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
        RP.SETPATHdll(os.environ['RPPREFIX'])
        r = RP.SETUPdll(1, self.getfluidname(), "HMX.BNC", "DEF")
        assert (r.ierr == 0)
        P = P / 1000
        x = RP.TPFLSHdll(T, P, [1.0])
        h = x.h/self.getMW()
        return h

# computeallproperties
# Computes selected properties at given conditions.
#
# INPUTS
# T: Temperature in K as float
# P: Pressure in Pa as float
#
# OUTPUT
# [rho, eta, lambda, cp, Pr]: An array containing by order: [ rho [kg/m³], eta [Pa s], lambda [W/(m K)], cp, Pr [-] ]

    def computeallproperties(self, T, P):
        rho = self.computedensity(T, P)  # kg/m³ -- density
        eta = self.computeeta(T, P)  # Pas -- dynamic viscosity
        lambda_ = self.computelambda(T, P)  # W/(m K) -- thermal conductivity
        cp = self.computecp(T, P)  # J/(kg*K) -- isobaric coefficient of thermal expansion
        Pr = self.computePr(T, P)  # - -- Prandtl number
        return rho, eta, lambda_, cp, Pr

# Note
# Refrop has upper limit of 3000 K and 900,000,000 Pa , so if we go over that refrop will not work that is why
# one need other way to calculate the properties .
# If it ever happens one can calculate properties by normal formulas or approximate manually using fucntions such as
# polyval and polyfit

# getcpconstants
# returns constants to compute c_p manually.
#
# INPUTS
# None
#
# OUTPUT
# cts: array of doubles. Contains A,B,C,D,E,F,G in descending order.

    def getcpconstants(self):
        name = self.name

        if name == 'methane':
            A = 2.00500e3
            B = -6.81428e-1
            C = 7.08589e-3
            D = -4.71368e-6
            E = 8.5131e-10
            F = 0
            G = 0
            cts = [G, F, E, D, C, B, A]
        elif name == 'hydrogen':
            A = 1.4147e4
            B = 1.73720e-1
            C = 6.9e-4
            D = 0
            E = 0
            F = 0
            G = 0
            cts = [G, F, E, D, C, B, A]
        elif name == 'CO2':
            A = 2.00500e3
            B = -6.81428e-1
            C = 7.08589e-3
            D = -4.71368e-6
            E = 8.51317e-10
            F = 0
            G = 0
            cts = [G, F, E, D, C, B, A]
        elif name == 'CO':
            A = 1.04669e3
            B = -1.56841e-1
            C = 5.39904e-4
            D = -3.01061e-7
            E = 5.05048e-11
            F = 0
            G = 0
            cts = [G, F, E, D, C, B, A]
        elif name == 'water':
            A = 1.93780e3
            B = -1.18077
            C = 3.64357e-3
            D = -2.86327e-6
            E = 7.59578e-10
            F = 0
            G = 0
            cts = [G, F, E, D, C, B, A]
        elif name == 'nitrogen':
            A = 1.02705e3
            B = 2.16182e-2
            C = 1.48638e-4
            D = -4.48421e-8
            E = 0
            F = 0
            G = 0
            cts = [G, F, E, D, C, B, A]
        elif name == 'oxygen':
            A = 8.76317e2
            B = 1.22828e-1
            C = 5.58304e-4
            D = -1.20247e-6
            E = 1.14741e-9
            F = 5.12377e-13
            G = 8.56597e-17
            cts = [G, F, E, D, C, B, A]
        else:
            raise ValueError('Component not found')

        return cts
