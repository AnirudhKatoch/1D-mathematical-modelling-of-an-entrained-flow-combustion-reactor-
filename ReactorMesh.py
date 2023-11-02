# ReactorMesh

# Note
# Inspired from the MATLAB code written by  authors Netter, Tobias and 
# Ceruti, Amedeo within the project "Entrained Flow Gasification Reactor Modeling with MATLAB"

# Mesh class containing properties and methods to compute the main calculations of the simulation. It stores arrays of properties and variables which are updated and can be accessed.

# Example
# RM = ReactorMesh(Parameters, 2000, T); # construct reactor mesh object RM.X # X axis coordinate in each mesh step.

# Contents

# ReactorMesh
# Properties
# stored parameters
# stored arrays and matrices
# variables compacted into matrices
# ReactorMesh
# computeElementalN
# computeGasConcCantera
# computeAnnealing
# computeFuelSurface
# computePyrolysis
# computeIntrinsicReactivity
# computeGasPartT
# Section III: Eigenschaften Gas
# Section IV: Wärmeübergang Gas-Wand
# Section V: Wärmeübergang Partikel-Wand und Particel-Gas
# Section VI: Berechnung von Temperatur
# writeToCsv
# customInitialization
# initMeshResTime

import numpy as np
import cantera as ct
import pandas as pd
import math
import warnings
import csv
from GasProperties import GasProperties
from MixtureProperties import MixtureProperties

class ReactorMesh:

    Params = None  # GlobalParameters object

    # variables compacted into matrices
    X = None  # Reactor steps in m
    tau = None  # Gas input Residence Time in s
    # T = None  # Reactor temperature in K
    annFuel1 = None  # Annealing_Kohle1 [-]
    # intCO21 = None  # Intrinsic Reactivity CO2_Kohle1 [g/(m3s)]
    # intH2O1 = None  # Intrinsic Reactivity H2O_Kohle1 [g/(m3s)]
    surfaceFuel1 = None  # Surface_Kohle1 [m3/g]
    dFuel1 = None  # Particle Diameter_Kohle1 [10^-6 m]
    rhoFuel1 = None  # Particle Density_Kohle1 [kg/m2]
    # effCO2Fuel1 = None  # Effectiveness CO2_Kohle1 [-]
    # effH2OFuel1 = None  # Effectiveness H2O_Kohle1 [-]
    # effO2Fuel1 = None  # Effectiveness O2_Kohle1 [-]
    dXdt = None  # dX/dt [1/s]
    charConv = None  # Char Conversion post [-]
    overallConv = None  # Overall Conversion [-]
    xi = None  # [H2, CO, H2O, CO2, N2, O2] [-]
    nGas = None  # Total Gas amount [mol/s]
    realTau = None  # Real residence Time [s]
    dVYdtFuel1 = None  # dVY/dt --> Pyrolysis_Kohle1 [1/s]
    totalXC = None  # C-Conversion-Gesamt [-] # Conversion total
    # intO21 = None  # Intrinsic Reactivity
    # O2_Kohle1 [g/(m3s)]
    Tgas = None  # Gas-Temperature [K]
    Tpart1 = None  # Particle-Temperature_Kohle1 [K]

    # obsReactCO21 = None  # observed reactivity CO2_Kohle1 [g/(g*s)]
    # obsReactH2O1 = None  # observed reactivity H2O_Kohle1 [g/(g*s)]
    # obsReactO21 = None  # observed reactivity O2_Kohle1 [g/(g*s)]
    annCO21 = None  # Annealing CO2_Kohle1 [-]
    Tpart2 = None  # Particle-Temperature_Kohle2 [K]

    dVYdtfuel2 = None  # dVY/dt --> Pyrolysis_Kohle2 [1/s]
    totalXfuel1 = None  # Overall Conversion_Kohle1 [-]
    totalXfuel2 = None  # Overall Conversion_Kohle2 [-]
    annFuel2 = None  # Annealing_Kohle2 [-]
    # intCO22 = None  # Intrinsic Reactivity CO2_Kohle2 [g/(m3s)]
    # intH2O2 = None  # Intrinsic Reactivity H2O_Kohle2 [g/(m3s)]
    surfaceFuel2 = None  # Surface_Kohle2 [m3/g]
    dFuel2 = None  # Particle Diameter_Kohle2 [10^-6 m]
    rhoFuel2 = None  # Particle Density_Kohle2 [kg/m3]
    # effCO2Fuel2 = None  # Effectiveness CO2_Kohle2 [-]
    # effH2OFuel2 = None  # Effectiveness H2O_Kohle2 [-]
    annCO22 = None  # Annealing CO2_Kohle2 [-]
    # effO2Fuel2 = None  # Effectiveness O2_Kohle2 [-]
    # intO22 = None  # Intrinsic Reactivity O2_Kohle2 [g/(m3s)]
    charXFuel1 = None  # Char Conversion post_Kohle1 [-]
    charXFuel2 = None  # Char Conversion post_Kohle2 [-]
    dX1dt = None  # dX/dt_Kohle1 [1/s]
    dX2dt = None  # dX/dt_Kohle2 [1/s]
    # obsReactCO22 = None  # observed reactivity CO2_Kohle2 [g/(g*s)]
    # obsReactH2O2 = None  # observed reactivity H2O_Kohle2 [g/(g*s)]
    # obsReactO22 = None  # observed reactivity O2_Kohle2 [g/(g*s)]
    volatileYield1 = None  # VolatileYieldFinal Kohle1 [-]'
    volatileYield2 = None  # VolatileYieldFinal Kohle2 [-]'
    w1 = None  # Massenanteil Kohle1 [-]'
    w2 = None  # Massenanteil Kohle1 [-]'
    annFuncFuel1 = None  # Annealing Funktion Kohle 1 [-]'
    annFuncFuel2 = None  # Annealing Funktion Kohle 2 [-]'
    Fpc1 = None  # F_pc PoreClosing Funktion Kohle 1 [-]'
    Fpc2 = None  # F_pc_2 PoreClosing Funktion Kohle 2 [-]'
    heatTransfer = None  # Waermeuebergang W/m3 # Heat transfer W/m3
    pyrolysisHeat = None  # durch Pyrolyse freigesetzte Waerme J/m^3*s # heat released by pyrolysis J/m^3*s
    dYVdt = None  # dYV/dt_gesamt
    YV = None  # YV_gesamt   # YV_Total
    YV1 = None  # dYV/dt_1
    YV2 = None  # dYV/dt_2
    CX1 = None  # C-Conversion-Kohle1 # C conversion coal1
    CX2 = None  # C-Conversion-Kohle2
    surfaceDry1 = None  # Surface_dry_1
    surfaceDry2 = None  # Surface_dry_2

    # variables compacted into matrices

    obsReact1 = [None, None, None]  # [obsReactCO21, obsReactH2O1, obsReactO21]
    obsReact2 = [None, None, None]  # [obsReactCO22, obsReactH2O2, obsReactO22]
    int1 = [None, None, None]  # [intCO21, intH2O1, intO21]
    int2 = [None, None, None]  # [intCO22, intH2O2, intO22]
    effFuel1 = [None, None, None]  # [effCO2Fuel1, effH2OFuel1, effO2Fuel1]
    effFuel2 = [None, None, None]  # [effCO2Fuel2, effH2OFuel2, effO2Fuel2]
    u_sink_1_result = None  # added to track u sink of the gas particle temperature calculation
    u_sink_2_result = None  #

    steps = None

# ReactorMesh
#
# Construct an instance of this class
#
# INPUT
#
# Parameters: GlobalParameters object. Needs to be structured in order to access its parameter values.
# steps: positive integer. amount of discretization steps along the reactor X axis.
# T: float Initial gas temperature  in C

    def __init__(self, Parameters, steps, T):
        args = [0] * 1
        if steps == int(steps) and steps > 0:  # check if positive integer
            args[0] = steps
        else:
            raise ValueError('steps is not a positive integer')

        self.steps = steps

        # # These are the annealing property  of coal, they start as 1 but they reduce over time from 1 to 0
        setattr(self, 'Fpc1', np.ones(args[0]))
        setattr(self, 'Fpc2', np.ones(args[0]))
        setattr(self, 'annFuncFuel1', np.ones(args[0]))
        setattr(self, 'annFuncFuel2', np.ones(args[0]))

        # These properties will give us compositions of all six elements over the length
        setattr(self, 'xi', np.zeros((args[0], 6)))

        # These properties will give us these properties for CO2,H2O and O2 hence three columns
        setattr(self,'obsReact1', np.zeros((args[0], 3)))
        setattr(self,'obsReact2', np.zeros((args[0], 3)))
        setattr(self,'effFuel1', np.zeros((args[0], 3)))
        setattr(self,'effFuel2', np.zeros((args[0], 3)))
        setattr(self,'int1', np.zeros((args[0], 3)))
        setattr(self,'int2', np.zeros((args[0], 3)))

        #Properties not that are left
        setattr(self, 'X', np.zeros(args[0]))
        setattr(self, 'tau', np.zeros(args[0]))
        setattr(self, 'surfaceFuel1', np.zeros(args[0]))
        setattr(self, 'annFuel1', np.zeros(args[0]))
        setattr(self, 'dFuel1', np.zeros(args[0]))
        setattr(self, 'rhoFuel1', np.zeros(args[0]))
        setattr(self, 'dXdt', np.zeros(args[0]))
        setattr(self, 'charConv', np.zeros(args[0]))
        setattr(self, 'overallConv', np.zeros(args[0]))
        setattr(self, 'nGas', np.zeros(args[0]))
        setattr(self, 'realTau', np.zeros(args[0]))
        setattr(self, 'dVYdtFuel1', np.zeros(args[0]))
        setattr(self, 'totalXC', np.zeros(args[0]))
        setattr(self, 'Tgas', np.zeros(args[0]))
        setattr(self, 'Tpart1', np.zeros(args[0]))
        setattr(self, 'annCO21', np.zeros(args[0]))
        setattr(self, 'Tpart2', np.zeros(args[0]))
        setattr(self, 'dVYdtfuel2', np.zeros(args[0]))
        setattr(self, 'totalXfuel1', np.zeros(args[0]))
        setattr(self, 'totalXfuel2', np.zeros(args[0]))
        setattr(self, 'annFuel2', np.zeros(args[0]))
        setattr(self, 'surfaceFuel2', np.zeros(args[0]))
        setattr(self, 'dFuel2', np.zeros(args[0]))
        setattr(self, 'rhoFuel2', np.zeros(args[0]))
        setattr(self, 'annCO22', np.zeros(args[0]))
        setattr(self, 'charXFuel1', np.zeros(args[0]))
        setattr(self, 'charXFuel2', np.zeros(args[0]))
        setattr(self, 'dX1dt', np.zeros(args[0]))
        setattr(self, 'dX2dt', np.zeros(args[0]))
        setattr(self, 'volatileYield1', np.zeros(args[0]))
        setattr(self, 'volatileYield2', np.zeros(args[0]))
        setattr(self, 'w1', np.zeros(args[0]))
        setattr(self, 'w2', np.zeros(args[0]))
        setattr(self, 'heatTransfer', np.zeros(args[0]))
        setattr(self, 'pyrolysisHeat', np.zeros(args[0]))
        setattr(self, 'dYVdt', np.zeros(args[0]))
        setattr(self, 'YV', np.zeros(args[0]))
        setattr(self, 'YV1', np.zeros(args[0]))
        setattr(self, 'YV2', np.zeros(args[0]))
        setattr(self, 'CX1', np.zeros(args[0]))
        setattr(self, 'CX2', np.zeros(args[0]))
        setattr(self, 'surfaceDry1', np.zeros(args[0]))
        setattr(self, 'surfaceDry2', np.zeros(args[0]))
        setattr(self, 'u_sink_1_result', np.zeros(args[0]))
        setattr(self, 'u_sink_2_result', np.zeros(args[0]))

        self.Params = Parameters  # overwrite the Params property with GlobalParameter object
        self.customInitialization(T, steps)  # Initialize T-dependent props

# computeElementalN
# compute elemental input of the fuels and output it in mol/s in this order: [nC, nH, nO, nN].
#
# INPUT
# coalFeedRateWF: float. Coal feed rate without moisture and ash in kg/h.
#
# OUTPUT
# elems: array of floats Elemental input of Gas phase before Pyrolysis in order [C, H, O, N].

    def computeElementalN(self, coalFeedRateWF):

        nMoistureH = (self.Params.Conditions['coalFeedRate'] - coalFeedRateWF) / 0.018 * 2 / 3600  # mol/s
        nMoistureO = (self.Params.Conditions['coalFeedRate'] - coalFeedRateWF) / 0.018 * 1 / 3600  # mol/s
        nC = self.Params.n0[3]
        nO = self.Params.n0[1] * 2 + self.Params.n0[3] * 2 + self.Params.n0[4]
        nH = self.Params.n0[2] * 2 + self.Params.n0[4] * 2
        nN = self.Params.n0[0] * 2

        # Elemental input of Gas phase before Pyrolysis
        elems = [nC, nH + nMoistureH, nO + nMoistureO, nN]
        return nC, nH + nMoistureH, nO + nMoistureO, nN

# computeGasConcCantera
# calculates the gas concentration at T and P given the elemental inputs ni in mol/s. and updated obj.xi and obj.nGas
#
# INPUT
# T: float. Temperature in K
# P: float. Pressure in Pa
# ni: array of floats. elemental molar flows in mol/s in order: [C, H, O, N]. i: integer. Number of reaction step.
#
# OUTPUT
# nGas: float. total gas flow in mol/s.
# xi: array of doubles. molar fractions of components in order: [H2, CO, H2O, CO2, N2, O2]
# gas: Cantera Solution object.

    def computeGasConcCantera(self,T, P, ni):

        N_C = ni[0]
        N_H = ni[1]
        N_O = ni[2]
        N_N = ni[3]

        Total = sum(ni)

        x_C = N_C / Total
        x_H = N_H / Total
        x_O = N_O / Total
        x_N = N_N / Total

        C = 'C:' + str(x_C)
        H = 'H:' + str(x_H)
        O = 'O:' + str(x_O)
        N = 'N:' + str(x_N)

        BlankSpace = ' '

        # gas = importPhase('gri30.yaml')
        gas = ct.Solution('gri30.yaml')
        x = C + BlankSpace + H + BlankSpace + O + BlankSpace + N

        # I don't there is any need for multipying P by 1e5 as we are already in pixel
        gas.TPX = T, P , x
        gas.equilibrate('TP')

        Mol_share = gas.mole_fraction_dict()
        species_Names = gas.species_names
        species_Names = [species for species in species_Names if species != 'AR']

        #species_Name has a extra component AR which is not present in the Mol_share

        lines = len(Mol_share)

        Gas_composition = [['Species', 'Molar fraction']]

        for i in range(lines):
            Gas_composition.append([species_Names[i], Mol_share[species_Names[i]]])

       # Share of respective gases

        Gas_H2 = Gas_composition[1][1]
        Gas_CO = Gas_composition[15][1]
        Gas_H2O = Gas_composition[7][1]
        Gas_CO2 = Gas_composition[16][1]
        Gas_N2 = Gas_composition[48][1]
        Gas_O2 = Gas_composition[4][1]

        n_N2 = 0.5 * N_N
        nGas = n_N2 / Gas_N2  # with assumption that no nitrogen is involved in the reactions!
        xi = [Gas_H2, Gas_CO, Gas_H2O, Gas_CO2, Gas_N2, Gas_O2]  # [H2, CO, H2O, CO2, N2, O2]
        #
        return nGas, xi, gas

# computeAnnealing
# computes annealing of the particle at a Temperature T, a given discretization i and Particle properties.
#
# INPUT
#
# particleTemp: floats. Temperature of the particle in K in iteration i-1.
# particleProps: GlobalParameters.Fuel dictionary
# i: integer. discretization step (index).
#
# OUTPUT
#
# ann: float. Annealing value
# f_an: float. Annealing increment with respect of the previous annealing function value.

    def computeAnnealingStep(self, particleTemp, particleProps, i):
        deltaTau = self.realTau[i] - self.realTau[i - 1]
        tstep = 0.00001
        f_an = self.annFuncFuel1[i - 1]

        for _ in range(int(deltaTau / tstep)):
            d_s_T = -particleProps['k0_Annealing'] * math.exp(-particleProps['EA_Annealing'] * 1000 / (8.314 * particleTemp)) * f_an * tstep
            f_an += d_s_T

        ann = 1 + f_an * (particleProps['A_max_Annealing'] - 1)
        return ann, f_an

# computeFuelSurface
# computes annealing od the particle at a Temperature T, a given discretization i and Particle properties.
#
# INPUT
#
# particleTemp: double. Temperature of the particle in K in iteration i-1.
# previousFpc: double. Previous fuel surface value in iteration i-1.
# fuelParams: GlobalParameters.Fuel structure.
# i: integer. discretization step (index).
#
# OUTPUT
#
# F_pc: double. fuel surface value.

    def computeFuelSurface(self, particleTemp, previousFpc, fuelParams, i):

        # So realTau is the property hence it needs to be defined initially and can't be None'
        # There the initialization of properties of __init__ needs to be done

        deltat = self.realTau[i] - self.realTau[i-1]

        f_pc_pre = (previousFpc - fuelParams['A_min']) / (1 - fuelParams['A_min'])
        f_pc = f_pc_pre * (1 - fuelParams['k0_pc'] *
                           math.exp((-fuelParams['EA_pc']*1E3) / (8.314 * particleTemp)) * deltat)
        F_pc = fuelParams['A_min'] + f_pc * (1 - fuelParams['A_min'])
        return F_pc

# computePyrolysis
# computes several volatileyield parameters of the fuels for the pyrolysis tract of the reactor.
#
# INPUT
#
# Coal_Feed_Rate_waf_1: float. Coal dry feed rate in kg/h.
# Coal_Feed_Rate_waf_2: float.
# i: integer. discretization step (index).
#
# OUTPUT
# VolatileYieldFinal: floats
# VolatileYieldFinal_2: floats
# VolatileYield_1: floats
# VolatileYield_2: floats
# dYV_dt_1: floats
# dYV_dt_2: floats
# dH_pyr: floats
# X_pyr_kombiniert: floats

    def computePyrolysis(self, Coal_Feed_Rate_waf_1, Coal_Feed_Rate_waf_2, i):

        # Molar masses
        M_C = 12.011 # g / mol ~ kg / kmol
        M_H = 1.008
        M_O = 15.999

        # Standard enthalpy of formation
        H_0_CO2 = -3.935e8  # J/kmol
        H_0_H20 = -2.418e8  # J/kmol
        H_0_C = -101.268  # J/kmol

        # Evaporation enthalpy
        h_lat_H2O = 2.26e6  # J/kg

        VolatileYieldMax_1 = self.Params.Conditions['YVtotal']
        VolatileYieldMax_2 = self.Params.Conditions['YVtotal2']
        # Lower heating value of pure carbon
        LCV_C = -(H_0_CO2 - H_0_C) / M_C  # J/kg

        # Conversion of upper heating value (determined experimentally) to lower heating value (LCV) of coal 1 and coal 2.
        def computeLCV(fuel, VYmax):
            LCV = ((fuel['HCV'] - h_lat_H2O * fuel['moisture'] / 100) -
                   (fuel['elementalH'] / 100 * (2 * M_H + M_O)) /
                   (2 * M_H) * h_lat_H2O) / \
                  (1 - fuel['ashContent'] / 100 - fuel['moisture'] / 100)

            # Lower calorific value of volatiles of coal 1 and coal 2
            LCV_FL = (LCV - (1 - VYmax) * LCV_C) / VYmax  # J/kg
            return LCV_FL

        LCV_FL_1 = computeLCV(self.Params.Fuel1, VolatileYieldMax_1)
        LCV_FL_2 = computeLCV(self.Params.Fuel2, VolatileYieldMax_2)

        VolatileYieldFinal = self.Params.Conditions['YVtotal']
        VolatileYieldFinal_2 = self.Params.Conditions['YVtotal2']

        # Average molar mass of coal 1 and coal 2 liquids

        def computeMFL(fuel, VYMax):

            M_FL = VYMax / (((fuel['elementalC'] / 100) /
                             (1 - fuel['ashContent'] / 100 - fuel['moisture'] / 100)) -
                            (1 - VYMax)) * M_C  # kg/kmol
            return M_FL

        M_FL_1 = computeMFL(self.Params.Fuel1, VolatileYieldMax_1)
        M_FL_2 = computeMFL(self.Params.Fuel2, VolatileYieldMax_2)

        # Calculation of stochiometry coefficients for the formation of CO2 (Cx) and H20 (Hy)

        def computeCxHy(fuel, VYMax, MFL):

            Cx = ((fuel['elementalC'] / 100 / (1 - fuel['ashContent'] / 100 - fuel['moisture'] / 100)) -
                  (1 - VYMax)) / VYMax * MFL / M_C

            Hy = (fuel['elementalH'] / 100) / (
                        1 - fuel['ashContent'] / 100 - fuel['moisture'] / 100) / VYMax * MFL / M_H

            return Cx, Hy

        [Cx_1, Hy_1] = computeCxHy(self.Params.Fuel1, VolatileYieldMax_1, M_FL_1);
        [Cx_2, Hy_2] = computeCxHy(self.Params.Fuel2, VolatileYieldMax_2, M_FL_2);

        # Standard enthalpy of formation of the volatiles of coal 1 and coal 2

        H_0_FL_1 = LCV_FL_1 * M_FL_1 + (H_0_CO2 * Cx_1 + H_0_H20 * 0.5 * Hy_1)  # J/kmol
        H_0_FL_2 = LCV_FL_2 * M_FL_2 + (H_0_CO2 * Cx_2 + H_0_H20 * 0.5 * Hy_2)  # J/kmol

        V_R = self.Params.Geometry['innerRadius'] ** 2 * math.pi * (self.X[i] - self.X[i - 1])  # m^3

        def computedYVdt(fuel, T, VYfinal, dYV):

            def computeK(k, T, EA):
                a = k * math.exp(-(EA * 1E3 / 8.3145 / T))
                return a

            k_1 = computeK(fuel['k_01'], T, fuel['E_A1'])
            k_2 = computeK(fuel['k_02'], T, fuel['E_A2'])
            dYV_dt = (fuel['a_1'] * k_1 + fuel['a_2'] * k_2) * (VYfinal - dYV)

            return dYV_dt

        dYV_dt_1 = computedYVdt(self.Params.Fuel1, self.Tpart1[i - 1], VolatileYieldFinal, self.YV1[i - 1])
        dYV_dt_2 = computedYVdt(self.Params.Fuel2, self.Tpart2[i - 1], VolatileYieldFinal_2, self.YV2[i - 1])

        if dYV_dt_1 < 0:
            dYV_dt_1 = 0

        if dYV_dt_2 < 0:
            dYV_dt_2 = 0

        # Pyrolysis sales coal 1
        VolatileYield_1 = dYV_dt_1 * (self.realTau[i] - self.realTau[i - 1]) + self.YV1[i - 1]

        # Pyrolysis sales coal 2
        VolatileYield_2 = dYV_dt_2 * (self.realTau[i] - self.realTau[i - 1]) + self.YV2[i - 1]

        if VolatileYield_1 >= self.Params.Conditions['YVtotal']:
            VolatileYield_1 = self.Params.Conditions['YVtotal']
            dYV_dt_1 = (VolatileYield_1 - self.YV1[i - 1]) / (self.realTau[i] - self.realTau[i - 1])

        if VolatileYield_2 >= self.Params.Conditions['YVtotal2']:
            VolatileYield_2 = self.Params.Conditions['YVtotal2']
            dYV_dt_2 = (VolatileYield_2 - self.YV2[i - 1]) / (self.realTau[i] - self.realTau[i - 1])

        if (1 - self.Params.Conditions['w_2']) == 0:
            VolatileYield_1 = 0
            dYV_dt_1 = 0

            VolatileYield_2 = 0
            dYV_dt_2 = 0

        m_dot_FL_1_total = Coal_Feed_Rate_waf_1 / 3600 * VolatileYieldMax_1  # kg/s
        m_dot_FL_2_total = Coal_Feed_Rate_waf_2 / 3600 * VolatileYieldMax_2  # kg/s

        # mass flow of coal available for flux release
        H_pyr_K_1_total = (-H_0_FL_1) * m_dot_FL_1_total / M_FL_1  # J/s
        H_pyr_K_2_total = (-H_0_FL_2) * m_dot_FL_2_total / M_FL_2  # J/s
        dH_pyr = (dYV_dt_1 / VolatileYieldMax_1 * H_pyr_K_1_total +dYV_dt_2 / VolatileYieldMax_2 * H_pyr_K_2_total) * (self.realTau[i] - self.realTau[i - 1]) / V_R  # J/m^3*s

        # Overall Conversion(pyrolysis) combined af - based
        X_pyr_kombiniert = (VolatileYield_1 * Coal_Feed_Rate_waf_1 + VolatileYield_2 * Coal_Feed_Rate_waf_2) / (Coal_Feed_Rate_waf_1 + Coal_Feed_Rate_waf_2)

        return VolatileYieldFinal, VolatileYieldFinal_2, VolatileYield_1, VolatileYield_2, dYV_dt_1, dYV_dt_2, dH_pyr, X_pyr_kombiniert

# computeIntrinsicReactivity
# Computes the intrinsic reactivity of CO2, H2O and O2 in that order.
#
# INPUT
# partT: float. Particle temperature in K.
# partialPressures: array of floats. partial pressures of the components in this order: [CO2, H2O, O2].
# fuelProps: GlobalParameters.Fuel dictionary.
#
# OUTPUT
# Reactivity: array of floats. Reactivity of the components in this order [CO2, H2O, O2].

    def computeIntrinsicReactivity(self, partT, partialPressures, fuelProps):

        k0 = np.array(fuelProps["k0"])
        Ea_i = np.array(fuelProps["Ea_i"])
        n = np.array(fuelProps["n"])
        partialPressures = np.array(partialPressures)
        R = 8.3145  # Universal gas constant in J/(mol*K)

        Reactivity = (k0 * 1000) * np.exp(-Ea_i * 1000 / (R * partT)) * (partialPressures * 1E5) ** n

        return Reactivity

# computeGasPartT
# computes the Gas - Particle temperatures at reactor step i.
#
# INPUT
#
# P: float Pressure in Pa
# Coal_Feed_Rate_waf_1: float Coal feed rate in kg / s
# Coal_Feed_Rate_waf_2: float
# m_FC_C_1: float
# m_FC_C_2: float
# m_ash_1: float
# m_ash_2: float
# i: integer Reactor discretization step
# canteraOutput: Cantera solution object of the gas mixture.
# ParticleDevelopment: ParticleDevelopment object.
# ParticleDevelopment_2: ParticleDevelopment object.
#
# OUTPUT
#
# Temperatur_Gas_Austritt: float Gas temperature at i in K.
# Waermestromdichte: double Heat flow in W / m2
# Temperatur_Partikel_Austritt_1: float.Particle temperature in K.
# Temperatur_Partikel_Austritt_2: float.

    def computeGasPartT(self,P,Coal_Feed_Rate_waf_1,Coal_Feed_Rate_waf_2,m_FC_C_1,m_FC_C_2, m_ash_1,m_ash_2,i,canteraOutput,ParticleDevelopment,ParticleDevelopment_2):

        import warnings

        warnings.filterwarnings("ignore")  # Turn off the warning messages

        Temperatur_Partikel_1 = self.Tpart1[i - 1] - 273.15
        Temperatur_Partikel_2 = self.Tpart2[i - 1] - 273.15
        Temperatur_Gas_Eintritt = self.Tgas[i - 1] - 273.15  # [K]-->[°C]
        Temperatur_Rohr_aussen = self.Params.Conditions['outerTemperatureTube']
        Rohrlaenge = self.Params.Geometry['tubeLength']
        Molenstrom_Gas = self.nGas[i - 1]
        Durchmesser_innen = 2 * self.Params.Geometry['innerRadius']
        Rohrlaenge_ein = self.X[i - 1]
        Rohrlaenge_aus = self.X[i]
        Char_Conversion_Kohle1 = self.charXFuel1[i - 1]
        Char_Conversion_Kohle2 = self.charXFuel2[i - 1]
        Verweilzeit_ein = self.realTau[i - 1]
        Verweilzeit_aus = self.realTau[i]

        r_obs_CO2 = self.obsReact1[i - 1, 0]
        r_obs_H2O = self.obsReact1[i - 1, 1]
        r_obs_O2 = self.obsReact1[i - 1, 2]
        r_obs_CO2_2 = self.obsReact2[i - 1, 0]
        r_obs_H2O_2 = self.obsReact2[i - 1, 1]
        r_obs_O2_2 = self.obsReact2[i - 1, 2]

        dH_pyr = self.pyrolysisHeat[i - 1]
        # [H2, CO, H2O, CO2, N2, O2]
        x_O2 = self.xi[i - 1, 5]
        x_N2 = self.xi[i - 1, 4]
        x_H2O = self.xi[i - 1, 2]
        x_H2 = self.xi[i - 1, 0]
        x_CO = self.xi[i - 1, 1]
        x_CO2 = self.xi[i - 1, 3]

        gas = canteraOutput

        rho_Partikel_Mittel_1 = self.rhoFuel1[i - 1]
        rho_Partikel_Mittel_2 = self.rhoFuel2[i - 1]

        Overall_Conversion_1 = self.totalXfuel1[i - 1]
        Overall_Conversion_2 = self.totalXfuel2[i - 1]

        g = 9.81  # m/s² Acceleration due to gravity
        Rohrquerschnittsflaeche = Durchmesser_innen ** 2 * math.pi / 4  # m²
        Re_krit1 = 1000  # SDY       2300  % critical Reynolds number lam

        # Section II: Heat generation by reaction
        M_O2 = 2 * 15.999
        M_N2 = 28.0134
        M_H2O = 18.01528
        M_H2 = 2.01588
        M_CO = 28.01
        M_CO2 = 44.01
        M = x_O2 * M_O2 + x_N2 * M_N2 + x_H2O * M_H2O + x_H2 * M_H2 + x_CO * M_CO + x_CO2 * M_CO2  # M_N2=28.01348;  # kg/kmol -- molare Masse Stickstoff
        Massenstrom_Gas = Molenstrom_Gas * M * 10 ** -3  # kg/s
        M_C = 12.0107  # kg/kmol -- molar mass of carbon

        H_R_O2_298K = -393.34 * 10 ** 3  # J/mol  Reaction enthalpies @ 298K, 1,013 bar; Lit. "Catalytic effects in heterogeneous combustion reactions of carbon with respect to thermal waste treatment".
        H_R_CO2_298K = 172.22 * 10 ** 3  # J/mol
        H_R_H2O_298K = 131.25 * 10 ** 3  # J/mol

        delta_h_molar_C = (get_values_look_up_ash_properties_fcn(1, Temperatur_Partikel_1, 1, obj.Params.Fuel1) * (Temperatur_Partikel_1 + 273.15) - get_values_look_up_ash_properties_fcn(1, 25, 1, obj.Params.Fuel1) * (25 + 273.15)) * M_C / 10 ** 3
        delta_h_molar_C_2 = (get_values_look_up_ash_properties_fcn(1, Temperatur_Partikel_2, 2, obj.Params.Fuel2) * (Temperatur_Partikel_2 + 273.15) - get_values_look_up_ash_properties_fcn(1, 25, 2, obj.Params.Fuel2) * ( 25 + 273.15)) * M_C / 10 ** 3

        O2 = GasProperties('oxygen')
        h_O2_std = O2.computeH(25 + 273.15, P )
        delta_h_O2 = O2.computeH(Temperatur_Partikel_1 + 273.15, P ) - h_O2_std
        delta_h_O2_2 = O2.computeH(Temperatur_Partikel_2 + 273.15, P ) - h_O2_std
        delta_h_molar_O2 = delta_h_O2 * M_O2 / 10 ** 3  # J/mol
        delta_h_molar_O2_2 = delta_h_O2_2 * M_O2 / 10 ** 3

        CO2 = GasProperties('CO2')
        h_CO2_std = CO2.computeH(25 + 273.15, P )
        delta_h_CO2 = CO2.computeH(Temperatur_Partikel_1 + 273.15, P ) - h_CO2_std
        delta_h_CO2_2 = CO2.computeH(Temperatur_Partikel_2 + 273.15, P ) - h_CO2_std
        delta_h_molar_CO2 = delta_h_CO2 * M_CO2 / 10 ** 3
        delta_h_molar_CO2_2 = delta_h_CO2_2 * M_CO2 / 10 ** 3

        CO = GasProperties('CO')
        h_CO_std = CO.computeH(25 + 273.15, P )
        delta_h_CO = CO.computeH(Temperatur_Partikel_1 + 273.15, P ) - h_CO_std
        delta_h_CO_2 = CO.computeH(Temperatur_Partikel_2 + 273.15, P ) - h_CO_std
        delta_h_molar_CO = delta_h_CO * M_CO / 10 ** 3
        delta_h_molar_CO_2 = delta_h_CO_2 * M_CO / 10 ** 3

        CO = GasProperties('CO')
        h_CO_std = CO.computeH(25 + 273.15, P)
        delta_h_CO = CO.computeH(Temperatur_Partikel_1 + 273.15, P) - h_CO_std
        delta_h_CO_2 = CO.computeH(Temperatur_Partikel_2 + 273.15, P) - h_CO_std
        delta_h_molar_CO = delta_h_CO * M_CO / 10 ** 3
        delta_h_molar_CO_2 = delta_h_CO_2 * M_CO / 10 ** 3

        H2O = GasProperties('water')
        h_H2O_std = H2O.computeH(25 + 273.15, P)
        delta_h_H2O = H2O.computeH(Temperatur_Partikel_1 + 273.15, P) - h_H2O_std
        delta_h_H2O_2 = H2O.computeH(Temperatur_Partikel_2 + 273.15, P) - h_H2O_std
        delta_h_molar_H2O = delta_h_H2O * M_H2O / 10 ** 3
        delta_h_molar_H2O_2 = delta_h_H2O_2 * M_H2O / 10 ** 3

        H2 = GasProperties('hydrogen')
        h_H2_std = H2.computeH(25 + 273.15, P)
        delta_h_H2 = H2.computeH(Temperatur_Partikel_1 + 273.15, P) - h_H2_std
        delta_h_H2_2 = H2.computeH(Temperatur_Partikel_2 + 273.15, P) - h_H2_std
        delta_h_molar_H2 = delta_h_H2 * M_H2 / 10 ** 3
        delta_h_molar_H2_2 = delta_h_H2_2 * M_H2 / 10 ** 3

        H_R_O2 = H_R_O2_298K + delta_h_molar_CO2 - delta_h_molar_C - delta_h_molar_O2  # J / mol
        H_R_CO2 = H_R_CO2_298K + 2 * delta_h_molar_CO - delta_h_molar_C - delta_h_molar_CO2
        H_R_H2O = H_R_H2O_298K + delta_h_molar_CO + delta_h_molar_H2 - delta_h_molar_C - delta_h_molar_H2O

        H_R_O2_2 = H_R_O2_298K + delta_h_molar_CO2_2 - delta_h_molar_C_2 - delta_h_molar_O2_2
        H_R_CO2_2 = H_R_CO2_298K + 2 * delta_h_molar_CO_2 - delta_h_molar_C_2 - delta_h_molar_CO2_2
        H_R_H2O_2 = H_R_H2O_298K + delta_h_molar_CO_2 + delta_h_molar_H2_2 - delta_h_molar_C_2 - delta_h_molar_H2O_2

        m_dot_C_k1 = (1 - Char_Conversion_Kohle1) * m_FC_C_1 / 3600  # kg/s  # Fuel 1
        m_C_k1 = m_dot_C_k1 * (Verweilzeit_aus - Verweilzeit_ein)  # kg
        m_dot_C_k2 = (1 - Char_Conversion_Kohle2) * m_FC_C_2 / 3600  # kg/s  # Fuel 2
        m_C_k2 = m_dot_C_k2 * (Verweilzeit_aus - Verweilzeit_ein)  # kg

        m_dot_Particle_1 = (Coal_Feed_Rate_waf_1 * (1 - Overall_Conversion_1) + m_ash_1) / 3600  # kg/s  # Fuel 1
        m_Particle_1 = m_dot_Particle_1 * (Verweilzeit_aus - Verweilzeit_ein)
        m_dot_Particle_2 = (Coal_Feed_Rate_waf_2 * (1 - Overall_Conversion_2) + m_ash_2) / 3600  # kg/s  # Fuel  2
        m_Particle_2 = m_dot_Particle_2 * (Verweilzeit_aus - Verweilzeit_ein)  # kg

        w_1 = m_Particle_1 / (m_Particle_1 + m_Particle_2)
        w_2 = m_Particle_2 / (m_Particle_1 + m_Particle_2)

        Waermeerzeugung_gas = dH_pyr

        # Conversion of reaction rates from [g/(g*m²)] to [mol/(m³*s)].

        V_R_k1 = Rohrquerschnittsflaeche * (Rohrlaenge_aus - Rohrlaenge_ein)  # m³
        r_obs_molar_O2_1 = r_obs_O2 * m_C_k1 / (M_C * V_R_k1) / 10 ** 3  # [mol/(m³*s)]
        r_obs_molar_CO2_1 = r_obs_CO2 * m_C_k1 / (M_C * V_R_k1) / 10 ** 3  # [mol/(m³*s)]
        r_obs_molar_H2O_1 = r_obs_H2O * m_C_k1 / (M_C * V_R_k1) / 10 ** 3  # [mol/(m³*s)]

        V_R_k2 = Rohrquerschnittsflaeche * (Rohrlaenge_aus - Rohrlaenge_ein)  # m³
        r_obs_molar_O2_2 = r_obs_O2_2 * m_C_k2 / (M_C * V_R_k2) / 10 ** 3  # [mol/(m³*s)]
        r_obs_molar_CO2_2 = r_obs_CO2_2 * m_C_k2 / (M_C * V_R_k2) / 10 ** 3  # [mol/(m³*s)]
        r_obs_molar_H2O_2 = r_obs_H2O_2 * m_C_k2 / (M_C * V_R_k2) / 10 ** 3  # [mol/(m³*s)]

        Waermeerzeugung_particle_1 = r_obs_molar_O2_1 * (-H_R_O2) + r_obs_molar_CO2_1 * (-H_R_CO2) + r_obs_molar_H2O_1 * (-H_R_H2O)
        Waermeerzeugung_particle_2 = r_obs_molar_O2_2 * (-H_R_O2_2) + r_obs_molar_CO2_2 * (-H_R_CO2_2) + r_obs_molar_H2O_2 * (-H_R_H2O_2)

        # end Section II

        # Section III: Properties gas

        Temperatur_Rohr_innen = Temperatur_Rohr_aussen # °C -- Temperatur Rohrwandinnenseite

        # for temperature change from gas and convective heat transfer gas-wall, convective heat transfer gas-particles
        T_bezug = Temperatur_Gas_Eintritt  # C
        comps = [['nitrogen', x_N2], ['oxygen',x_O2], ['CO',x_CO], ['water',x_H2O], ['CO2',x_CO2], ['hydrogen', x_H2]]
        mixture = MixtureProperties(comps, T_bezug, P, 'stream')

        rho_N2, eta_N2, _, c_p_N2, _ = mixture.computeallproperties()

        lambda_N2 = thermalConductivity(gas)  # bc
        Pr_N2 = eta_N2 * c_p_N2 / lambda_N2;

        mixtureRohr = MixtureProperties(comps, Temperatur_Rohr_innen, P, 'rohr')
        _, eta_N2_W, _, c_p_N2_W, _ = mixtureRohr.computeallproperties()

        Pr_N2_wand = eta_N2_W * c_p_N2_W / lambda_N2

        # for free convection with wall

        T_bezug = (Temperatur_Gas_Eintritt + Temperatur_Rohr_innen) / 2  # °C  -- Reference temperature for substance properties
        mixture2 = MixtureProperties(comps, T_bezug, P, 'stream2')
        lambda_N2_W = lambda_N2

        rho_N2_W, eta_N2_T_W, _, c_p_T_N2_W, _ = mixture2.computeallproperties()
        beta_W = mixture2.computebeta()
        Pr_N2_T_W = eta_N2_T_W * c_p_T_N2_W / lambda_N2

        # for free convection with particle 1
        if w_1 > 0:
            T_bezug = (Temperatur_Gas_Eintritt + Temperatur_Partikel_1) / 2  # °C  -- Reference temperature for substance properties
            mixtureP = MixtureProperties(comps, T_bezug, P, 'P1')
            rho_N2_P1 = mixtureP.computedensity()
            eta_N2_P1 = mixtureP.computeeta()
            cp_N2_P1 = mixtureP.computecp()
            lambda_N2_P1 = lambda_N2
            Pr_N2_T_P1 = eta_N2_P1 * cp_N2_P1 / lambda_N2_P1
            beta_P1 = mixtureP.computebeta()

            mixtureP1W = MixtureProperties(comps, Temperatur_Partikel_1, P, 'P1 wand')
            eta_N2_P1_W = mixtureP1W.computeeta()
            cp_N2_P1_W = mixtureP1W.computecp()
            lambda_N2_P1_W = lambda_N2
            Pr_N2_T_W_P1 = eta_N2_P1_W * cp_N2_P1_W / lambda_N2_P1_W

        # for free convection with particle 2

        if w_2 > 0:
            T_bezug = (Temperatur_Gas_Eintritt + Temperatur_Partikel_2) / 2  # °C  -- Reference temperature for substance properties
            mixtureP2 = MixtureProperties(comps, T_bezug, P, 'P2')
            rho_N2_P2 = mixtureP2.computedensity()
            eta_N2_P2 = mixtureP2.computeeta()
            cp_N2_P2 = mixtureP2.computecp()
            lambda_N2_P2 = lambda_N2
            Pr_N2_T_P2 = eta_N2_P2 * cp_N2_P2 / lambda_N2_P2
            beta_P2 = mixtureP2.computebeta()

            mixtureP2W = MixtureProperties(comps, Temperatur_Partikel_2, P, 'P2 wand')
            eta_N2_P2_W = mixtureP2W.computeeta()
            cp_N2_P2_W = mixtureP2W.computecp()
            lambda_N2_P2_W = lambda_N2
            Pr_N2_T_W_P2 = eta_N2_P2_W * cp_N2_P2_W / lambda_N2_P2_W

        sigma_s = 5.67 * 10 ** -8  # W/(m^2 K^4) -- Stefan-Boltzmann radiation constant # Thermal radiation
        emissionsgrad_N2 = 0.03  # Emissivity N2

        # end Section III

        # Section IV: Heat transfer gas-wall

        Stroemungsgeschwindigkeit = Massenstrom_Gas / (rho_N2 * Rohrquerschnittsflaeche)  # m/s  # Flow velocity
        Re = Stroemungsgeschwindigkeit * Durchmesser_innen * rho_N2 / eta_N2  # Reynolds number

        # convective gas wall (SSC 10.2 p.37)
        if Re <= Re_krit1:  # LAMINAR
            L_e_h = 0.056 * Re * Durchmesser_innen  # ca 9 m -- hydrodynamic inlet length
            # hydrodynamic and thermal start-up with laminar pipe flow
            Pe = Re * Pr_N2 * Durchmesser_innen / Rohrlaenge
            Nu_t = 3.657 / math.tanh(2.264 * (1 / Pe) ** (1 / 3) + 1.7 * (1 / Pe) ** (2 / 3)) + 0.0499 * Pe * math.tanh(1 / Pe)  # Baehr (WSÜ 10.2.1 S.38)
            Nu = Nu_t / math.tanh(2.432 * (Rohrlaenge / Durchmesser_innen / Re) ** (1 / 6))

        else:  # TURBULENT
            # Correction factors
            K_L = 1 + (Durchmesser_innen / Rohrlaenge) ** (2 / 3)  # -  -- Correction factor: thermal & hydrodynamic start-up
            K_Pr = (Pr_N2 / Pr_N2_wand) ** 0.11  # -  -- Correction factor: temperature-dependent substance values
            # Pipe friction coefficient (Petukhov)
            zeta = (0.790 * math.log(Re) - 1.64) ** (-2)
            # Nusseltzahl (Gnielinski)
            Pr = Pr_N2
            Nu = (zeta / 8 * (Re - 1000) * Pr) / (1 + 12.7 * math.sqrt(zeta / 8) * (Pr ** (2 / 3) - 1)) * K_L * K_Pr

        alpha_konv = Nu * lambda_N2 / Durchmesser_innen  # W/(m² K) -- convective heat transfer coefficient (WSÜ 10.2.2 p.40)

        # free convection gas wall
        Gr = g * beta_W * abs(Temperatur_Rohr_innen - Temperatur_Gas_Eintritt) * Durchmesser_innen ** 3 / ((eta_N2_W / rho_N2_W) ** 2)  # - -- Grashof-Zahl
        Ra = Gr * Pr_N2_T_W

        C1 = 1 / 16
        C2 = 0.52
        Nu_FK = (1 / (C1 * Ra) ** (3 / 2) + 1 / (C2 * Ra ** (1 / 4)) ** (3 / 2)) ** (-2 / 3)
        alpha_FK = Nu_FK * lambda_N2_W / Durchmesser_innen

        # Radiation wall gas
        alpha_Strahlung = emissionsgrad_N2 * sigma_s * ((Temperatur_Gas_Eintritt + 273.15) + (Temperatur_Rohr_innen + 273.15)) * ((Temperatur_Gas_Eintritt + 273.15) ** 2 + (Temperatur_Rohr_innen + 273.15) ** 2) # W / (m ^ 2 K)
        U_innen = alpha_konv + alpha_FK + alpha_Strahlung # radius_inside/lambda_ceramics*log(radius_outside/radius_inside))
        Waermestromdichte = U_innen * (Temperatur_Rohr_aussen - Temperatur_Gas_Eintritt) * 10 ** -3

        # end Section IV

        # Section V: Particle-Wall and Particle-Gas Heat Transfer

        epsillon_Partikel = 0.8  # Emissivity particles (Tremel)
        sigma_s = 5.67e-8  # W/(m^2 K^4) # Stefan-Boltzmann radiation constant


        c_C_1 = get_values_look_up_ash_properties_fcn(1, Temperatur_Partikel_1, 1)
        c_C_2 = get_values_look_up_ash_properties_fcn(1, Temperatur_Partikel_2, 2)
        c_Ash_1 = get_values_look_up_ash_properties_fcn(2, Temperatur_Partikel_1, 1, obj.Params.Fuel1)
        c_Ash_2 = get_values_look_up_ash_properties_fcn(2, Temperatur_Partikel_2, 2, obj.Params.Fuel2)

        if m_dot_C_k1 > 0:
            c_Partikel_1 = (c_C_1 * (m_dot_Particle_1 - m_ash_1 / 3600) + c_Ash_1 * m_ash_1 / 3600) / m_dot_Particle_1
        else:
            c_Partikel_1 = c_Ash_1

        if m_dot_C_k2 > 0:
            c_Partikel_2 = (c_C_2 * (m_dot_Particle_2 - m_ash_2 / 3600) + c_Ash_2 * m_ash_2 / 3600) / m_dot_Particle_2
        else:
            c_Partikel_2 = c_Ash_2

#------------------------------------------------------------------------------------------------------------------------

        # PARTICLE DISTRIBUTION @TODO

        d_Partikel_Mittel_1 = sum(ParticleDevelopment.D[:, i - 1] * 10 ** (-6) * ParticleDevelopment.xFrac[:, i - 1] / 100)
        Oberflaeche_1 = sum(math.pi * (ParticleDevelopment.D[:, i - 1] * 10 ** (-6)) ** 2 * ParticleDevelopment.xFrac[:,i - 1] / 100)  # *((1-ParticleDevelopment.(['GasificationStep' num2str(Reactionstep-1)]){j,34}));
        Masse_1 = sum(math.pi / 6 * (ParticleDevelopment.D[:, i - 1] * 10 ** (-6)) ** 3 * ParticleDevelopment.xFrac[:,i - 1] / 100 * ParticleDevelopment.rho[:, i - 1])

        d_Partikel_Mittel_2 = sum(ParticleDevelopment_2.D[:, i - 1] * 10 ** (-6) * ParticleDevelopment_2.xFrac[:, i - 1] / 100)
        Oberflaeche_2 = sum(math.pi * (ParticleDevelopment_2.D[:, i - 1] * 10 ** (-6)) ** 2 * ParticleDevelopment_2.xFrac[:,i - 1] / 100)  # *((1-ParticleDevelopment.(['GasificationStep' num2str(Reactionstep-1)]){j,34}));
        Masse_2 = sum(math.pi / 6 * (ParticleDevelopment_2.D[:, i - 1] * 10 ** (-6)) ** 3 * ParticleDevelopment_2.xFrac[:,i - 1] / 100 * ParticleDevelopment_2.rho[:, i - 1])

        S_m_1 = Oberflaeche_1 / Masse_1  # m²/kg --> specific surface of the particle cluster
        S_m_2 = Oberflaeche_2 / Masse_2  # m²/kg --> specific surface of the particle cluster

        Temperatur_Gas_Eintritt = Temperatur_Gas_Eintritt + 273.15  # in K
        Temperatur_Rohr_innen = Temperatur_Rohr_aussen + 273.15  # in K
        Temperatur_Partikel_Eintritt_1 = Temperatur_Partikel_1 + 273.15  # in K
        Temperatur_Partikel_Eintritt_2 = Temperatur_Partikel_2 + 273.15  # in K

#--------------------------------------------------------------------------------------------------------------------------
        # Radiation
        # Wall-->Particle
        alpha_str_WP_1 = epsillon_Partikel * sigma_s * (Temperatur_Partikel_Eintritt_1 + Temperatur_Rohr_innen) * ( Temperatur_Partikel_Eintritt_1 ** 2 + Temperatur_Rohr_innen ** 2)
        alpha_str_WP_2 = epsillon_Partikel * sigma_s * (Temperatur_Partikel_Eintritt_2 + Temperatur_Rohr_innen) * ( Temperatur_Partikel_Eintritt_2 ** 2 + Temperatur_Rohr_innen ** 2)
        # Gas<-->Particle
        alpha_str_GP_1 = epsillon_Partikel * sigma_s * (Temperatur_Partikel_Eintritt_1 + Temperatur_Gas_Eintritt) * (Temperatur_Partikel_Eintritt_1 ** 2 + Temperatur_Gas_Eintritt ** 2)
        alpha_str_GP_2 = epsillon_Partikel * sigma_s * (Temperatur_Partikel_Eintritt_2 + Temperatur_Gas_Eintritt) * (Temperatur_Partikel_Eintritt_2 ** 2 + Temperatur_Gas_Eintritt ** 2)

#-------------------------------------------------------------------------------------------------------------------------

        # Circulation of gas around particles - convective VT gas-particles

        if i == 2:
            u_sink_1_pre = 0.5  # SDY 6.2 Boundary conditions
            u_sink_2_pre = 0.5
        else:
            u_sink_1_pre = self.u_sink_1_result[i - 1]
            u_sink_2_pre = self.u_sink_2_result[i - 1]

        Delta_t = min(0.001, (Verweilzeit_aus - Verweilzeit_ein) / 2)  # Time step of the iteration

        if w_1 == 0: #Overall_Conversion_1 < 0.99
            u_sink_1 = Stroemungsgeschwindigkeit
            Re_1 = 0
        elif rho_Partikel_Mittel_1 < rho_N2:  # Prevents particle density smaller than gas density at high conversions
            u_sink_1 = Stroemungsgeschwindigkeit # if particle density is less than gas density particle has same velocity as gas
            Re_1 = 0
        else:
            u_sink_0_1 = u_sink_1_pre
            Re_1 = abs(u_sink_1_pre - Stroemungsgeschwindigkeit) * d_Partikel_Mittel_1 * rho_N2 / eta_N2
            for InterrationStep in range(int((Verweilzeit_aus - Verweilzeit_ein) / Delta_t)):
                if Re_1 <= 10 ** 3:
                    while True:
                        index = np.sign(Stroemungsgeschwindigkeit - u_sink_0_1)
                        Re_1 = abs(u_sink_0_1 - Stroemungsgeschwindigkeit) * d_Partikel_Mittel_1 * rho_N2 / eta_N2
                        zeta_1 = 24 / Re_1
                        a_1 = g - rho_N2 / rho_Partikel_Mittel_1 * (g - index * 3 / 4 * zeta_1 / d_Partikel_Mittel_1 * abs(u_sink_0_1 - Stroemungsgeschwindigkeit) ** 2)
                        u_sink_1_1 = u_sink_1_pre + a_1 * Delta_t
                        # Abort criterion
                        if abs(u_sink_1_1 - u_sink_0_1) / u_sink_0_1 * 100 <= 1: # (corresponds to 1% change)
                            break
                        u_sink_0_1 = u_sink_1_1
                    u_sink_1 = u_sink_1_1
                else:
                    while True:
                        index = np.sign(Stroemungsgeschwindigkeit - u_sink_0_1)
                        Re_1 = abs(u_sink_0_1 - Stroemungsgeschwindigkeit) * d_Partikel_Mittel_1 * rho_N2 / eta_N2
                        zeta_1 = 0.4
                        a_1 = g - rho_N2 / rho_Partikel_Mittel_1 * (g - index * 3 / 4 * zeta_1 / d_Partikel_Mittel_1 * abs(u_sink_0_1 - Stroemungsgeschwindigkeit) ** 2)
                        u_sink_1_1 = u_sink_1_pre + a_1 * (Verweilzeit_aus - Verweilzeit_ein)
                        if abs(u_sink_1_1 - u_sink_0_1) / u_sink_0_1 * 100 <= 1:  # (corresponds to 1% change)
                            break
                        u_sink_0_1 = u_sink_1_1
                    u_sink_1 = u_sink_1_1
        if Re_1 > 3.5:
            Nu_konv_K_1 = 2 + (0.4 * Re_1 ** 0.5 + 0.06 * Re_1 ** 0.067) * Pr_N2 ** 0.4 * (Pr_N2 / Pr_N2_T_W_P1) ** 0.25
            alpha_konv_K_1 = Nu_konv_K_1 * lambda_N2 / d_Partikel_Mittel_1
        else:
            alpha_konv_K_1 = 0
        # Flow around coal 2

        if w_2 == 0 : #Overall_Conversion_1 < 0.99
            u_sink_2 = Stroemungsgeschwindigkeit
            Re_2 = 0
        elif rho_Partikel_Mittel_2 < rho_N2: # Prevents particle density smaller than gas density at high conversions
            u_sink_2 = Stroemungsgeschwindigkeit # if particle density is less than gas density particle has same velocity as gas
            Re_2 = 0
        else:
            u_sink_0_2 = u_sink_2_pre
            Re_2 = abs(u_sink_2_pre - Stroemungsgeschwindigkeit) * d_Partikel_Mittel_2 * rho_N2 / eta_N2
            for InterrationStep in range(int((Verweilzeit_aus - Verweilzeit_ein) / Delta_t)):
                if Re_2 <= 10 ** 3:
                    while True:
                        index = np.sign(Stroemungsgeschwindigkeit - u_sink_0_2)
                        Re_2 = abs(u_sink_0_2 - Stroemungsgeschwindigkeit) * d_Partikel_Mittel_2 * rho_N2 / eta_N2
                        zeta_2 = 24 / Re_2
                        a_2 = g - rho_N2 / rho_Partikel_Mittel_2 * (
                                    g - index * 3 / 4 * zeta_2 / d_Partikel_Mittel_2 * abs(
                                u_sink_0_2 - Stroemungsgeschwindigkeit) ** 2)
                        u_sink_1_2 = u_sink_2_pre + a_2 * Delta_t
                        if abs(u_sink_1_2 - u_sink_0_2) / u_sink_0_2 * 100 <= 1:
                            break
                        u_sink_0_2 = u_sink_1_2
                    u_sink_2 = u_sink_1_2
                else:
                    while True:
                        index = np.sign(Stroemungsgeschwindigkeit - u_sink_0_2)
                        Re_2 = abs(u_sink_0_2 - Stroemungsgeschwindigkeit) * d_Partikel_Mittel_2 * rho_N2 / eta_N2
                        zeta_2 = 0.4
                        a_2 = g - rho_N2 / rho_Partikel_Mittel_2 * (g - index * 3 / 4 * zeta_2 / d_Partikel_Mittel_2 * abs(u_sink_0_2 - Stroemungsgeschwindigkeit) ** 2)
                        u_sink_1_2 = u_sink_2_pre + a_2 * (Verweilzeit_aus - Verweilzeit_ein)
                        u_sink_0_2 = u_sink_1_2
                    u_sink_2 = u_sink_1_2
        if Re_2 > 3.5:
            Nu_konv_K_2 = 2 + (0.4 * Re_2 ** 0.5 + 0.06 * Re_2 ** 0.067) * Pr_N2 ** 0.4 * (Pr_N2 / Pr_N2_T_W_P2) ** 0.25
            alpha_konv_K_2 = Nu_konv_K_2 * lambda_N2 / d_Partikel_Mittel_2
        else:
            alpha_konv_K_2 = 0

#---------------------------------------------------------------------------------------------------------------------------

        # free WÜ gas particle_1

        if w_1 > 0 and abs(Temperatur_Partikel_Eintritt_1 - Temperatur_Gas_Eintritt) > 0:
            Gr = g * beta_P1 * d_Partikel_Mittel_1 ** 3 * abs(Temperatur_Partikel_Eintritt_1 - Temperatur_Gas_Eintritt) / (eta_N2_P1 / rho_N2_P1) ** 2  # Grashof-Zahl
            Ra = Gr * Pr_N2_T_P1
            C1 = 1 / 16
            C2 = 0.52
            Nu_FK_K_1 = (1 / (C1 * Ra) ** (3 / 2) + 1 / (C2 * Ra ** (1 / 4)) ** (3 / 2)) ** (-2 / 3)
            alpha_FK_K_1 = Nu_FK_K_1 * lambda_N2_P1 / Durchmesser_innen
        else:
            alpha_FK_K_1 = 0

        # free WÜ gas particle_2

        if w_2 > 0 and abs(Temperatur_Partikel_Eintritt_2 - Temperatur_Gas_Eintritt) > 0:
            Gr = g * beta_P2 * d_Partikel_Mittel_2 ** 3 * abs(Temperatur_Partikel_Eintritt_2 - Temperatur_Gas_Eintritt) / (eta_N2_P2 / rho_N2_P2) ** 2  # Grashof-Zahl
            Ra = Gr * Pr_N2_T_P2
            C1 = 1 / 16
            C2 = 0.52
            Nu_FK_K_2 = (1 / (C1 * Ra) ** (3 / 2) + 1 / (C2 * Ra ** (1 / 4)) ** (3 / 2)) ** (-2 / 3)
            alpha_FK_K_2 = Nu_FK_K_2 * lambda_N2_P2 / Durchmesser_innen
        else:
            alpha_FK_K_2 = 0

        self.u_sink_1_result[i] = u_sink_1
        self.u_sink_2_result[i] = u_sink_2

#-------------------------------------------------------------------------------------------------------------------------

        U_Partikel_1 = alpha_konv_K_1 + alpha_str_GP_1 + alpha_FK_K_1
        U_Partikel_2 = alpha_konv_K_2 + alpha_str_GP_2 + alpha_FK_K_2

        # end Section V

        # Section VI: Calculation of temperature

        Tk = Temperatur_Gas_Eintritt  # in K
        TP1 = Temperatur_Partikel_Eintritt_1  # in K
        TP2 = Temperatur_Partikel_Eintritt_2  # in K

        Zellenlaenge = (Rohrlaenge_aus - Rohrlaenge_ein)  # Discretization
        Laengenschritt = min(Zellenlaenge / 100, 0.000001)  # [m]
        V_R_cell = Rohrquerschnittsflaeche * Laengenschritt
        cell_number = Zellenlaenge / Laengenschritt

        for i in range(int(Zellenlaenge / Laengenschritt)):
            Tk1 = Tk + ((Waermeerzeugung_gas * V_R_cell) + \
                        (U_innen * (Temperatur_Rohr_innen - Tk) * math.pi * Durchmesser_innen * Laengenschritt) + \
                        (U_Partikel_1 * (TP1 - Tk) * S_m_1 * m_Particle_1 / cell_number) + \
                        (U_Partikel_2 * (TP2 - Tk) * S_m_2 * m_Particle_2 / cell_number)) /  (Massenstrom_Gas * c_p_N2)  # K

            if w_1 == 0:
                TP11 = TP1
            else:
                TP11 = TP1 + ((Waermeerzeugung_particle_1 * V_R_cell) + (U_Partikel_1 * (Tk - TP1) + alpha_str_WP_1 * (Temperatur_Rohr_innen - TP1)) * S_m_1 * m_Particle_1 / cell_number) / (m_dot_Particle_1 * c_Partikel_1)  # K
            if w_2 == 0:
                TP21 = TP2
            else:
                TP21 = TP2 + ((Waermeerzeugung_particle_2 * V_R_cell) + (U_Partikel_2 * (Tk - TP2) + alpha_str_WP_2 * (Temperatur_Rohr_innen - TP2)) * S_m_2 * m_Particle_2 / cell_number) / (m_dot_Particle_2 * c_Partikel_2)  # K

            Tk = Tk1
            TP1 = TP11
            TP2 = TP21

        Temperatur_Gas_Austritt = Tk; # °C --> K
        Temperatur_Partikel_Austritt_1 = TP1;
        Temperatur_Partikel_Austritt_2 = TP2;

        return Temperatur_Gas_Austritt, Waermestromdichte, Temperatur_Partikel_Austritt_1, Temperatur_Partikel_Austritt_2


# end Section VI

# writeToCsv

# writes the object properties to a .csv file.

# INPUTS
# filepath: string . Complete file path with file name.
#
# OUTPUT
# None

    def writeToCsv(self, filepath):
        # Iterate over props
        props = dir(self)
        numProps = len(props)
        Names = []  # Empty list to store property names
        M = []  # Empty list to store property value

        for iProp in range(len(props)):

            if any(item in props[iProp] for item in ['xi']):
                Names.extend(["x_H2", "x_CO", "x_H2O", "x_CO2", "x_N2", "x_O2"])
                values = getattr(self, props[iProp])[:, :]
                for i, name in enumerate(["x_H2", "x_CO", "x_H2O", "x_CO2", "x_N2", "x_O2"]):
                    M.append(values[:, i])

            elif any(item in props[iProp] for item in ['obsReact1', 'obsReact2', 'effFuel1', 'effFuel2', 'int1', 'int2']):
                Names.append(f"{props[iProp]}_CO2")
                Names.append(f"{props[iProp]}_H2O")
                Names.append(f"{props[iProp]}_O2")
                values = getattr(self, props[iProp])[:, :]
                for i in range(3):
                    M.append(values[:, i])


            elif any(item in props[iProp] for item in ['tau', 'surfaceFuel1','dFuel1','dXdt','charConv'\
                        'overallConv','nGas','realTau','dVYdtFuel1','totalXC','Tgas','Tpart1',\
                        'annCO21','Tpart2','dVYdtfuel2','totalXfuel1','totalXfuel2','surfaceFuel2','dFuel2', \
                        'annCO22','charXFuel1','charXFuel2','dX1dt','dX2dt','volatileYield1','volatileYield2', \
                        'w1','w2','annFuncFuel1','annFuncFuel2','Fpc1','Fpc2','heatTransfer', \
                        'pyrolysisHeat', 'dYVdt','YV','YV1','YV2','CX1','CX2','surfaceDry1',\
                        'surfaceDry2','u_sink_1_result','u_sink_2_result','annFuel1','annFuel2 ','rhoFuel1','rhoFuel2']):
                Names.append(props[iProp])
                values = getattr(self, props[iProp])[:self.steps]
                M.append(values)
            #
            elif props[iProp]  == 'X' :
                Names.append(props[iProp])
                values = getattr(self, props[iProp])[:self.steps]
                M.append(values)

        # Write the data to a CSV file
        with open('filepathcsv.csv', mode='w', newline='') as file:
            writer = csv.writer(file)

            # Write the header row with Names
            writer.writerow(Names)

            # Write the data rows with values from M
            for i in range(len(M[0])):  # Assuming M contains rows of values
                row_values = [M[j][i] for j in range(len(M))]
                writer.writerow(row_values)

# customInitialization
# constructs an object of the class ReactorMesh and initializes its values given a certain amount of discretization steps. obj.Params need to be stored before calling it.
#
# INPUTS
# T: float. input temperature in K.
# steps: float. reactor discretization steps.

# OUTPUT
# None

    def customInitialization(self, T, steps):

        self.Tgas[0] = T + 273.15  # set gas temp to K
        self.Tpart1[0] = self.Params.Conditions['Tparticlein'] + 273.15 # Particle temperature in K
        self.Tpart2[0] = self.Params.Conditions['Tparticlein'] + 273.15

        self.w2[0] = self.Params.Conditions['w_2']  # Massenanteil Fuel 2
        self.w1[0] = 1 - self.w2[0]  # Massenanteil Fuel 1

        self.rhoFuel1 = self.Params.Fuel1['dryDensity']
        self.rhoFuel2 = self.Params.Fuel2['dryDensity']
        self.annFuel1 = self.Params.Fuel1['A_max_Annealing']
        self.annFuel2 = self.Params.Fuel2['A_max_Annealing']

        _, self.tau = self.initMeshResTime(T, steps)
        self.X = np.arange(0, self.Params.Geometry['pyrolysisLength'] + self.Params.Geometry['stepLengthPyrolysis'],\
                          self.Params.Geometry['stepLengthPyrolysis']).reshape(-1, 1)

#initMeshResTime
# updates the residence time of the gas in the reactor mesh.
#
# INPUTS
# T: float. input temperature in K.
# steps: float. reactor discretization steps.
#
# OUTPUTS
# x: array of floats. x coordinate of each incremental step in %.
# tau: array of floats. residence time in seconds at each step.

    def initMeshResTime(self, T, steps):

        # create vector of x
        x = np.arange(0, self.Params.Geometry['pyrolysisLength'] + self.Params.Geometry['stepLengthPyrolysis'],\
                          self.Params.Geometry['stepLengthPyrolysis']) / (steps - 1)

        # # input mols

        # # residence time (ideal gas assumption)
        totalTau = (self.Params.Geometry['tubeLength'] * self.Params.Geometry['innerRadius']**2 * np.pi * \
                    self.Params.Conditions['pressure'] * 10**5) / (np.sum(self.Params.n0) * 8.314 * T) # s
        # # vector with numPyrolysisSteps
        tau = np.linspace(0, totalTau, steps+1).reshape(-1,1)

        return  x,tau









