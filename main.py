# Note
# Derived from the MATLAB code main.m written by  authors Netter, Tobias (tobias.netter@tum.de) ,
# Ceruti, Amedeo (amedeo.ceruti@tum.de) within the project "Entrained Flow Gasification Reactor Modeling with MATLAB"

# This whole code needs to updated for a combustor
# Some parts of this code can be utilized
# Everthing is avalaile in other classes, just need to achieve a right flow based on theory

import os
import numpy as np
import pandas as pd
import warnings
from GlobalParameters import GlobalParameters
from GasProperties import GasProperties
from MixtureProperties import MixtureProperties
from ParticleDevelopment import ParticleDevelopment
from ReactorMesh import ReactorMesh
import computeSwellingRatio
import computeVM

jCase = 66 # Fuel input

Parameters = GlobalParameters("Parameters/db_combustion_reactor.xlsx","Parameters/reactor_geometry_combustion_reactor.xlsx",jCase ,"Parameters/reactor_conditions_combustion_reactor.xlsx")

FOLDERNAME = "Results"
os.mkdir(FOLDERNAME)

P = Parameters.Conditions['pressure']
Tin = Parameters.Conditions['Tin']
PARTICLEDISCRETIZATION = 99 # made smaller for the sake for less computation time
Parameters.Geometry['stepLengthPyrolysis'] = 1e-3
Parameters.Geometry['stepLengthCombustion'] = 1e-2
REACTORDISCRETIZATION = round(Parameters.Geometry['tubeLength'] / Parameters.Geometry['stepLengthPyrolysis'])
FILENAME = FOLDERNAME + "/" + Parameters.Fuel1['Name'] + "_" + "case_" + str(jCase) + ".csv"
flag = False  # flag to track if mesh was widened or not

Tmixing = Parameters.Conditions['Tmixing']

# Initialize RM
RM = ReactorMesh(Parameters, REACTORDISCRETIZATION, Tmixing)  # Construct reactor mesh object

# Init ParticleDevelopment
PD1 = ParticleDevelopment(PARTICLEDISCRETIZATION, REACTORDISCRETIZATION, Parameters.Fuel1, Parameters.Conditions)

# Init ParticleDevelopment for Fuel2
PD2 = ParticleDevelopment(PARTICLEDISCRETIZATION, REACTORDISCRETIZATION, Parameters.Fuel2, Parameters.Conditions)

# Compute dry fuel feed rate
coalFeedRateWF = Parameters.Conditions['coalFeedRate'] * ((1 - Parameters.Conditions['w_2']) * \
             (1 - Parameters.Fuel1['moisture'] / 100) + Parameters.Conditions['w_2'] * (1 - Parameters.Fuel2['moisture'] / 100))  # in kg/h
coalFeedRatewaf1 = Parameters.Conditions['coalFeedRate'] * (1 - Parameters.Conditions['w_2']) * \
                   (1 - Parameters.Fuel1['ashContent'] / 100 - Parameters.Fuel1['moisture'] / 100)  # in kg/h
coalFeedRatewaf2 = Parameters.Conditions['coalFeedRate'] * Parameters.Conditions['w_2'] * \
            (1 - Parameters.Fuel2['ashContent'] / 100 - Parameters.Fuel2['moisture'] / 100) # in kg/h
niGas = RM.computeElementalN(coalFeedRateWF)
RM.nGas[0], RM.xi[0, :], canteraOutput = RM.computeGasConcCantera(Tmixing+273.15, P, niGas) # update RM.xi, RM.nGas

RM.dFuel1[0] = PD1.computeMeanD()
RM.dFuel2[0] = PD2.computeMeanD()

# assign diffusion coefficients to ParticleDevelopment
PD1.Deff[:, :, 0], PD1.DM[:, :, 0] = PD1.computeDeff(RM.Tpart1[0], P, PD1.rpore[:, 0]*2, PD1.porosity[:, 0])
PD2.Deff[:, :, 0], PD2.DM[:, :, 0] = PD2.computeDeff(RM.Tpart2[0], P, PD2.rpore[:, 0]*2, PD2.porosity[:, 0])

RM.surfaceFuel1[0] = np.sum(PD1.xFrac[:, 0] / 100 * PD1.surface[:, 0])  # fuel surface
RM.surfaceFuel2[0] = np.sum(PD2.xFrac[:, 0] / 100 * PD2.surface[:, 0])

RM.surfaceDry1[0] = np.sum(PD1.xFrac[:, 0] / 100 * PD1.surfaceInDry[:, 0])  # fuel inner dry surface
RM.surfaceDry2[0] = np.sum(PD2.xFrac[:, 0] / 100 * PD2.surfaceInDry[:, 0])

# mass of ash in fuels
mAsh1 = Parameters.Conditions['coalFeedRate'] * (1 - Parameters.Conditions['w_2']) * Parameters.Fuel1['ashContent'] / 100  # in kg/h
mAsh2 = Parameters.Conditions['coalFeedRate'] * Parameters.Conditions['w_2'] * Parameters.Fuel2['ashContent'] / 100  # in kg/h

m_FC_C_1 = coalFeedRatewaf1 * (1 - Parameters.Conditions['YVtotal'])  # in kg/h
m_FC_C_2 = coalFeedRatewaf2 * (1 - Parameters.Conditions['YVtotal2'])  # in kg/h
m_FC_C = m_FC_C_1 + m_FC_C_2  # in kg/h

mVolatileMatter_1 = coalFeedRatewaf1 * Parameters.Conditions['YVtotal']; # in kg / h
mVolatileMatter_2 = coalFeedRatewaf2 * Parameters.Conditions['YVtotal2'] # in kg / h

# computeVM
# Computes elemental mole flow of volatile matter contained in a given fuel.
#
# INPUT
# waf:double.
# fuelDictionary: ordered dictionary.
# mVolatileMatter: double.
#
# OUTPUT
# [nH, nN, nS, nO, nC]: doubles. Molar flows in mol/s of each element.

def computeVM(waf, fuelDictionary, mVolatileMatter):
    H = waf * fuelDictionary["elementalH"] / 100  # in kg/h
    N = waf * fuelDictionary["elementalN"] / 100  # in kg/h
    S = waf * fuelDictionary["elementalS"] / 100  # in kg/h
    O = waf * fuelDictionary["elementalO"] / 100  # in kg/h
    mVM = mVolatileMatter - (H + N + S + O)  # in kg/h

    nH = H / 3600 * 1000 / 1   # mol/s
    nN = N / 3600 * 1000 / 14  # mol/s
    nS = S / 3600 * 1000 / 32  # mol/s
    nO = O / 3600 * 1000 / 16  # mol/s
    nC = mVM / 3600 * 1000 / 12  # mol/s

    return nH, nN, nS, nO, nC

n_VM_H_1, n_VM_N_1, n_VM_S_1, n_VM_O_1, n_VM_C_1 = computeVM(coalFeedRatewaf1, Parameters.Fuel1, mVolatileMatter_1)
n_VM_H_2, n_VM_N_2, n_VM_S_2, n_VM_O_2, n_VM_C_2 = computeVM(coalFeedRatewaf2, Parameters.Fuel2, mVolatileMatter_2)
m_VM_C_1 = n_VM_C_1 * 3600 / 1000 * 12
m_VM_C_2 = n_VM_C_2 * 3600 / 1000 * 12

# computeSwellingRatio
# Computes swelling ratio given previous step conditions.
#
# INPUT
# dYV: array of doubles. Example: ReactorMesh.dYV1
# VolatileYieldMax: double. Example: 0.663
# SW_max: double. Example: Parameters.Fuel2.SW
# i: integer. Current Reactor Step.
#
# OUTPUT
# SW_ratio_1: double. Swelling ratio.

def computeSwellingRatio(dYV, SW_1_max, VolatileYieldMax, i):
    if i == 2:
        SW_ratio_1 = 1
    else:
        SW_1_pre = 1 + (SW_1_max - 1) * dYV[i - 2] / VolatileYieldMax
        SW_1 = 1 + (SW_1_max - 1) * dYV[i - 1] / VolatileYieldMax
        SW_ratio_1 = SW_1 / SW_1_pre

    return SW_ratio_1

#Note
# This values are assumed as array of size depenedent on  REACTORDISCRETIZATION
# Need to look more into it
RM.annFuel1 = np.array(REACTORDISCRETIZATION * [RM.annFuel1])
RM.annFuel2 = np.array(REACTORDISCRETIZATION * [RM.annFuel2])
RM.rhoFuel1 = np.array(REACTORDISCRETIZATION * [RM.rhoFuel1])
RM.rhoFuel2 = np.array(REACTORDISCRETIZATION * [RM.rhoFuel2])

for i in range(1, REACTORDISCRETIZATION):
    # Real Residence Time
    RM.realTau[i] = (RM.X[i] - RM.X[i - 1]) * Parameters.Geometry['tubeDiameter']**2 * 0.25 * \
                   np.pi * P * 1E5 / (RM.nGas[i - 1] * 8.314 * RM.Tgas[i - 1]) + RM.realTau[i - 1]

    # Anirudh Just for running the code making a assumption
    RM.realTau = np.array([0, 0.00116285, 0.00116285, 0.003488099, 0.004650499, 0.00581275])

    # Annealing
    RM.annFuel1[i], RM.annFuncFuel1[i] = RM.computeAnnealingStep(RM.Tpart1[i - 1], Parameters.Fuel1, i)
    RM.annFuel2[i], RM.annFuncFuel2[i] = RM.computeAnnealingStep(RM.Tpart2[i - 1], Parameters.Fuel2, i)

    # Surface f_pc by pyrolysis
    RM.Fpc1[i] = RM.computeFuelSurface(RM.Tpart1[i - 1], RM.Fpc1[i - 1], Parameters.Fuel1, i)
    RM.Fpc2[i] = RM.computeFuelSurface(RM.Tpart2[i - 1], RM.Fpc2[i - 1], Parameters.Fuel2, i)

    # Pyrolysis
    if RM.YV[i - 1] < Parameters.Conditions['YV'] * 0.999:
        # Pyrolysis-Parameter
        (RM.volatileYield1[i], RM.volatileYield2[i], RM.YV1[i], RM.YV2[i], RM.dVYdtFuel1[i], RM.dVYdtfuel2[i],
         RM.pyrolysisHeat[i], RM.YV[i]) = RM.computePyrolysis(coalFeedRatewaf1, coalFeedRatewaf2, i)

        # dYV/dt_Total
        RM.dYVdt[i] = (RM.YV[i] - RM.YV[i - 1]) / (RM.realTau[i] - RM.realTau[i - 1])

        # Swelling through pyrolysis
        SW_ratio_1 = computeSwellingRatio(RM.YV1, Parameters.Fuel1['SW'], Parameters.Conditions['YVtotal'], i)
        PD1.D[:, i] = SW_ratio_1 * PD1.D[:, i]

        SW_ratio_2 = computeSwellingRatio(RM.YV2, Parameters.Fuel2['SW'], Parameters.Conditions['YVtotal2'], i)
        PD2.D[:, i] = SW_ratio_2 * PD2.D[:, i]
    else:
        RM.YV1[i] = RM.YV1[i - 1]  # Pyrolysis conversion 1. coal
        RM.YV2[i] = RM.YV2[i - 1]  # Pyrolysis conversion 2. coal
        RM.YV[i] = RM.YV[i - 1]  # total pyrolysis turnover
        SW_ratio_1 = 1  # reset SW_ratios
        SW_ratio_2 = 1


# Note
# So now it is gasification in the old main matlab file. There are certain conditions that need to be fullfill for
# gasification to start,that means that the pyrolysis should end , so at what point does gasification starts depends on our research
#
# Same logic has to be applied to matlab like at what point pyrolysis end and combustion start, and on the basis of
#that below written code should me modified

    # Combustion
    # Combustion -Start-Condition (this needs to be found)
    if RM.YV1[i - 1] < Parameters.Conditions['YVtotal'] * 0.:
        combustionFuel1 = False
    else:
        combustionFuel1 = True # Gasification of coal 1 starts

    if RM.YV2[i - 1] < Parameters.Conditions['YVtotal2'] * 0.:
        combustionFuel2 = False
    else:
        combustionFuel2 = True  # # Gasification of coal 2 starts

    RM.int1[i, :] = RM.computeIntrinsicReactivity(RM.Tpart1[i - 1], P * RM.xi[i - 1, [3, 2, 5]], Parameters.Fuel1)
    RM.int2[i, :] = RM.computeIntrinsicReactivity(RM.Tpart2[i - 1], P * RM.xi[i - 1, [3, 2, 5]], Parameters.Fuel2)

    PD1.Deff[:, :, i], PD1.DM[:, :, i] = PD1.computeDeff(RM.Tpart1[i - 1], P, PD1.rpore[:, i - 1] * 2,PD1.porosity[:, i - 1])
    PD2.Deff[:, :, i], PD2.DM[:, :, i] = PD2.computeDeff(RM.Tpart2[i - 1], P, PD2.rpore[:, i - 1] * 2,PD2.porosity[:, i - 1])

    for jStep in range(0, PARTICLEDISCRETIZATION ):

        # Effectiveness of [CO2, H2O, O2] for PD1
        PD1.eff[:, jStep, i] = (PD1.computeEffectiveness(RM.Tpart1[i - 1], P * RM.xi[i - 1, [3, 2, 5]], RM.annFuel1[i],RM.int1[i, :], i, jStep)).flatten()
        # Effectiveness of [CO2, H2O, O2] for PD2
        PD2.eff[:, jStep, i] = (PD2.computeEffectiveness(RM.Tpart2[i - 1], P * RM.xi[i - 1, [3, 2, 5]], RM.annFuel2[i],RM.int2[i, :], i, jStep)).flatten()

        # If combustionFuel1 for PD1
        (PD1.robs[:, jStep, i], PD1.robspore[:, jStep, i], PD1.robsfilm[:, jStep, i], PD1.robsin[:, jStep, i],PD1.regime[:, jStep, i]) = PD1.computeMasstransportLimit(RM.Tgas[i - 1], P, RM.xi[i - 1, [3, 2, 5]],RM.annFuel1[i], RM.int1[i, :], i, jStep)
        # If gasificationFuel1 for PD2
        (PD2.robs[:, jStep, i], PD2.robspore[:, jStep, i], PD2.robsfilm[:, jStep, i], PD2.robsin[:, jStep, i],PD2.regime[:, jStep, i]) = PD2.computeMasstransportLimit(RM.Tgas[i - 1], P, RM.xi[i - 1, [3, 2, 5]],RM.annFuel2[i], RM.int2[i, :], i, jStep)

        # Update char conversion for PD1
        PD1.dXdt[jStep, i], PD1.charConv[jStep, i] = PD1.updateCharConv(RM.YV1[i - 1], Parameters.Conditions['YVtotal'],RM.realTau, jStep, i)
        # Update char conversion for PD2
        PD2.dXdt[jStep, i], PD2.charConv[jStep, i] = PD2.updateCharConv(RM.YV2[i - 1], Parameters.Conditions['YVtotal2'],RM.realTau, jStep, i)

        # Update particle size and density for PD1
        (PD1.D[jStep, i], PD1.rho[jStep, i], PD1.surface[jStep, i], PD1.surfaceIn[jStep, i], PD1.rpore[jStep, i],PD1.porosity[jStep, i], PD1.surfaceDry[jStep, i], PD1.surfaceInDry[jStep, i],PD1.vol2C[jStep, i]) = PD1.updateParticleSize(RM.Fpc1, RM.YV1[i], Parameters.Conditions['rf'],Parameters.Conditions['YVtotal'], Parameters.Conditions['YVtotal'],jStep, i)
        # Update particle size and density for PD2
        (PD2.D[jStep, i], PD2.rho[jStep, i], PD2.surface[jStep, i], PD2.surfaceIn[jStep, i], PD2.rpore[jStep, i],PD2.porosity[jStep, i], PD2.surfaceDry[jStep, i], PD2.surfaceInDry[jStep, i],PD2.vol2C[jStep, i]) = PD2.updateParticleSize(RM.Fpc2, RM.YV2[i], Parameters.Conditions['rf_2'],Parameters.Conditions['YVtotal2'], Parameters.Conditions['YVtotal2'],jStep, i)

    PD1.vFrac[:, i] = ((PD1.D[:, i] ** 3) * PD1.qtyFrac / np.sum((PD1.D[:, i] ** 3) * PD1.qtyFrac)) * 100
    dichteMischung = np.sum(PD1.vFrac[:, i]* PD1.rho[:, i]/ 100)
    PD1.xFrac[:, i] = PD1.rho[:, i]/ dichteMischung * PD1.vFrac[:, i]

    # mass-average values
    RM.effFuel1[i, :] = np.dot(PD1.xFrac[:, i].T / 100, PD1.eff[:, :, i].T)
    RM.dX1dt[i] = np.sum(PD1.xFrac[:, i] / 100 * PD1.dXdt[:, i])
    RM.charXFuel1[i] = np.sum(PD1.xFrac[:, i] / 100 * PD1.charConv[:, i])
    RM.dFuel1[i] = np.sum(PD1.xFrac[:, i] / 100 * PD1.D[:, i])
    RM.rhoFuel1[i] = np.sum(PD1.xFrac[:, i] / 100 * PD1.rho[:, i])
    RM.surfaceFuel1[i] = np.sum(PD1.xFrac[:, i] / 100 * PD1.surface[:, i])
    RM.surfaceDry1[i] = np.sum(PD1.xFrac[:, i] / 100 * PD1.surfaceInDry[:, i])
    RM.obsReact1[i, :] = np.dot(PD1.xFrac[:, i].T / 100, PD1.robs[:, :, i].T)

    PD2.vFrac[:, i] = ((PD2.D[:, i] ** 3) * PD2.qtyFrac / np.sum((PD2.D[:, i] ** 3) * PD2.qtyFrac)) * 100
    dichteMischung = np.sum(PD2.vFrac[:, i]* PD2.rho[:, i]/ 100)
    PD2.xFrac[:, i] = PD2.rho[:, i]/ dichteMischung * PD2.vFrac[:, i]

    # mass-average values
    RM.effFuel2[i, :] = np.dot(PD2.xFrac[:, i].T / 100, PD2.eff[:, :, i].T)
    RM.dX2dt[i] = np.sum(PD2.xFrac[:, i] / 100 * PD2.dXdt[:, i])
    RM.charXFuel2[i] = np.sum(PD2.xFrac[:, i] / 100 * PD2.charConv[:, i])
    RM.dFuel2[i] = np.sum(PD2.xFrac[:, i] / 100 * PD2.D[:, i])
    RM.rhoFuel2[i] = np.sum(PD2.xFrac[:, i] / 100 * PD2.rho[:, i])
    RM.surfaceFuel2[i] = np.sum(PD2.xFrac[:, i] / 100 * PD2.surface[:, i])
    RM.surfaceDry2[i] = np.sum(PD2.xFrac[:, i] / 100 * PD2.surfaceInDry[:, i])
    RM.obsReact2[i, :] = np.dot(PD2.xFrac[:, i].T / 100, PD2.robs[:, :, i].T)

    RM.YV[i] = max(RM.YV[i], RM.YV[i - 1])
    RM.YV1[i] = max(RM.YV1[i], RM.YV1[i - 1])
    RM.YV2[i] = max(RM.YV2[i], RM.YV2[i - 1])
    # Char_Conversion

    RM.charConv[i] = (m_FC_C_1 / m_FC_C) * RM.charXFuel1[i] + (m_FC_C_2 / m_FC_C) * RM.charXFuel2[i]

    # Note
    # Pyrolysis revenue increase due to gasification.
    # These all calculatuions are for gasification, needs to be updated for combustion

    Delta_YV_1 = (Parameters.Conditions['YVtotal'] - RM.YV1[i]) * (RM.charXFuel1[i] - RM.charXFuel1[i - 1])
    Delta_YV_2 = (Parameters.Conditions['YVtotal2'] - RM.YV2[i]) * (RM.charXFuel2[i] - RM.charXFuel2[i - 1])
    Delta_YV = (Delta_YV_1 * coalFeedRatewaf1 + Delta_YV_2 * coalFeedRatewaf2) / (coalFeedRatewaf1 + coalFeedRatewaf2) # Add additional pyrolysis heat
    RM.pyrolysisHeat[i] = RM.pyrolysisHeat[i] / RM.YV[i] * (RM.YV[i] + Delta_YV)

    # Pyrolysis turnover with additional release of volatiles due to gasification.
    RM.YV1[i] = RM.YV1[i] + Delta_YV_1
    RM.YV2[i] = RM.YV2[i] + Delta_YV_2
    RM.YV[i] = (RM.YV1[i] * coalFeedRatewaf1 + RM.YV2[i] * coalFeedRatewaf2) / (coalFeedRatewaf1 + coalFeedRatewaf2) # Total pyrolysis conversion

    # Overall_Conversion
    RM.totalXfuel1[i] = RM.YV1[i] + RM.charXFuel1[i] * (1 - Parameters.Conditions['YVtotal'])
    RM.totalXfuel2[i] = RM.YV2[i] + RM.charXFuel2[i] * (1 - Parameters.Conditions['YVtotal'])
    RM.overallConv[i] = (1 - Parameters.Conditions['w_2']) * RM.totalXfuel1[i] + Parameters.Conditions['w_2'] * RM.totalXfuel2[i]

    # Carbon-Conversion
    RM.CX1[i] = max((RM.YV1[i] / Parameters.Conditions['YVtotal'] * m_VM_C_1 + RM.charXFuel1[i] * m_FC_C_1) / (coalFeedRatewaf1 * Parameters.Fuel1['elementalC'] / 100), 0) # (fugitive carbon+gasifiedFixedCarbon)/ total carbon
    RM.CX2[i] = max((RM.YV2[i] / Parameters.Conditions['YVtotal'] * m_VM_C_2 + RM.charXFuel2[i] * m_FC_C_2) / (coalFeedRatewaf2 * Parameters.Fuel2['elementalC'] / 100), 0) # (fugitive carbon+gasifiedFixedCarbon)/ total carbon
    RM.totalXC[i] = (RM.CX1[i] * coalFeedRatewaf1 * Parameters.Fuel1['elementalC'] + RM.CX2[i] * coalFeedRatewaf2 * Parameters.Fuel2['elementalC']) / (coalFeedRatewaf1 * Parameters.Fuel1['elementalC'] + coalFeedRatewaf2 * Parameters.Fuel2['elementalC'])

    # GasParticle temperature

    # RM.Tgas[i], RM.heatTransfer[i], RM.Tpart1[i], RM.Tpart2[i] = RM.computeGasPartT(P, coalFeedRatewaf1, coalFeedRatewaf2,m_FC_C_1, m_FC_C_2, mAsh1, mAsh2, i,canteraOutput, PD1, PD2)
    # Assuming values for above as too many values are missing for this function

    RM.Tgas = np.array([400,500,650,700,800,950])
    RM.heatTransfer = np.array([1000,2000,3000,4000,500,6000])
    RM.Tpart1 = np.array([350,487.23,550.23,600,750,850,659])
    RM.Tpart2 = np.array([300,450,550,650,750,900])

    # CANTERA
    n_C_NextStep = niGas[0] + RM.totalXC[i] * (coalFeedRatewaf1 * Parameters.Fuel1['elementalC'] / 100 + coalFeedRatewaf2 * Parameters.Fuel2['elementalC'] / 100) / 3600 / (12 / 10 ** 3)
    n_H_NextStep = niGas[1] + RM.YV1[i] / Parameters.Conditions['YVtotal'] * n_VM_H_1 + RM.YV2[i] / Parameters.Conditions['YVtotal2'] * n_VM_H_2
    n_O_NextStep = niGas[2] + RM.YV1[i] / Parameters.Conditions['YVtotal'] * n_VM_O_1 + RM.YV2[i] / Parameters.Conditions['YVtotal2'] * n_VM_O_2
    n_N_NextStep = niGas[3] + RM.YV1[i] / Parameters.Conditions['YVtotal'] * n_VM_N_1 + RM.YV2[i] / Parameters.Conditions['YVtotal2'] * n_VM_N_2


# Note
# RM.porosity is a matrix and not an array but here it is used as an array.Its wrong

    # if RM.YV1[i] < Parameters.Conditions['YVtotal'] * 0.995:
    #     PD1.surfaceIn[i] = PD1.surfaceIn[i] / RM.Fpc1[i - 1]
    # else:
    #     PD1.porosity[i] = PD1.porosity[i - 1]
    #     PD1.surfaceIn[i] = PD1.surfaceIn[i - 1]
    #
    # if RM.YV2[i] < Parameters.Conditions['YVtotal2'] * 0.995:
    #     PD2.surfaceIn[i] = PD2.surfaceIn[i] / RM.Fpc2[i - 1]
    # else:
    #     PD2.porosity[i] = PD2.porosity[i - 1]
    #     PD2.surfaceIn[i] = PD2.surfaceIn[i - 1]

    # All this code beow is specifically for pyrolysis
    # Just copied it from the main matlab code , it definitely needs to be modified for combsution
    # Adjustment of length when turnover decreases (gasification
    # tract is ten times faster)

    # if not flag and (((RM.overallConv[i] - RM.overallConv[i - 1]) / RM.overallConv[i - 1]) / (RM.X[i] - RM.X[i - 1])) < 10:
    #     PyrolysisSteps = i
    #     PyrolysisLength = RM.X[i]
    #     GasificationSteps = int((Parameters.Geometry['tubeLength'] - PyrolysisLength) / Parameters.Geometry['stepLengthGasification'])
    #     resTime = Parameters.computeResidenceTime(Tmixing + 273.15, P, RM.nGas[0])
    #     RM.X[PyrolysisSteps:(PyrolysisSteps + GasificationSteps)] = np.arange(PyrolysisLength,Parameters.Geometry['tubeLength'],Parameters.Geometry['stepLengthGasification'])
    #
    #     flag = true
    #
    #     if RM.X[PyrolysisSteps + GasificationSteps] < Parameters.Geometry['tubeLength']:
    #         REACTORDISCRETIZATION = PyrolysisSteps + GasificationSteps + 1  # force to end at length
    #         RM.X[REACTORDISCRETIZATION] = Parameters.Geometry['tubeLength']
    #     else:
    #         REACTORDISCRETIZATION = PyrolysisSteps + GasificationSteps  # force to end at length
    #
    #     for index in range(PyrolysisSteps, PyrolysisSteps + GasificationSteps):
    #         # Residence time (Gas input)
    #         RM.tau[index] = (RM.X[index] - RM.X[index - 1]) / Parameters.Geometry['tubeLength'] * resTime + RM.tau[index - 1]  # Residence Time (Aufteilung der Verweilzeit in gleichgroße Segmente)
    #
    #     if i > REACTORDISCRETIZATION or RM.X[i] >= Parameters.Geometry['tubeLength']:
    #         break
    #
    #     if i > 50 and abs((RM.overallConv[i - 10] - RM.overallConv[i - 1]) / RM.overallConv[i - 10]) / max((RM.X[i] - RM.X[i - 10]), 0.1) < 1e-4:  # converged, ends up breaking the code
    #         print('Converged, ending simulation to avoid crash.')
    #         break
    #
    #     if i % 50 == 0:
    #         print(f'iteration {i}')

# Section IX: Mass fraction via pyrolysis
# Evolution of the mass fraction over the course of pyrolysis+gasification

w_1_neu = (1 - RM.totalXfuel1) * (1 - Parameters.Conditions['w_2'])
w_2_neu = (1 - RM.totalXfuel2) * Parameters.Conditions['w_2']

RM.w1 = w_1_neu / (w_1_neu + w_2_neu)
RM.w2 = w_2_neu / (w_1_neu + w_2_neu)

RM.writeToExcel("filepathcsv.csv")