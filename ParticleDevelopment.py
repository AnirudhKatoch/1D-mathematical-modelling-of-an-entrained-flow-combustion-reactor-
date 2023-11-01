# ParticleDevelopment

# Note
# Derived from the MATLAB class ParticleDevelopment written by  authors Netter, Tobias (tobias.netter@tum.de) ,
# Ceruti, Amedeo (amedeo.ceruti@tum.de) within the project "Entrained Flow Gasification Reactor Modeling with MATLAB"

# Contains properties as arrays of doubles with size discretization x steps. These can be accessed and modified by the
# main script in order to procede with the simulation of each particle with the information contained in this class.

# Example
# PD1 = ParticleDevelopment(99, 2000, Parameters.Fuel1, Parameters.Conditions)

# Contents
#
# Properties
# ParticleDevelopment
# computeMeanD
# computeDeff
# computeEffectiveness
# computeMasstransportLimit
# updateCharConv
# updateParticleSize
# initParticleDev
# initParticleStructure

import numpy as np
import math

class ParticleDevelopment:

    fuelParams = None # GlobalParameters.Fuel dictionary
    discretization = None # integer: discretization num of the particle
    steps = None # integer: discretization num of the reactor

    # stored props in matrices of discretization x steps
    D = None  # 'Diameter [10^-6 m]'
    vFrac = None  # 'Volume fraction [%]'
    xFrac = None  # 'Mass fraction [%]'
    surface = None  # 'Surface [m²/g]'
    rho = None  # 'Density [kg/m³]'
    eff = None  # [CO2, H2O, O2] -> 'EffCO2 [-]'
    dXdt = None  # 1/s
    charConv = None  # [-]
    robspore = None  # [CO2, H2O, O2] -> 'r_obs_pore_O2 [g/g*s]'
    robsfilm = None  # idem -> 'r_obs_film_CO2 [g/g*s]'
    robs = None  # idem -> 'r_obs_H2O [g/g*s]'
    regime = None  # [CO2, H2O, O2]; cell array
    surfaceIn = None
    robsin = None  # 'r_obs_intern_CO2 [g/g*s]'
    Deff = None  # [CO2, H2O, O2] -> 'D_eff_CO2'
    DM = None  # [CO2, H2O, O2] -> 'D_M_CO2'
    rpore = None  # 'r_pore'
    porosity = None  # 'Porosity_1'
    surfaceDry = None  # 'Surface_dry'
    surfaceInDry = None  # 'Surface_in_dry'
    vol2C = None  # 'Volatile/Carbon'

    # from particlestructurevar_init
    qtyFrac = None  # Partikelgroessenverteilung_diskret_anzahl_2{1,2}='Quantity fraction [%]'
    dAshParticle = None  # {1,6}= d_AshParticle_init_1 [10^6m]'
    mFixedCarbon = None  # ParticleStructure_init{1,7}='m_FixedCarbon_init_1 [kg]';
    mVolatile = None  # m_Volatile_init [kg]
    # Po = None  # Po_init [-]

# Construct an object of the ParticleDevelopment class.
#
# INPUTS
#
# discretization: integer  discretization num of the particle
# steps: integer discretization num of the reactor
# fueldictionary: # GlobalParameters.Fuel dictionary
#  Conditions: # GlobalParameters.Conditions

    def __init__(self,discretization, steps, fueldictionary, Conditions):

        self.fuelParams = fueldictionary  # fuel parameters
        self.discretization = discretization  # particle discretization
        self.steps = steps  # num of reactor length steps
        # initialize the values with the initialization function.
        self.initParticleDev(fueldictionary, Conditions)  # initialize step 0

# computeMeanD
# computes the mean particle diameter given the stored mass fractions in self.xFrac and the cumulative diameter in self.D
#
# INPUT
# None
#
# OUTPUT
# D: float. Mean particle diameter.

    def computeMeanD(self):
        D = np.sum(self.xFrac[:, 0] * self.D[:, 0] / 100)
        return D

# computeDeff
# Compute effective diffusion coefficients
#
# INPUT
# temperature: float. Temp in K.
# pressure: float. Pressure in bar
# dPore: array of floats. Pore Diameter in m.
# epsilon: array of floats

    def computeDeff(self, temperature, pressure, dPore, epsilon):

        dPore = np.array(dPore)
        epsilon = np.array(epsilon)

        # OPERATING PARAMETERS
        # Print
        # [bar]
        p = pressure / 1.01325  # [atm]! --> Mass transfer script (Tremel: in bar)
        # Temperature
        T = temperature  # [K]
        # SUBSTANCE PARAMETERS & CONSTANTS
        # molar masses
        M_CO2 = 44  # [kg/kmol] = g/mol
        M_H2O = 18  # [kg/kmol]
        M_N2 = 28  # [kg/kmol]
        M_O2 = 32  # [kg/kmol]

        # Boltzmann constant
        k_b = 1.3806504e-23  # [J/K]
        # universal gas constant
        R = 8.3145  # [J/(mol K)]

        # LOOK-UP-TABLES: Collision Integral & Lennard-Jones Parameters for Molecular Diffusion
        s = np.loadtxt('parameters/Look_up_Diffusionkoefficient')

        # MOLECULAR DIFFUSION
        # CO2
        sigma_CO2 = (s.Lennard_Jones_Parameter[1, 1] + s.Lennard_Jones_Parameter[3, 1]) / 2  # Lennard-Jones-Parameter [Angström]
        epsilon_CO2 = (s.Lennard_Jones_Parameter[1, 2] * s.Lennard_Jones_Parameter[3, 2] * k_b ** 2) ** 0.5  # Lennard-Jones-Parameter [J]
        # H2O
        sigma_H2O = (s.Lennard_Jones_Parameter[2, 1] + s.Lennard_Jones_Parameter[3, 1]) / 2  # Lennard-Jones-Parameter [Angström]
        epsilon_H2O = (s.Lennard_Jones_Parameter[2, 2] * s.Lennard_Jones_Parameter[3, 2] * k_b ** 2) ** 0.5  # Lennard-Jones-Parameter [J]
        # O2
        sigma_O2 = (s.Lennard_Jones_Parameter[4, 1] + s.Lennard_Jones_Parameter[3, 1]) / 2  # Lennard-Jones-Parameter [Angström]
        epsilon_O2 = (s.Lennard_Jones_Parameter[4, 2] * s.Lennard_Jones_Parameter[5, 2] * k_b ** 2) ** 0.5  # Lennard-Jones-Parameter [J]

        for i in range(2, 37):
            # CO2
            if s.Kollisionsintegral[i, 0] <= k_b * T / epsilon_CO2 and k_b * T / epsilon_CO2 <= s.Kollisionsintegral[i+1, 0]:
                Omega_CO2 = (k_b * T / epsilon_CO2 - s.Kollisionsintegral[i, 0]) / (s.Kollisionsintegral[i + 1, 0] - \
                            s.Kollisionsintegral[i, 0]) * (s.Kollisionsintegral[i + 1, 1] - s.Kollisionsintegral[i, 1]) + \
                            s.Kollisionsintegral[i, 1]

            # H2O
            if s.Kollisionsintegral[i, 0] <= k_b * T / epsilon_H2O and k_b * T / epsilon_H2O <= s.Kollisionsintegral[i+1, 0]:
                Omega_H2O = (k_b * T / epsilon_H2O - s.Kollisionsintegral[i, 0]) / (s.Kollisionsintegral[i + 1, 0] - \
                            s.Kollisionsintegral[i, 0]) * (s.Kollisionsintegral[i + 1, 1] - s.Kollisionsintegral[i, 1]) + \
                            s.Kollisionsintegral[i, 1]

            # O2
            if s.Kollisionsintegral[i, 0] <= k_b * T / epsilon_O2 and k_b * T / epsilon_O2 <= s.Kollisionsintegral[i+1, 0]:
                Omega_O2 = (k_b * T / epsilon_O2 - s.Kollisionsintegral[i, 0]) / (s.Kollisionsintegral[i + 1, 0] - \
                           s.Kollisionsintegral[i, 0]) * (s.Kollisionsintegral[i + 1, 1] - s.Kollisionsintegral[i, 1]) + \
                           s.Kollisionsintegral[i, 1]

        # deviation between tremel fct. and look - up - table:
        MW = np.array([M_CO2,M_H2O,M_O2])
        sigma = np.array([sigma_CO2,sigma_H2O,sigma_O2])
        omega = np.array([Omega_CO2,Omega_H2O,Omega_O2])

        DM = 1.8583e-7 * T ** 1.5 * (1 / MW + 1 / M_N2) ** 0.5 / (p * sigma ** 2 * omega)

        # KNUDSEN-DIFFUSION
        DK = dPore / 3 * (8 * R * T * 10 ** 3 / (math.pi * MW)) ** 0.5
        tau = 1 / epsilon
        Deff = (epsilon / tau) * (1 / (1 / DM + 1 / DK))
        DM = DM * np.ones(Deff.shape)
        # O2 [m²/s]
        DM = DM.transpose()
        DM = DM.tolist()
        Deff = Deff.transpose()
        Deff = Deff.tolist()

        return Deff, DM

# computeEffectiveness
# Computes effectiveness-factor of both fuel types at iteration i and particle discretization step j.
#
# INPUTS
#
# T:double. Temperature in K.
# Pi: 3x1 array of floats. vector of partial pressures in bar [CO2, H2O, O2].
# ann: array of floats. Annealing values. Ex: ReactorMesh.ann1.
# int: array of floats. Intrinsic reactivity. Ex: ReactorMesh.ann1.
# i: integer. discretization step (index).
# j: integer. Particle discretization step
#
# OUTPUT
# Eta: array of doubles. Order: [CO2, H2O, O2].

    def computeEffectiveness(self, T, Pi, ann, int, i, j):

        M_C = 12 / (10 ** 3)  # kg/mol
        nu = np.array([1, 1, 0.5])

        def eta(T, ann, int):

            m_Carbon = (1 - self.charConv[j, i - 1]) * self.mFixedCarbon[j]
            V_p = (math.pi / 6) * (self.D[j, i - 1] / 1E6) ** 3

            # Thiele modulus

            partial_pressures = np.array(self.fuelParams['n'])
            Phi = self.D[j, i - 1] / (10 ** 6) / 6 * ((1 + partial_pressures) *(ann * int) * self.surfaceIn[j, i - 1] * m_Carbon * 8.3145 * T * nu /(2 * M_C *self.Deff[:, j, i - 1].transpose() * V_p * Pi * 1E5)) # multiplied by 1E5 because partial pressures is in bar

            # Effectiveness factor

            if self.surfaceIn[j, i - 1] == 0:
                Eta = np.zeros((3, 1))
            else:
                Eta = (1 / Phi) * (1 / np.tanh(3 * Phi) - 1 / (3 * Phi))  # Tremel

            # Loop to check the values
            for iComp in range(Eta.shape[0]):
                Eta[iComp] = max(min(Eta[iComp], 1), 0)
            return Eta
        Eta = eta(T, ann, int)
        return Eta

# computeMasstransportLimit
# Computes effectiveness-factor of both fuel types at iteration i and particle discretization step j.

# INPUT
# T: double. Temperature in K.
# P: double. pressure in Pa
# xi: 3x1 array of floats. molar fracs of components CO2, H2O, O2 in this order at the previopus reactor step. Example: ReactorMesh(i-1, [4,3,6]).
# ann: float. Annealing values at reactor step i and particle discretization j. Ex: ReactorMesh.ann1.
# int: 3x1 array of floats. Intrinsic reactivity. Ex: ReactorMesh.ann1.
# i: integer. discretization step (index).
# j: integer. Particle discretization step.

# OUTPUT
# robs: array of floats. Order: [CO2, H2O, O2]
# robspore: array of floats. Order: [CO2, H2O, O2]
# robsfilm: array of floats. Order: [CO2, H2O, O2]
# robsin: array of floats. Order: [CO2, H2O, O2]
# Regime: cell array of strings. Rate determining mass transport regime in order: [CO2, H2O, O2]. Possible strings: 'Film', 'Pore'.

    def computeMasstransportLimit(self, T, P, xi, ann, int, i, j):

        # Old mass transport limitation function.
        # Calculation of the mass transport limitation: differentiation between
        # regime I&II (pore diffusion) and regime III(film diffusion)

        # Assign variable names
        Surface = self.surface[j, i - 1]
        Annealing = ann
        eta_CO2 = self.eff[0, j, i]
        eta_H2O = self.eff[1, j, i]
        eta_O2 = self.eff[2, j, i]
        r_int_CO2 = int[0]
        r_int_H2O = int[1]
        r_int_O2 = int[2]
        diameter = self.D[j, i - 1]
        density = self.rho[j, i - 1]
        D_M_CO2 = self.DM[0, j, i - 1]
        D_M_H2O = self.DM[1, j, i - 1]
        D_M_O2 = self.DM[2, j, i - 1]
        y_CO2 = xi[0]
        y_H2O = xi[1]
        y_O2 = xi[2]
        Surface_in = self.surfaceIn[j, i - 1]

        # Universelle Gaskonstante
        R = 8.314  # [J/(mol K)]

        # Gravimetric-stoichiometric coefficient
        N_CO2 = 6 / 11  # [-] Waldemar 3/11
        # N_CO2 = 3 / 11  # [-] Wa
        N_H2O = 2 / 3  # [-]
        N_O2 = 3 / 8  # [-]

        # Molare Masses
        M_CO2 = 0.044  # [kg/mol]
        M_H2O = 0.018  # [kg/mol]
        M_O2 = 0.032  # [kg/mol]

        # Unit conversion
        # pressure = P * 10 ** 5  # [bar]-->[Pa]
        pressure = P
        diameter = diameter * 10 ** -6  # [um]-->[m]

        # Regime I&II (pore diffusion)
        # r_obs_pore_CO2 = Surface * Annealing_CO2 * eta_CO2 * r_int_CO2 # [g/g*s] resp. [1/s] %%JH: formula commented out because of not considering different anealing
        r_obs_pore_CO2 = Surface * Annealing * eta_CO2 * r_int_CO2  # [g/g*s] bzw. [1/s]
        r_obs_pore_intern_CO2 = Surface_in * Annealing * eta_CO2 * r_int_CO2  # [g/g*s] bzw. [1/s]
        r_obs_pore_H2O = Surface * Annealing * eta_H2O * r_int_H2O  # [g/g*s] bzw. [1/s]
        r_obs_pore_intern_H2O = Surface_in * Annealing * eta_H2O * r_int_H2O  # [g/g*s] bzw. [1/s]
        r_obs_pore_O2 = Surface * Annealing * eta_O2 * r_int_O2  # [g/g*s] bzw. [1/s]
        r_obs_pore_intern_O2 = Surface_in * Annealing * eta_O2 * r_int_O2  # [g/g*s] bzw. [1/s]

        # Regime III (Filmdiffusion)
        r_obs_film_CO2 = 12 / (diameter ** 2 * density) * pressure / (R * T) * N_CO2 * D_M_CO2 * M_CO2 * y_CO2  # [g/g*s] or [1/s]
        r_obs_film_intern_CO2 = 0  # [g/g*s] or [1/s]
        r_obs_film_H2O = 12 / (diameter ** 2 * density) * pressure / (R * T) * N_H2O * D_M_H2O * M_H2O * y_H2O  # [g/g*s] or [1/s]
        r_obs_film_intern_H2O = 0
        r_obs_film_O2 = 12 / (diameter ** 2 * density) * pressure / (R * T) * N_O2 * D_M_O2 * M_O2 * y_O2  # [g/g*s] or [1/s]
        r_obs_film_intern_O2 = 0  # [g/g*s] or [1/s]

        if np.all(r_obs_film_CO2 <= r_obs_pore_CO2):
            r_obs_CO2 = r_obs_film_CO2  # For mass transport limitation by film diffusion: Initialize!
            r_obs_intern_CO2 = r_obs_film_intern_CO2
            Regime_CO2 = 'Film'
        else:
            r_obs_CO2 = r_obs_pore_CO2  # For mass transport limitation by film diffusion: Initialize!
            r_obs_intern_CO2 = r_obs_pore_intern_CO2
            Regime_CO2 = 'Pore'

        # H2O
        if np.all(r_obs_film_H2O <= r_obs_pore_H2O):
            r_obs_H2O = r_obs_film_H2O  # For mass transport limitation by film diffusion: Initialize!
            r_obs_intern_H2O = r_obs_film_intern_H2O
            Regime_H2O = 'Film'
        else:
            r_obs_H2O = r_obs_pore_H2O  # For mass transport limitation by film diffusion: Initialize!
            r_obs_intern_H2O = r_obs_pore_intern_H2O
            Regime_H2O = 'Pore'

        # O2
        if np.all(r_obs_film_O2 <= r_obs_pore_O2):
            r_obs_O2 = r_obs_film_O2 # For mass transport limitation by film diffusion: Initialize!
            r_obs_intern_O2 = r_obs_film_intern_O2
            Regime_O2 = 'Film'
        else:
            r_obs_O2 = r_obs_pore_O2 # For mass transport limitation by film diffusion: Initialize!
            r_obs_intern_O2 = r_obs_pore_intern_O2
            Regime_O2 = 'Pore'

        robs = [r_obs_CO2, r_obs_H2O, r_obs_O2]
        robsfilm = [r_obs_film_CO2, r_obs_film_H2O, r_obs_film_O2]
        robsin = [r_obs_intern_CO2, r_obs_intern_H2O, r_obs_intern_O2]
        robspore = [r_obs_pore_CO2, r_obs_pore_H2O, r_obs_pore_O2]
        Regime = [Regime_CO2, Regime_H2O, Regime_O2]

        return robs, robspore, robsfilm, robsin, Regime

# updateCharConv
# update delta X / dt and charConv of a given particle number at particle step j and reactor step i.

# INPUTS
#
# dYV: float. dYV of the particle. Example: RM.dYV(j,i).
# YV: float. YVtotal of the particle. Example: Parameters.Conditions.YVtotal.
# Tau: array of floats. Real residence time in the reactor.
# i: integer. discretization step (index).
# j: integer. Particle discretization step.

# OUTPUT
# dXdt: double. delta charConv / dt
# charConv: double. total char conversion

    def updateCharConv(self,dYV, YV, Tau, j, i):

        dXdt = (1 - self.charConv[j, i - 1]) * np.sum(self.robs[:, j, i]) * (min((dYV / YV), 1)) ** 2
        charConv = dXdt * (Tau[i] - Tau[i - 1]) + self.charConv[j, i - 1]
        if charConv >= 0.999:
            charConv = 1
        return dXdt, charConv

# updateParticleSize
# Computes several particle variables at particle step j and reactor step i.

# INPUTS
#
# Fpc: array of floats.
# VolatileYield: floats.
# rf: float.
# VolatileYieldMax: float. Parameters.Conditions['YVtotal'].
# YV_total: float. Parameters.Conditions['YVtotal'].
# j: integer. Discretization step of the particle
# i: integer. Discretization step of the reactor.

# Outputs
# d_Particle_post : float
# density_post : float
# Surface_daf : float
# Surface_in_daf : float
# r_pore : float
# Po : float
# Surface_dry : float
# Surface_in_dry : float
# Volatile_Carbon_ratio : float

    def updateParticleSize(self,Fpc, VolatileYield, rf, VolatileYieldMax, YV_total, j, i):

        # Assign variable names
        Density_ash = self.fuelParams['ashDensity']
        density_carbon = self.fuelParams['carbonDensity']
        S0_max = self.fuelParams['S0max']
        # from this self
        d_Particle = self.D[j, i - 1] * 1e-6
        d_AshParticle = self.dAshParticle[j] * 1E-6
        m_Carbon_init = self.mFixedCarbon[j]
        m_Carbon = (1 - self.charConv[j, i]) * m_Carbon_init
        m_Carbon_pre = (1 - self.charConv[j, i - 1]) * m_Carbon_init
        V_AshParticle = math.pi / 6 * d_AshParticle ** 3
        m_Ash = Density_ash * V_AshParticle
        m_Volatile = self.mVolatile[j] * (1 - VolatileYield / VolatileYieldMax)
        m_Particle = m_Carbon + m_Ash + m_Volatile
        Po_init = self.porosity[j, i - 1]

        if m_Carbon > 0:
            if sum(self.robs[:, j, i]) > 0:
                beta = 1 / 3 * (1 - (sum(self.eff[:, j, i] * self.robsin[:, j, i]) / sum(self.robs[:, j, i])))
                if beta > 1 / 3:
                    beta = 1 / 3
            else:
                beta = 0

            d_Particle_post = (d_AshParticle ** 3 + (d_Particle ** 3 - d_AshParticle ** 3) * (m_Carbon / m_Carbon_pre) ** (3 * beta)) ** (1 / 3)
            V_Particle_post = math.pi / 6 * d_Particle_post ** 3
            V_Void = V_Particle_post - m_Carbon / density_carbon - V_AshParticle
            Po = V_Void / (V_Void + m_Carbon / density_carbon)

            if Po > 0.999:
                Po = 0.999

            Surface_ex_daf = 6 * d_Particle_post ** 2 * rf / ((d_Particle_post ** 3 - d_AshParticle ** 3) * density_carbon) / 1000  # in m^2/g

            if VolatileYield < YV_total * 0.99:  # Surface_dry erreicht S0_Max bei YV_total
                Surface_dry = S0_max * (m_Carbon_init + m_Ash + self.mVolatile[j] * (1 - YV_total / VolatileYieldMax)) / (m_Carbon_init + m_Volatile + m_Ash)  # ohne Pore-Closing
                Surface_in_daf = Fpc[i] * (Surface_dry * (m_Carbon + m_Volatile + m_Ash) / (m_Carbon + m_Volatile) - Surface_ex_daf) * math.sqrt(math.log(1 - Po) / math.log(1 - Po_init))
            else:
                Surface_in_daf = Fpc[i] * self.surfaceIn[j, i - 1] / Fpc[i - 1] * math.sqrt(math.log(1 - Po) / math.log(1 - Po_init))

            if Surface_in_daf <= 0:
                Surface_in_daf = 0.0001
                Surface_daf = Surface_ex_daf
                Surface_ex_daf = Surface_daf - Surface_in_daf

            r_pore = -2 * math.log(1 - Po) * rf / (Surface_in_daf * 10 ** 3 * density_carbon)
            # L_pore=Surface_in_daf*10^3*m_Carbon/(2*pi()*r_pore*rf)
            d_Particle_post = d_Particle_post * 10 ** 6
            density_post = (m_Carbon + m_Ash + m_Volatile) / V_Particle_post
            Surface_daf = Surface_in_daf + Surface_ex_daf

        else:
            d_Particle_post = d_AshParticle * 10 ** 6
            Surface_in_daf = 0
            Surface_daf = 0
            density_post = Density_ash
            r_pore = 0.0001
            Po = 1

        Surface_dry = Surface_daf * (m_Carbon + m_Volatile) / m_Particle  # mit Pore-Closing
        Surface_in_dry = Surface_in_daf * (m_Carbon + m_Volatile) / m_Particle

        Volatile_Carbon_ratio = max(m_Volatile / m_Carbon, 0)  # vermeide 0/0=NaN

        return d_Particle_post, density_post, Surface_daf, Surface_in_daf, r_pore, Po, Surface_dry, Surface_in_dry, Volatile_Carbon_ratio

# initParticleDev
# initializes the step 0 of particledevelopment.
#
# computes the particle size distribution and updates the obj properties accordingly. adapted from the function Partikelgroessenverteilung_Rosin_Rammler.m.
#
# INPUTS
# fuelParams: GlobalProperties.Fuel dictionary
# conditions: float. GlobalProperties.Conditions dictionary
#
# OUTPUT
# None

    def initParticleDev(self, fuelParams, conditions):

        d_RR_1 = fuelParams['d_RR_1']
        n_RR_1 = fuelParams['n_RR_1']
        d_min_1 = fuelParams['d_min_1']
        d_max_1 = fuelParams['d_max_1']

        AbbruchKriterium = 0.0001  # Change of Q3 less than 0.0001
        if d_min_1 == 0:
            d_min_1 = 1

        Partikelgroessenverteilung_Matrix = np.zeros((self.discretization + 1 , 2))  # preallocate

        if d_max_1 != 0:
            Partikelgroessenverteilung_Matrix[:, 0] = np.linspace(d_min_1, d_max_1, self.discretization + 1 )
            Partikelgroessenverteilung_Matrix[:, 1] = 100 - np.exp(-((Partikelgroessenverteilung_Matrix[:, 0] / d_RR_1) ** n_RR_1)) * 100

        else:
            Partikelgroessenverteilung_Matrix[0, 0] = d_min_1
            Partikelgroessenverteilung_Matrix[0, 1] = 100 - np.exp(-((Partikelgroessenverteilung_Matrix[0, 0] / d_RR_1) ** n_RR_1)) * 100

            for i in range(1, self.discretization + 1):
                Partikelgroessenverteilung_Matrix[i, 0] = Partikelgroessenverteilung_Matrix[i - 1, 0] * 1.1
                Partikelgroessenverteilung_Matrix[i, 1] = 100 - np.exp(-((Partikelgroessenverteilung_Matrix[i, 0] / d_RR_1) ** n_RR_1)) * 100
                if ((Partikelgroessenverteilung_Matrix[i, 1] - Partikelgroessenverteilung_Matrix[i - 1, 1]) /Partikelgroessenverteilung_Matrix[i - 1, 1]) < AbbruchKriterium:
                    break
        Partikelgroessenverteilung_Matrix[self.discretization, 1] = 100

        # Creation of "Particle size distribution_discrete_volume" with quantity type
        # "Volume

        x_u = np.flipud(Partikelgroessenverteilung_Matrix[:, 0])
        x_u = (x_u[:-1])
        x_o = np.flipud(Partikelgroessenverteilung_Matrix[:-1, 0])
        Q3_u = np.flipud(Partikelgroessenverteilung_Matrix[:, 1])
        Q3_u = (Q3_u[:-1])
        Q3_o = np.flipud(Partikelgroessenverteilung_Matrix[:-1, 1])
        cumDiam = (x_u + x_o) / 2
        v = Q3_u - Q3_o

        sumProduct = cumDiam ** (-3) * v
        qty = sumProduct / np.sum(sumProduct) * 100

        # Extension of the array for further parameters
        self.D = np.zeros((self.discretization, self.steps))
        self.D[:, 0] = cumDiam  # 'Diameter [10^-6 m]'
        self.vFrac = np.zeros((self.discretization, self.steps))
        self.vFrac[:, 0] = v  # 'Volume fraction [%]'
        # self.xFrac = fuelParams.density  # 'Mass fraction [%]

        self.surface = np.zeros((self.discretization, self.steps))  # 'Surface [m²/g]'
        self.rho = fuelParams['dryDensity'] * np.ones((self.discretization, self.steps))  # 'Density [kg/m³]'
        self.eff = np.ones((3, self.discretization, self.steps))  # [CO2, H2O, O2]  ->  'EffCO2 [-]'
        self.dXdt = np.zeros((self.discretization, self.steps))  # 1/s
        self.charConv = np.zeros((self.discretization, self.steps))  # [-]
        self.robspore = np.zeros((3, self.discretization, self.steps))  # [CO2, H2O, O2] -> 'r_obs_pore_O2 [g/g*s]'
        self.robsfilm = np.zeros((3, self.discretization, self.steps))  # idem ->  'r_obs_film_CO2 [g/g*s]'
        self.robs = np.zeros((3, self.discretization, self.steps))  # idem ->  'r_obs_H2O [g/g*s]'
        self.regime = np.empty((3, self.discretization, self.steps), dtype=object)  # [CO2, H2O, O2] ; cell array
        self.surfaceIn = np.zeros((self.discretization, self.steps))  # Surface_in_init_1_daf [m^2/g]'
        self.robsin = np.zeros((3, self.discretization, self.steps))  # 'r_obs_intern_CO2 [g/g*s]'
        self.Deff = np.zeros((3, self.discretization, self.steps))  # [CO2, H2O, O2] -> 'D_eff_CO2'
        self.DM = np.zeros((3, self.discretization, self.steps))  # [CO2, H2O, O2] -> 'D_M_CO2'
        self.rpore = np.zeros((self.discretization, self.steps))  # 'r_pore'
        self.porosity = np.zeros((self.discretization, self.steps))  # 'Porosity_1'
        self.surfaceDry = np.zeros((self.discretization, self.steps))  # 'Surface_dry'
        self.surfaceInDry = np.zeros((self.discretization, self.steps))  # 'Surface_in_dry'

        self.vol2C = np.zeros((self.discretization, self.steps))  # 'Volatile/Carbon'
        self.vol2C[:,0] =  (conditions['YVtotal'] / (1 - conditions['YVtotal']) * np.ones((self.discretization, 1))).flatten()  # from line 432

        self.xFrac = np.zeros((self.discretization, self.steps))
        Mass_fraction = v * fuelParams['density'] / (np.sum(v * fuelParams['density']) / 100) * np.ones((self.discretization, 1)) # 'Mass fraction [%]' from line 392 particle size distribution
        self.xFrac[:, 0]= (Mass_fraction[0,:])
        self.qtyFrac = qty  #particle_size_distribution_discrete_number_2{1,2}='Quantity fraction [%]'

        # Run particlestructure initialization and allocate to the variables.
        PS = self.initParticleStructure(self.discretization, fuelParams, conditions)

        self.dAshParticle = PS[:, 5]  # {1,5} = d_AshParticle_init_1 [10^6m]'
        self.mFixedCarbon = PS[:, 6]  # ParticleStructure_init{1,6} = 'm_FixedCarbon_init_1 [kg]';
        self.mVolatile = PS[:, 16]  # m_Volatile_init [kg]

        self.porosity = np.zeros((self.discretization, self.steps))  # Po_init [-]
        self.porosity[:, 0] = PS[:, 8]
        self.surfaceIn[:, 0] = PS[:, 10]
        self.surface[:, 0] = PS[:, 10] + PS[:, 9]
        self.surfaceDry[:, 0] = PS[:, 14]
        self.surfaceInDry[:, 0] = PS[:, 15]
        self.rpore[:, 0] = PS[:, 11]
        self.porosity[:, 0] = PS[:, 8]

# initParticleStructure
# This fucntion initializes the values and variables which are relevant.

#  INPUT
#
# discretization : number of the particle conditions: fuelStruct of the fuel.
# structFuel: GlobalProperties.Fuel dictionary
# conditions: float. GlobalProperties.Conditions dictionary
#
# OUTPUT
# Returns a matrix which is equivalent to ParticleStructure_init.

    def initParticleStructure(self, discretization, structFuel, conditions):

        # preallocate variables
        PS = np.zeros((discretization, 17))
        # compute some parameters
        VolatileYieldMax = conditions['YVtotal']
        Carbon_ar = (1 - VolatileYieldMax) * (1 - structFuel['ashContent'] / 100 - structFuel['moisture'] / 100) # with ash and moisture (ar=As-received)
        Volatile_ar = VolatileYieldMax * (1 - structFuel['ashContent'] / 100 - structFuel['moisture'] / 100) # with ash and moisture (ar=As-received) % remaining volatile after pyrolysis, based on ar-mass

        d_Particle_init = self.D[:, 0] * 10 ** -6  # in m
        V_Particle_init = np.pi / 6 * d_Particle_init ** 3  # in m^3
        m_Particle_init = structFuel['dryDensity'] * V_Particle_init  # in kg
        m_Ash_init = m_Particle_init * (structFuel['ashContent'] / 100) / (Carbon_ar + Volatile_ar + structFuel['ashContent'] / 100)  # = m_Ash=const
        d_AshParticle = (6 * m_Ash_init / (np.pi * structFuel['ashDensity'])) ** (1 / 3)  # = d_Ash = constant (no reaction of the ash)
        m_Carbon_init = m_Particle_init * Carbon_ar / (Carbon_ar + Volatile_ar + structFuel['ashContent'] / 100)  # = m_Carbon=const during pyrolysis
        m_Volatile_init = m_Particle_init * Volatile_ar / (Carbon_ar + Volatile_ar + structFuel['ashContent'] / 100)
        Char_AshContent_init = (structFuel['ashContent'] / 100) / ((structFuel['ashContent'] / 100) + Carbon_ar + Volatile_ar)
        V_Void_init = V_Particle_init - m_Carbon_init / structFuel['carbonDensity'] - m_Ash_init / structFuel['ashDensity']  # refers to carbon only
        Po_init = V_Void_init / (V_Void_init + (m_Carbon_init / structFuel['carbonDensity']))  # refers to carbon only
        Surface_init_dry = structFuel['S0max'] * (m_Carbon_init + m_Ash_init) / m_Particle_init  # m^2/g
        Surface_ex_init_daf = 6 * d_Particle_init ** 2 * conditions['rf'] / ((d_Particle_init ** 3 - d_AshParticle ** 3) * structFuel['carbonDensity']) / 1000  # in m^2/g
        Surface_in_init_daf = Surface_init_dry * m_Particle_init / (m_Carbon_init + m_Volatile_init) - Surface_ex_init_daf
        r_pore_init = -2 * np.log(1 - Po_init) * conditions['rf'] / (Surface_in_init_daf * structFuel['carbonDensity'])
        L_pore_init = Surface_in_init_daf * m_Carbon_init / (2 * np.pi * r_pore_init)

        for j in range(discretization):
            if Surface_in_init_daf[j] <= 0:
                Surface_ex_init_daf[j] = Surface_ex_init_daf[j] - abs(Surface_in_init_daf[j])
                Surface_in_init_daf[j] = 0.0001
                r_pore_init[j] = 0.0001
                L_pore_init[j] = Surface_in_init_daf[j] * 10 ** 3 * m_Carbon_init[j] / (2 * np.pi * r_pore_init[j])

        Surface_in_init_dry = Surface_in_init_daf * (m_Carbon_init + m_Volatile_init) / m_Particle_init

        # store values in the preallocated matrix
        PS[:, 0] = d_Particle_init * 10 ** 6
        PS[:, 1] = V_Particle_init
        PS[:, 2] = m_Particle_init
        PS[:, 3] = Char_AshContent_init
        PS[:, 4] = m_Ash_init
        PS[:, 5] = d_AshParticle * 10 ** 6
        PS[:, 6] = m_Carbon_init
        PS[:, 7] = V_Void_init
        PS[:, 8] = Po_init
        PS[:, 9] = Surface_ex_init_daf
        PS[:, 10] = Surface_in_init_daf
        PS[:, 11] = r_pore_init
        PS[:, 12] = L_pore_init
        PS[:, 13] = structFuel['dryDensity'] * np.ones(discretization)
        PS[:, 14] = Surface_init_dry
        PS[:, 15] = Surface_in_init_dry
        PS[:, 16] = m_Volatile_init

        return PS















