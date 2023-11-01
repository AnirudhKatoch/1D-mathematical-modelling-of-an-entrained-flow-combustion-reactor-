# GlobalParameters

# Note
# Inspired from the MATLAB class GlobalParameters  written by  authors Netter, Tobias (tobias.netter@tum.de) ,
# Ceruti, Amedeo (amedeo.ceruti@tum.de) within the project "Entrained Flow Gasification Reactor Modeling with MATLAB"
#
# Class containing global parameters needed for the entrained flow combustion Reactor.
# This class stores diverse parameters and is input in everal functions to provide needed parameters to the model.
#
# Example:
# Parameters = GlobalParameters('Parameters\db.xlsx','Parameters\reactor_geometry.xlsx', 1,'Parameters\reactor_conditions.xlsx')
#
# Contents
#
# Class properties
# GlobalParameters
# setfueldictionary
# computedensityoffuel
# setreactorconditions
# setreactorgeometry
# computeYV
# computevolumetricflow
# getInputMixture
# molarflow
# computeResidenceTime

from GasProperties import GasProperties # Another class  required to run fucntion molarflow in GLobal Paramters class
import pandas as pd
import numpy as np

class GlobalParameters:

    Fuel1 = None # dictionary containing parameters of fuels, set with setfueldictionary
    Fuel2 = None # dictionary containing parameters of fuels, set with setfueldictionary
    Geometry = None # dictionary containing reactor geometry paramters, setreactorgeometry
    Conditions = None # dictionary containing reaction conditions, setreactorconditions
    n0 = None # initial molar flows entering the reaction in an array
    R = 8.3145 # gas constant

# Construct an object of the GlobalParameters class.
#
# INPUTS
#
# pathToFuelExcelFile: char array containing the path to the database with fuel parameters. Example: 'Parameters\db_combustion_reactor.xlsx'
# pathToReactorExcelFile: Path to the excel file containing the reactor geometry information. Example: 'Parameters\reactor_geometry_combustion_reactor.xlsx'
# caseNumber: integer of wanted case (row) in the excel file that is inout fuel , can be seen in excel file reactor_conditions_combustion_reactor
# pathToConditionsExcelFile: char array containing the path to the excel with fuel parameters. Example: 'Parameters\reactor_conditions_combustion_reactor.xlsx'.

    def __init__(self,pathToFuelExcelFile, pathToReactorExcelFile, caseNumber, pathToConditionsExcelFile):
        import pandas as pd
        args = [0] * 6
        if isinstance(pathToFuelExcelFile, str):
            args[0] = pathToFuelExcelFile
        else:
            raise ValueError('pathToFuelExcelFile is not a string.')

        if isinstance(pathToConditionsExcelFile, str):
            args[3] = pathToConditionsExcelFile

            if isinstance(caseNumber, int):
                args[4] = caseNumber
                df = pd.read_excel(args[3])
                requiredName = df.iloc[args[4] - 1]
                args[1] = requiredName['Name1']
                args[2] = requiredName['Name2']
            else:
                raise ValueError('caseNumber is not an integer.')
        else:
            raise ValueError('pathToReactorExcelFile is not a string.')

        if isinstance(pathToReactorExcelFile, str):
            args[5] = pathToReactorExcelFile
        else:
            raise ValueError('pathToReactorExcelFile is not a string.')

        self.Fuel1 = self.setfueldictionary(args[0], args[1])  # construct Fuel dictionary
        self.Fuel2 = self.setfueldictionary(args[0], args[2])  # construct Fuel dictionary
        self.Conditions = self.setreactorconditions(args[3], args[4])  # construct Reactor conditions dictionary
        self.Geometry = self.setreactorgeometry(args[5])  # construct Reactor geometry dictionary
        self.n0 = self.molarflow(298.15, 101325)
        self.Conditions['YVtotal'], self.Conditions['YVtotal2'], self.Conditions['YV'] = self.computeYV()

# setfueldictionary
# Sets a dictionary containing all parameters in an excel file with an appropioate format.
#
# INPUTS
# path: Path to table.
# fuelType: The fuel type should be as a string. Example: 'B-1'.
#
# OUTPUT
# s: Structured dictionary

    def setfueldictionary(self, path, fuelType):
        import pandas as pd
        # Read the excel sheets
        df_coal = pd.read_excel(path, sheet_name='Parameters Coal')
        df_char = pd.read_excel(path, sheet_name='Parameters Char')
        df_ash = pd.read_excel(path, sheet_name='Parameters Ash')

        parameters_coal = df_coal.iloc[1:, 1].tolist()
        column_name_coal = fuelType
        if column_name_coal in df_coal.columns:
            value_coal = df_coal.loc[2, column_name_coal]
            values_coal = df_coal[column_name_coal].iloc[1:].tolist()
            dictionary_coal = dict(zip(parameters_coal, values_coal))

        keys_to_change_coal = {
            'Moisture, %wt (ar)': 'moisture',
            'Ash, %wt (wf)': 'ashContent',
            'Volatile Matter, %wt (wf)': 'volatiles',
            'Fixed Carbon, %wt (wf)': 'fixedC',
            'Carbon, %wt (wf)': 'elementalC',
            'Hydrogen, %wt (wf)': 'elementalH',
            'Nitrogen, %wt (wf)': 'elementalN',
            'Sulfur, %wt (wf)': 'elementalS',
            'ϑ, -': 'Phi',
            'ρ, -': 'rho',
            'T_set, °C': 'Tset',
            'p_set': 'pset',
            'X_Tset_pset, -': 'XTsetpset',
            'X_max, -': 'Xmax',
            'T_max, °C': 'Tmax',
            'd (Rosin Rammler), m': 'd_RR_1',
            'n (Rosin Rammler), -': 'n_RR_1',
            'd_10, m': 'd_min_1',
            'd_90, m': 'd_max_1',
            'Swelling Factor (SF), -': 'SW',
        }

        s = {}

        for key, value in dictionary_coal.items():
            if key in keys_to_change_coal:
                new_key = keys_to_change_coal[key]
                s[new_key] = value

        s['k0_Pyrolysis'] = 293  # kJ/mol
        s['Ea_Pyrolysis'] = 51
        s['k_01'] = 200000  # De Young page 37
        s['E_A1'] = 104.7
        s['a_1'] = 0.56  # Volatile Matter, %wt (wf)
        s['k_02'] = 1.3E7
        s['E_A2'] = 167.5
        s['a_2'] = 1

        s['d_RR_1'] *= 1E6
        s['d_min_1'] *= 1E6
        s['d_max_1'] *= 1E6

        parameters_char = df_char.iloc[1:, 1].tolist()
        column_name_char = fuelType
        if column_name_char in df_char.columns:
            value_char = df_char.loc[1, column_name_char]
            values_char = df_char[column_name_char].iloc[1:].tolist()
            dictionary_char = dict(zip(parameters_char, values_char))

        keys_to_change_char = {
            'E_A_CO2, J/kmol required': 'Ea_CO2',
            'k_0_CO2, kg_C/(m2 s Pa^n) required': 'k0_CO2',
            'n_CO2, required': 'n_CO2',
            'E_A_O2, J/kmol required': 'Ea_O2',
            'k_0_O2, kg_C/(m2 s Pa^n) required': 'k0_O2',
            'n_O2, required': 'n_O2',
            'E_A_H2O, J/kmol required': 'Ea_H2O',
            'k_0_H2O, kg_C/(m2 s Pa^n) required': 'k0_H2O',
            'n_H2O, required': 'n_H2O',
            'A_max, annealing': 'A_max_Annealing',
            'k_0, 1/s annealing ': 'k0_Annealing',
            'E_A, J/kmol annealing': 'EA_Annealing',
            'S_0, m2/g (waf) pc': 'S0max',
            'k_0, 1/s pc ': 'k0_pc',
            'E_A, J/kmol pc ': 'EA_pc',
            'A_min, pc ': 'A_min',
        }

        for key, value in dictionary_char.items():
            if key in keys_to_change_char:
                new_key = keys_to_change_char[key]
                s[new_key] = value

        s['Ea_CO2'] *= 1E-6
        s['Ea_O2'] *= 1E-6
        s['Ea_H2O'] *= 1E-6
        s['EA_Annealing'] *= 1E-6
        s['EA_pc'] *= 1E-6
        s['carbonDensity'] = 2270  # constant

        parameters_ash = df_ash.iloc[1:, 1].tolist()
        column_name_ash = fuelType
        if column_name_ash in df_ash.columns:
            value_ash = df_ash.loc[1, column_name_ash]
            values_ash = df_ash[column_name_ash].iloc[1:].tolist()
            dictionary_ash = dict(zip(parameters_ash, values_ash))

        keys_to_change_char = {
            'Na2O, %wt required': 'Na2O',
            'MgO, %wt required': 'MgO',
            'Al2O3, %wt required': 'Al2O3',
            'SiO2, %wt required': 'SiO2',
            'P2O5, %wt required': 'P2O5',
            'SO3, %wt required': 'SO3',
            'K2O, %wt required': 'K2O',
            'CaO, %wt required': 'CaO',
            'TiO2, %wt required': 'TiO2',
            'MnO, %wt required': 'MnO',
            'Fe2O3, %wt required': 'Fe2O3',
        }

        for key, value in dictionary_ash.items():
            if key in keys_to_change_char:
                new_key = keys_to_change_char[key]
                s[new_key] = value

        s['Group_1'] = 0
        s['Group_2'] = 0
        s['NiO'] = 0
        s['CuO'] = 0
        s['ZnO'] = 0
        s['SrO'] = 0
        s['PbO'] = 0
        s['V2O5'] = 0
        s['MnO'] = 0
        s['cp'] = 1680  # constant
        s['HCV'] = 21.98E6  # HHV

        # Additional stuff
        s['Name'] = fuelType
        s['Group_3'] = 100 - s['Group_1'] - s['Group_2']
        s['ashDensity'] = self.computedensityoffuel(s)
        # Elemental Composition
        s['elementalO'] = 100 - (s['elementalC'] + s['elementalH'] + s['elementalN'] + s['elementalS'])

        # fitting Surface_dry S0_max für TWB (?)
        if s['Name'] == "TWB":
            s['S0max'] = 320
            s['dryDensity'] = 1370
            s['density'] = 1460
        elif s['Name'] == "B-5":
            s['dryDensity'] = 1510
            s['density'] = 1670
        elif s['Name'] == "L-3":
            s['dryDensity'] = 1620
            s['density'] = 1820
        else:
            s['dryDensity'] = 1450
            s['density'] = 1650

        s['n'] = [s['n_CO2'], s['n_H2O'], s['n_O2']]  # reaction orders
        s['k0'] = [s['k0_CO2'], s['k0_H2O'], s['k0_O2']]  # reaction orders
        s['Ea_i'] = [s['Ea_CO2'], s['Ea_H2O'], s['Ea_O2']]  # activation energies of each component

        return s

#computedensityoffuel
#Computes density of the fuel based on composition.
#
#INPUTS
#fueldictionary: dictionary containing the fuel information.
#
#OUTPUT
#density: float. fuel density.

    def computedensityoffuel(self,fuelDictionary):
        denom = (fuelDictionary['Na2O'] / 2.27 + fuelDictionary['MgO'] / 3.58 +
                 fuelDictionary['Al2O3'] / 3.95 + fuelDictionary['SiO2'] / 2.65 +
                 fuelDictionary['P2O5'] / 2.72 + fuelDictionary['SO3'] / 1.92 +
                 fuelDictionary['K2O'] / 2.35 + fuelDictionary['CaO'] / 3.34 +
                 fuelDictionary['TiO2'] / 4.23 + fuelDictionary['V2O5'] / 3.36 +
                 fuelDictionary['MnO'] / 5.37 + fuelDictionary['Fe2O3'] / 5.24 +
                 fuelDictionary['NiO'] / 6.67 + fuelDictionary['CuO'] / 6.32 +
                 fuelDictionary['ZnO'] / 5.61 + fuelDictionary['SrO'] / 4.7 +
                 fuelDictionary['PbO'] / 9.53)
        density = 1 / (denom / 100) * 1000
        return density

#setreactorconditions
#Setter for the Reactor property containing a dictionary with the conditions of the reactor.
#
#INPUTS
#excelTable: Excel table read from provided excel file path.
#caseNumber: integer of wanted case (row) in the excel file.
#
#OUTPUT
#s_set_reactor_conditions :  Structured dictionary

    def setreactorconditions(self, excelTable, caseNumber):

        # Read the Excel file
        df = pd.read_excel(excelTable)

        # Find rows with the specified number in the first column
        matching_rows = df[df.iloc[:, 0] == caseNumber]

        # Convert the matching rows to a list of lists
        rows_list = matching_rows.values.tolist()

        # Get the first row as a list
        first_row = df.columns.tolist()

        # Create a dictionary
        s_set_reactor_conditions = {}

        # Iterate over the first row elements
        for col_name in first_row:
            # Get the column index of the current element
            col_index = first_row.index(col_name)

            # Get the values from matching rows for the current column
            values = [row[col_index] for row in rows_list]

            # Remove square brackets and convert to float if only one value exists
            if len(values) == 1:
                values = values[0]

            # Add the values to the dictionary
            s_set_reactor_conditions[col_name] = values

        s_set_reactor_conditions['Tin'] = 25  #Degree celcius
        s_set_reactor_conditions['Tparticlein'] = 25 #Degree celcius
        s_set_reactor_conditions['Tmixing'] = 25 #Degree celcius

        return s_set_reactor_conditions

#setreactorgeometry
#Fucntions for setting geometry of the reactor for the property Geometry
#
#INPUTS
#excelTable: Excel table read from provided excel file path.
#
#OUTPUT
#s_setreactorgeometry:dictionary containing geometry parameters

    def setreactorgeometry(self,excelTable):
        import pandas as pd
        s_setreactorgeometry = {}

        # Read the Excel file into a DataFrame
        excel_data = pd.read_excel(excelTable)

        # Select columns with values
        data = excel_data.iloc[0, 1:].values

        # Convert to dictionary
        s_setreactorgeometry = {col_name: value for col_name, value in zip(excel_data.columns[1:], data)}

        # Calculate additional values
        # Calculate additional values
        s_setreactorgeometry['heatedTubeLength'] = s_setreactorgeometry['numHeatElems'] * s_setreactorgeometry['heightHeatElem']  # m :Heated pipe length
        s_setreactorgeometry['tubeDiameter'] = s_setreactorgeometry['innerRadius'] * 2  # m :Pipe inner diameter reaction zone

        s_setreactorgeometry['stepLengthPyrolysis'] = 0.001
        s_setreactorgeometry['pyrolysisLength'] = s_setreactorgeometry['tubeLength']
        s_setreactorgeometry['numPyrolysisSteps'] = round(s_setreactorgeometry['pyrolysisLength'] / s_setreactorgeometry['stepLengthPyrolysis'])


        s_setreactorgeometry['stepLengthCombustion'] = s_setreactorgeometry['stepLengthPyrolysis'] * 50
        s_setreactorgeometry['numCombustionSteps'] = s_setreactorgeometry['numPyrolysisSteps'] +\
                                                     round((s_setreactorgeometry['tubeLength'] - s_setreactorgeometry['pyrolysisLength']) / s_setreactorgeometry['stepLengthCombustion']);


        # Adjust conditions depending on external temp and fuel
        if (self.Fuel1['Name'] == "L-3" and self.Conditions['outerTemperatureTube'] == 1200) or \
                (self.Fuel1['Name'] == "T-1" and self.Conditions['outerTemperatureTube'] == 1400):
            s_setreactorgeometry['stepLengthPyrolysis'] = s_setreactorgeometry['stepLengthPyrolysis'] / 2
        elif (self.Fuel1['Name'] == "L-3" and self.Conditions['outerTemperatureTube'] == 1400):
            s_setreactorgeometry['stepLengthPyrolysis'] = s_setreactorgeometry['stepLengthPyrolysis'] / 4
        elif (self.Fuel1['Name'] == "L-3" and self.Conditions['outerTemperatureTube']) == 1600:
            s_setreactorgeometry['stepLengthPyrolysis'] = s_setreactorgeometry['stepLengthPyrolysis'] / 8

        return s_setreactorgeometry

#computeYV
# Compute maximum volatile conversion depending on a given pressure and temperature stored in the self.Conditions dictionary
#
# INPUTS
# None
#
# OUTPUT
# YV: double.Total volatilization conversion at T and p of reactor operation.
# X_V1: double.Total volatilization conversion at T and p of fuel 1.
# X_V2: double.Total volatilization conversion at T and p of fuel 2.

    def computeYV(self):
        import numpy as np
        T = self.Conditions['outerTemperatureTube'] + 273.15  # K
        p = self.Conditions['pressure'] * 1E5  # Pa

        Tset1 = self.Fuel1['Tset'] + 273.15
        Tset2 = self.Fuel2['Tset'] + 273.15
        pset1 = self.Fuel1['pset'] * 1E5
        pset2 = self.Fuel2['pset'] * 1E5

        X_Vp1 = self.Fuel1['XTsetpset'] / 100 - np.log(p / pset1) / self.Fuel1['rho']
        X_Vp2 = self.Fuel2['XTsetpset'] / 100 - np.log(p / pset2) / self.Fuel2['rho']

        X_VT1 = self.Fuel1['XTsetpset'] / 100 + (self.Fuel1['Xmax'] / 100 - self.Fuel1['XTsetpset'] / 100) * (
                    1 - np.exp(-self.Fuel1['Phi'] * (T - Tset1)))
        X_VT2 = self.Fuel2['XTsetpset'] / 100 + (self.Fuel2['Xmax'] / 100 - self.Fuel2['XTsetpset'] / 100) * (
                    1 - np.exp(-self.Fuel2['Phi'] * (T - Tset2)))

        X_V1 = self.Fuel1['XTsetpset'] / 100 * (X_Vp1 / (self.Fuel1['XTsetpset'] / 100)) * (
                    X_VT1 / (self.Fuel1['XTsetpset'] / 100))
        X_V2 = self.Fuel2['XTsetpset'] / 100 * (X_Vp2 / (self.Fuel1['XTsetpset'] / 100)) * (
                    X_VT2 / (self.Fuel1['XTsetpset'] / 100))

        YV = X_V1 * (1 - self.Conditions['w_2']) + X_V2 * self.Conditions['w_2']

        return YV, X_V1, X_V2

#computevolumetricflow
#Computes volumetric gas flow of the mixture of primary and secondary gases.
#
#INPUTS
#None
#
#Output
#gasflow: volumetric gas flow in a double and in m3/h.

    def computevolumetricflow(self):
        gasflow = self.Conditions['inputprimary'] + self.Conditions['inputsecondary']  # m3/h
        return gasflow

#getInputMixture
#Returns a MixtureProperties object at given conditions for the input components specified in Conditions.
#
#INPUTS
#T: float in K
#P: float in Pa
#
#OUTPUT
#mix: MixtureProperties object with the input components at T and P.

    def getInputMixture(self, T, P):
        n = self.molarflow(T, P)
        xi = n / sum(n)

        mixture1 = {'nitrogen': xi[0], 'oxygen': xi[1], 'hydrogen': xi[2], 'CO2': xi[3], 'water': xi[4]}
        mix = MixtureProperties(mixture1, T, P)
        return mix

#molarflow
#Converts input gas mass flows to mol/s computing its properties.
#
#INPUTS
#T:double. Temperature in K
#P: double. Pressure in Pa.
#
#OUTPUT
#n: ordered array of doubles. Molar flows in order [xN2, xO2, xH2, xCO2, xH2O] and in mol/s

    def molarflow(self, T, P):

        N2 = GasProperties('nitrogen')  # m3/h
        xN2 = ((self.Conditions['inputprimary'] + self.Conditions['inputsecondary']) * 0.78084) / 3600 * N2.computedensity(T, P) / N2.getMW()

        O2 = GasProperties('oxygen')  # m3/h
        xO2 = ((self.Conditions['inputprimary'] + self.Conditions['inputsecondary']) * 0.20946) / 3600 * O2.computedensity(T, P) / O2.getMW()

        H2 = GasProperties('hydrogen')  # m3/h
        xH2 = ((self.Conditions['inputprimary'] + self.Conditions['inputsecondary']) * 5e-7) / 3600 * H2.computedensity(T, P) / H2.getMW()

        CO2 = GasProperties('CO2')  # m3/h
        xCO2 = ((self.Conditions['inputprimary'] + self.Conditions['inputsecondary']) * 0.0004) / 3600 * CO2.computedensity(T, P) / CO2.getMW()

        H2O = GasProperties('water')  # m3/h
        xH2O = ((self.Conditions['inputprimary'] + self.Conditions['inputsecondary']) * 0.0092995) / 3600 * CO2.computedensity(T, P) / H2O.getMW()

        n = [xN2, xO2, xH2, xCO2, xH2O]
        return n

#computeResidenceTime
#computes gas residence time in the reactor in s.
#
#INPUT
#T: double. Temperature in K.
#P: double. Pressure in bar.
#molarFlow: double. total molar flow in mol/s.
#
#OUTPUT
#tau: double. Residence time in s.

    def computeResidenceTime(self, T, P, molarFlow):
        tau = self.Geometry['tubeLength'] * self.Geometry['tubeDiameter'] ** 2 * 0.25 * pi * P * 1E5 / (molarFlow * 8.314 * T)  # s
        return tau

