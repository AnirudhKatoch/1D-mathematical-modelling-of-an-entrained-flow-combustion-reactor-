# 1D mathematical modelling of an entrained flow combustion reactor






FLow chart - https://drive.google.com/file/d/1ZzINDuxjELkepz-B3d3XgW1K6srqNra2/view?usp=sharing

Tasks carried out

* Gained expertise in object-oriented programming and mathematical modeling on Python and MATLAB.
* Skillfully extracted and analyzed complex Excel data using Pandas, SciPy and NumPy in Python.
* Designed, automated, and simulated a Python-based mathematical model for an entrained flow combustor.
* Successfully integrated the developed model with REFPROP to enhance fluid behavior modeling and analysis.


## Class definitions

`GlobalParameters.py` is a Python class designed for entrained flow combustion reactor simulations. It provides essential global parameters, ensuring accuracy and consistency in the simulations. The class follows object-oriented principles, encapsulating properties like `Fuel1`, `Fuel2`, `Geometry`, `Conditions`, and `n0`. Well-documented and error-free, it allows for easy extension and maintenance, serving as a fundamental component of the simulation.

`GasProperties.py` is a Python class crucial for entrained flow combustion reactor simulations, providing access to essential gas properties at specific temperature and pressure conditions. It encapsulates functions to retrieve properties such as heat capacity, heat conductivity, dynamic viscosity, Prandtl number, and enthalpy for various gases like methane, hydrogen, nitrogen, and more. This modular class follows object-oriented principles, making it easy to add new gas properties. It is well-documented, error-free, and seamlessly integrates into the reactor simulation framework.

`MixtureProperties.py` is a Python class for entrained flow combustion reactor simulations, offering gas mixture property access at specified temperature and pressure using REFPROP data. It inherits a modular, object-oriented design, functioning as a subclass of `GasProperties.py` for property calculations like heat capacity, heat conductivity, dynamic viscosity, Prandtl number, and enthalpy. It can easily integrate new properties, maintains well-documented error-free code, and ensures seamless interaction with `GasProperties.py`, enhancing adaptability and expandability for future simulation needs.

`ParticleDevelopment.py` is a Python class integral to entrained flow gasification reactor simulations, focusing on particle behavior modeling. Key features include property arrays for simulation, functions for particle development aspects, and the creation of a `fuelParams` dictionary for storing fuel parameters. The class's integration with reactor discretization and object-oriented structure makes it flexible and extensible. However, some issues, like missing files impacting "MOLECULAR DIFFUSION" calculations and a mathematical problem with matrix dimensions, need attention for improved functionality. This class offers valuable insights into particle dynamics within the reactor simulation framework.

The `ReactorMesh.py` class is a vital part of a Python-based reactor simulation. It performs essential calculations by integrating data from other classes. This object-oriented class creates a one-dimensional mesh with user-defined steps. However, it faces challenges due to missing properties and dependencies, which need resolution for accurate simulation results. Clear documentation and comments within the code are essential for understanding its functionality.

The existing Python code, 'main.py,' is structured to simulate a reactor model. It comprises classes for global parameters, particle development, reactor mesh, gas properties, and mixture properties. The code introduces functions to compute the swelling ratio and 'VM.' There are comments denoting potential issues. A significant overhaul is needed to accommodate new insights on ignition, combustion, and the transition from pyrolysis to combustion reactions within the reactor. The code serves as the foundation for this updated simulation, and it requires restructuring and rewriting to align with the latest understanding of the reactor's behavior.




