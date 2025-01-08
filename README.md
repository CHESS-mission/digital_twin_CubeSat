# A Digital Twin of the CHESS CubeSat

This repository hosts a Python-based simulation framework developed to support the mission design and operations of the CHESS CubeSat. It has two main purposes:

### Purposes

1. *Simulation tool for mission design*: Provides tools for modeling and simulating different subsystems of the CHESS CubeSat to inform decisions and help in the development, testing, and validation tasks.

2. *Foundation for a digital twin*: Serves as the groundwork for a digital twin capable of simulating real-time satellite operations, including anomaly detection.

This framework was initially developed as part of a semester project with the EPFL Spacecraft Team in collaboration with the EPFL Space Center.

## Table of Contents
- [Libraries and Requirements](#libraries_and_requirements)
- [Repository Structure](#repository_structure)
- [Usage](#usage)
- [Test Case](#test_case)
- [Documentation](#documentation)
- [Limitations and Future Work](#limitations_and_future_work)
- [Authors and Acknowledgement](#authors_acknowledgement)

## Libraries and Requirements

This application runs with Python version 3.9 or above. All dependencies are listed in the file *digital_twin_env.yaml*. This file can be used to initiate the appropriate Conda environment and install necessary libraries with the command on the terminal, at the root of the repository: 

```bash
conda env create -f digital_twin_env.yml
```

The created Conda environment is named *digital_twin_env*. To activate it, use the following command:

```bash
conda activate digital_twin_env
```

You are ready to run a simulation!

## Repository Structure

```bash
├── README.md
├── data
│   ├── atmosphere_data					# Files for atmosphere models
│   ├── ground_station					# Config files for ground stations
│   ├── mission_design					# Config files for mission design parameters
│   ├── orbit							# Config files for orbit parameters
│   ├── simulation						# Config files for simulation parameters
│   └── spacecraft						# Config files for spacecraft parameters
├── digital_twin_env.yml 				# Environment configuration file for setting up dependencies
├── docs								# Documentation
│   ├── UML_diagram.png					# UML class diagram of the framework
│   ├── doc_generation					# Folder used to generate the HTML documentation
│   ├── html							# HTML documentation
│   └── parameters.xlsx					# Default parameters sources
├── results								# Stores results
│   ├── data							# Stores output numpy files (can be used for more visualizations with jupyter notebooks)
│   └── figures							# Stores generated figures
└── src									# Source code
    ├── digital_twin
    │   ├── constants.py				# General constants for the simulation
    │   ├── ground_station				# Ground station module
    │   │   └── ground_station.py
    │   ├── mode_switch.py				# Implementation of the mode switch decision tree
    │   ├── orbit_propagator			# Orbit propagator module
    │   │   ├── atmosphere_model.py		# Implementation of the atmosphere models
    │   │   ├── constants.py			# Constants used for propagation
    │   │   └── orbit_propagator.py		# Implementation of the propagator
    │   ├── plotting.py					# Functions used to plot results
    │   ├── report.py					# Functions used to generate outputs
    │   ├── simulation.py				# Implementation of the Simulation class
    │   ├── spacecraft					# Spacecraft module
    │   │   ├── adcs					# ADCS module
    │   │   │   └── adcs.py
    │   │   ├── eps						# EPS module
    │   │   │   └── eps.py
    │   │   ├── obc						# OBC module
    │   │   │   └── obc.py
    │   │   ├── payload					# Payload module
    │   │   │   └── payload.py
    │   │   ├── spacecraft.py			# Implementation of the spacecraft class
    │   │   ├── subsystem.py			# Subsystem interface
    │   │   └── telecom					# Telecom module
    │   │       └── telecom.py
    │   └── utils.py					# Utils functions
    ├── main.py						# Main file to run
    └── notebooks					# Additional notebooks used for plots, experiments, and analysis of results
```

### Notes

- *Notebooks/* folder: description of notebooks are provided at the beginning of these notebooks. They were generally used to either produce more complicated plots which plotting functions are not general enough to add them in the framework or do additional analysis, or to do experiments and test different functions that are suibsequently added to the framework.

- *data/atmosphere_data/* folder: this folder does not need to be updated. It contrains data for the **NRLMSISE00**, **JB2008** and **solar_activity** atmosphere models. The data for the first two are frequently updated by the library using them (*ATMOS*) at the beginning of the simulations. Regarding the $3^{rd}$ model, the data does not need to be updated. It is solar activity data ($F10.7$ index and $Ap$) which were exrtracted from the software DRAMA.

- *__init__.py* files are not written in the tree for clarity. However they exist in every folder that constitutes a module.

- *parameters.xlsx* file: this excel document gathers descriptions and sources for the default parameters for the CHESS CubeSat provided in the configuration files. If this framework has to be used for another satellite, those will need to be updated.

- *mode_switch.py*: the mode switch decision tree used for the implementation of the mode switch algorithm in the class `ModeSwitch` can be found in the Appendix of the report located in the *docs/* folder.

- The organization of the code as well as the main design choices are explained in details in the report located in the *docs/* folder.


## Usage

### Input
User parameters are divided into five categories, each with its own JSON configuration file:

- **Spacecraft**: Defines general parameters like mass, cross-sectional area, initial operating mode, and subsystem-specific details.
- **Orbit**: Specifies initial orbital elements, epoch, and orbit type (e.g., SSO).
- **Simulation**: Includes general simulation settings, such as duration, timestep, and the atmospheric model to use.
- **Ground station**: Lists ground stations with their respective locations and elevation angles.
- **Mission design**: Gathers parameters to generate a report (e.g., data to save, figures to generate) and additional user input, such as commands to initiate SAFE mode at specific times.

Templates filled with default values for the configuration files can be found in the folder *data/* and subfolders *spacecraft/*, *orbit/*, *simulation/*, *ground_station/* and *mission_design*. It is important to place each of the 5 config files listed above in their respective folder.

If for any reason want to change where to place them, change the path constants in "src/constands.py".

By convention all the values provided in the  in *data/* folder are precised in kilograms ($kg$), seconds ($s$), kilometers ($km$) and degrees ($deg$), Watts ($W$), or any combination of those. However, there are some exceptions:
- in the config files, when the used is asked to specified the unit for a specific parameters, they are free to use any. This occurs in ..... show
- the consumption rates are expressed in $Watt-Hour ($Wh$) and not in Watt seconds. This is because they are provided with this unit in the EPFL Spacecraft Team power budget.
- solar cell area in $m^2$ (because otherwise it is too small)

=> check if there is any other!

Arguments in the config files
- arguments written with init are specific to the spacecraft state 



### Running a Simulation

In order to run a simulation, follow these steps:

1. Place the 5 configuration files in their respective folder
2. Activate the conda environment
3. At the root directory run:

```bash
python3 -W"ignore" src/main.py simulation_template.json orbit_template.json spacecraft_template.json ground_station_template.json mission_design_template.json

```

where the example was provided with the files with default values.

### Flow

Once the above command has been executed, the following steps are performed:
- If running the application for the first time or after a long gap, it automatically downloads the updated files required by the `astropy.utils.iers` module (tables with Earth orientation data). Otherwise, it displays a message indicating that the tables are already updated
- if the parameter `"verbose"` is set to `"yes"` in the "simulation.json" configuration file, it displays the initialization of classes and the timesteps every 1000 timestep and the time it took to run.
- if the parameter `"print_parameters"` is set to `"yes"` in the "simulation.json" configuration file, displays initial parameters before starting the main simulation loop.
- run the simulation
- add results in the folder provided in the "mission_design.json" configuration file.


### Output

The output of the simulation depends on the parameters provided in the "mission_design.json" configuration file
=> explain what all of those parameters mean!!!
explain what kind of output can be given!

possibility to save state!

2 kinds of outputs can be saved: figures, and data


## Test Case

For more background about this test case, refer to the report.



## Documentation

Detailed information about each function and class can be found in the **documentation**. To access it, open the file "index.html" in the folder *docs/html/*, where you can explore the available classes, functions, and modules

If documentation needs to be generated again, follow these commands in the *build/* folder:

   ```bash
   cmake .. -DDOCUMENTATION=0N
   make
   ```


doc_generation

In addition, the default parameters for the config files are in the excel file


## Limitations and Future Work

Design assumptions and limitations are documented in detail within *docs/* in the file "report.pdf". This document needs to be updated when more is added to the project.

The main design choices, as well as the assumptions and limitations of this application are detailed in the report, which can be found in the `docs/` folder.

## Authors and Acknowledgement

### Initial development (fall 2024)

Author:
- Mathilde Simoni, MSc in Computational Science and Engineering at EPFL

Supervisors:
- Andrew Price
- Mathieu Udriot

Referent Professor:
- Jean-Paul Kneib

Future authors are encouraged to add their names and details as they contribute to this  project.
