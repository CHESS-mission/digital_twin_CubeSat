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

- Python 3.9 or above (for advanced typing hints)
- All dependencies are listed in the file *environment.yaml* which can be used to initiate the appropriate Conda environment.


## Repository Structure
atmosphere data/ is something filled and automatically updated by the application usign a library

- solar activity needs to be added manually
- other 2 are updated automatically


## Usage

### User Parameters
User parameters are divided into five categories, each with its own JSON configuration file:

- **Spacecraft**: Defines general parameters like mass, cross-sectional area, initial operating mode, and subsystem-specific details.
- **Orbit**: Specifies initial orbital elements, epoch, and orbit type (e.g., SSO).
- **Simulation**: Includes general simulation settings, such as duration, timestep, and the atmospheric model to use.
- **Ground station**: Lists ground stations with their respective locations and elevation angles.
- **Mission design**: Gathers parameters to generate a report (e.g., data to save, figures to generate) and additional user input, such as commands to initiate SAFE mode at specific times.

Templates filled with default values for the configuration files can be found in the folder *data/* and subfolders *spacecraft/*, *orbit/*, *simulation/*, *ground_station/* and *mission_design*. It is important to place each of the 5 config files listed above in their respective folder.

If for any reason want to change where to place them, change the path constants in "main.py".


### Run a Simulation

In order to run a simulation, follow these steps:

1. Place the 5 configuration files in their respective folder
2. Activate the conda environment
3. At the root directory run:

```bash
python3 -W"ignore" src/main.py simulation.json CHESS_1.json cubeSat.json lausanne.json mission_design.json

```

where the example was provided with the files with default values.

Simulation 




- by convention, values in data folder are precised in kg, km, sec, deg.
	- except for the simulation user args where you can define your unit


- also, save state -> explain

- watt hour for electrical energy

- exceptions
	- m^2 for the cell area
	- Wh and not Ws for battery energy (since it is given in WH in power budget????)

Notes: most of the functions in the notebook are then copied and added to the framework. Notebooks can have 2 different purposes
- Analysis of some results of the simulation.
- Space to create functions that are then added to the main framework.

- uplink to ground station



=> commands displaying every 1000 timestep


## Test Case



## Documentation

Detailed information about each function and class can be found in the **documentation**. To access it, open the file "index.html" in the folder *docs/html/*.

If documentation needs to be generated again, follow these commands in the *build/* folder:

   ```bash
   cmake .. -DDOCUMENTATION=0N
   make
   ```

In addition, the default parameters for the config files are in the excel file


## Limitations and Future Work




The main design choices, as well as the assumptions and limitations of this application are detailed in the report, which can be found in the `docs/` folder.

## Authors and Acknowledgement

Authors:
- (Fall 2024 Semester) Mathilde Simoni, MSc in Computational Science and Engineering at EPFL
- To be continued!

I thank my supervisors Andrew Price and Mathieu Udriot for their help and time during the semester. I also thank Jean-Paul Kneib for being my referent professor.

