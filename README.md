# A Digital Twin of the CHESS CubeSat

This repository hosts a Python-based simulation framework developed to support the mission design and operations of the CHESS CubeSat.

### Purposes

1. *Simulation tool for mission design*: Provides tools for modeling and simulating different subsystems of the CHESS CubeSat to inform decisions and help in the development, testing, and validation tasks.

2. *Foundation for a digital twin*: Serves as the groundwork for a digital twin capable of simulating real-time satellite operations, including anomaly detection.

This framework was initially developed as part of a semester project with the EPFL Spacecraft Team in collaboration with the EPFL Space Center. The semester project report under the *docs/* folder is a valuable ressource in addition to this README file.

## Table of Contents
- [Libraries and Requirements](#libraries-and-requirements)
- [Repository Structure](#repository-structure)
- [Usage](#usage)
- [Documentation](#documentation)
- [TODOs](#todos)
- [Authors](#authors)

## Libraries and Requirements

This application requires **Python 3.9** or later. All necessary dependencies are specified in the *digital_twin_env.yaml* file. Follow these steps to create and activate the required **Conda** environment:

1. Create the environment by running the following command from the root of the repository:

```bash
conda env create -f digital_twin_env.yml
```

2. Activate the newly created environment:

```bash
conda activate digital_twin_env
```

Once the environment is activated, you are ready to run simulations!

## Repository Structure

The repository is organized as follows:


```bash
├── README.md
├── data
│   ├── atmosphere_data					# Files for atmosphere models
│   ├── ground_station					# Config files for ground stations
│   ├── mission_design					# Config files for mission design parameters
│   ├── orbit						# Config files for orbit parameters
│   ├── simulation					# Config files for simulation parameters
│   └── spacecraft					# Config files for spacecraft parameters
├── digital_twin_env.yml 				# Environment configuration file for setting up dependencies
├── docs						# Documentation and diagrams
│   ├── UML_diagram.png					# UML class diagram of the framework
│   ├── doc_generation					# Folder for generating HTML documentation
│   ├── html						# HTML documentation
│   └── parameters.xlsx					# Default parameters descriptions and sources
├── results						# Output directory for simulation results
│   ├── data						# Numpy output files 
│   └── figures						# Generated figures
└── src							# Source code
    ├── digital_twin
    │   ├── constants.py				# General constants for the simulation
    │   ├── ground_station				# Ground station module
    │   │   └── ground_station.py
    │   ├── mode_switch.py				#  Mode switch decision tree implementation
    │   ├── orbit_propagator				# Orbit propagator module
    │   │   ├── atmosphere_model.py			#  Atmosphere models implementation
    │   │   ├── constants.py				# Constants for propagation
    │   │   └── orbit_propagator.py			# Orbit propagator implementation
    │   ├── plotting.py					# Functions to plot results
    │   ├── report.py					# Functions to generate outputs
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
    │   │   ├── spacecraft.py				# Spacecraft class implementation
    │   │   ├── subsystem.py				# Subsystem interface
    │   │   └── telecom					# Telecom module
    │   │       └── telecom.py
    │   └── utils.py					# Utility functions
    ├── main.py						# Main entry point for running the simulation
    └── notebooks					# Notebooks for visualizations, experiments, and analysis
```

### Notes


- *Notebooks/*: Each notebook includes a description at the beginning. They are used for advanced plotting, specific analyses, experiments, and testing functions later integrated into the framework.

- *data/atmosphere_data/*: Contains data for **NRLMSISE00**, **JB2008** and **solar_activity** atmosphere models. Data for the first two are updated automatically by the *ATMOS* library during simulations. Solar activity data ($F10.7$ index and $Ap$) for the third model are pre-extracted from DRAMA software and do not require updating.

- "\_\_init\_\_.py": These files exist in all module directories but are omitted for clarity in the tree.

- "mode_switch.py": the mode switch decision tree used for the `ModeSwitch` class implementation can be found in the Appendix of the report under the *docs/* folder.

- The full code organization and design decisions are explained in detail in the report available in *docs/*.

## Usage

### Input
User parameters are divided into five categories, each with its own JSON configuration file:

- **Spacecraft**: Defines general parameters like mass, cross-sectional area, initial operating mode, and subsystem-specific details.
- **Orbit**: Specifies initial orbital elements, epoch, and orbit type (e.g., SSO).
- **Simulation**: Includes general simulation settings, such as duration, timestep, and the atmospheric model to use.
- **Ground station**: Lists ground stations with their respective locations and elevation angles.
- **Mission design**: Gathers parameters to generate a report (e.g., data to save, figures to generate) and additional user input, such as commands to initiate SAFE mode at specific times.

Default configuration templates for each category are located in the *data/* folder and its subdirectories (*spacecraft/*, *orbit/*, *simulation/*, *ground_station/*, and *mission_design/*). The file "parameters.xlsx" under *docs/* lists default parameter descriptions and sources for the CHESS CubeSat (for other satellite models, this file should be updated accordingly). Place the relevant configuration files in their corresponding subfolder before running a simulation.

**Units**

Unless otherwise specified, all values in the configuration files under *data/* follow these standard units:

- Mass: **kilograms (kg)**
- Time: **seconds (s)**
- Distance: **kilometers (km)**
- Angles: **degrees (deg)**
- Power: **Watts (W)**
- Data size: **Megabit (Mbit)**

Exceptions:

1. Custom Units:
	
	Some parameters allow unit customization. In these cases, specify the unit explicitly within the configuration file.
	
	- Time: Choose from `"second"`, `"hour"`, `"day"`, or `"year"`
	- Angle: Choose from `"degree"` or `"radian"`
	
	Example fields: `"units_duration_sim"` or `"unit_delta_t"` in the **simulation** config file, or `"elevation_angle_unit"` in the **ground station** config file.

2. Battery Energy: 

	Expressed in Watt-Hours (Wh) rather than Watt-seconds, following the EPFL Spacecraft Team's power budget conventions. This applies to initial battery energy, minimum and maximum thresholds, and energy-related parameters like UHF and X-band communication levels.

3. Solar Cell Area: 

	Measured in square meters (m²) for practicality due to small values.


### Starting a Simulation

In order to run a simulation, follow these steps:

1. Place the 5 configuration files in their respective folder
2. Activate the conda environment
3. Make sure you are connected to the internet (some files are automatically updated by libraries)
3. At the root directory, run:

	```bash
	python3 -W"ignore" src/main.py simulation_template.json orbit_template.json spacecraft_template.json ground_station_template.json mission_design_template.json

	```

This example provided uses the files with default values. It is important to keep the file arguments in the specified order.

### Verbose 

During the simulation, relevant information and progress updates can be displayed in the terminal:

1. **Earth Orientation Data Update**: If running the application for the first time or after a prolonged period, the necessary Earth orientation data tables required by the `astropy.utils.iers`z
 module are automatically downloaded and a message is displayed. If the tables are already up-to-date, a message confirms their current status.

2. **Simulation Logging (Optional)**: If `"verbose"` is set to `"yes"` in the **simulation** configuration file, detailed logs are displayed, including class initializations and progress updates at every 1,000 timesteps, along with the final runtime.

3. **Parameter Display (Optional)**: If `"print_parameters"` is set to `"yes"` in the **simulation** configuration file, the initial simulation parameters are printed before entering the main simulation loop.


### Output

Upon completing the simulation, all outputs are stored in the *results/* directory. The generated outputs depend on user choices specified in the **mission design** configuration file, where each field must be set to `"yes"` or `"no"` to enable or disable the respective output.

The outputs are divided into two categories:

- Data: Stores simulation results as numpy arrays (for each timestep), which can be used for further analysis or visualizations (refer to the notebooks for examples). Additionally, two JSON files are created to save the final states of the orbit and spacecraft, allowing them to be reused as starting configurations for subsequent simulations. The other three configuration files (**simulation**, **ground station**, **mission design**) are not directly linked to simulation states and are therefore not automatically saved.

- Figures: Generates visual plots of various parameters.

**Detailed Output Description**

- **data**:
	- `"telecom_data"`: Visibility windows (1 = visible, 0 = not visible), data storage (total, housekeeping, scientific data)
	- `"eps_data"`: Eclipse status (1 = eclipse, 0 = not in eclipse), battery energy level, power consumption, and generation
	- `"modes"`: Operating modes (refer to "constants.py" for mode definitions)
	- `"altitude_data"`: Satellite altitude
	- `"orbital_element_data"`: Orbital elements including altitude, RAAN, AOP, ECC, and INC
	- `"eclipse_data"`: Eclipse status (1 = eclipse, 0 = not in eclipse)
	- `"orbit_state"`: Final orbit state as a JSON file
	- `"spacecraft_state"`: Final spacecraft state as a JSON file
	- `"density"`: Air density at the satellite's position.

All data arrays are saved with an accompanying `"times.npy"` array for use in plotting.

- **figures**:
	- `"orbital_elem_evolution"`: Evolution of RAAN, AOP, ECC, INC, and altitude over time
	- `"trajectory_2d"`: 2D plot of the initial orbit
	- `"trajectory_3d"`: 3D plot showing the initial orbit and full trajectory
	- `"groundtrack"`: Ground track of the satellite
	- `"modes"`:  Horizontal bar plot of spacecraft operating modes during the simulation
	- `"dashboard"`:  Combined plot showing operating modes, visibility windows, and eclipse status
	- `"battery_energy"`:  Battery energy level over time
	- `"power_consumption"`: Power consumption over time
	- `"power_generation"`:  Power generation over time
	- `"power_balance"`:  Net power (generation minus consumption) over time
	- `"data_storage"`: Storage usage (total, scientific, and housekeeping data)
	- `"visibility_windows"`: Boolean bar plot of visibility windows
	- `"eclipse_windows"`: Boolean bar plot of eclipse status

## Documentation

Comprehensive information on all functions, classes, and modules is available in the generated documentation. To explore it, open the "index.html" file located in the *docs/html/* folder. This provides a user-friendly interface for navigating the code's structure and available features.

The documentation is created using [*Sphinx*](https://www.sphinx-doc.org/en/master/) mainly following a [step-by-step tutorial](https://www.youtube.com/watch?v=BWIrhgCAae0).  If changes are made to the codebase structure (such as adding or removing modules or packages), the documentation tree must be updated. Follow these steps to regenerate it:

1. Navigate to the *docs/doc_generation/* directory
2. Remove all **.rst** files EXCEPT "index.rst" 
3. Return to the projec's root repository
4. Activate the virtual environment (**digital_twin_env**)
5. Run the following command to rebuild the documentation structure:

	```bash
	sphinx-apidoc -o docs/doc_generation src/ --force
	```
5. Go back to the *docs/doc_generation/* directory
6. Generate the HMTL documentation
	```bash
	make html
	```

The updated documentation will be located in the *docs/html/* folder.

If only code changes are made (such as modifying functions or classes), and no new files are added or removed, simply run `make html` in the *docs/doc_generation/* directory to re-generate the documentation.

### Type hints

In this project, each function is defined with type hints to indicate the expected types of variables. While Python is not a statically typed language and does not enforce these types at runtime, type hints serve as helpful guidance for developers, providing insight into what type of data a function expects or returns. Additionally, tools like **Mypy** can be used for static type checking, allowing developers to detect type errors before runtime.

Commonly used types include basic ones like `float`, `int`, `str`, and `bool`, as well as more complex types such as `list` (e.g., `list[int]`), `dict`, `set`, and `tuple`, or *Numpy* arrays (`np.ndarray`). When a variable can be either `None` or another type, its type hint is written as `Optional[other_type]`.

In this application, variables are frequently passed as `Quantity`. These are special objects from the **Astropy** library that allows to associate a number with a physical unit from the `astropy.units` module.

The type hints for Quantities can be specified in two ways:

- Specific Unit: When the unit is known, for example:
	- `Quantity[u.km]` for a distance in kilometers.

- Physical Type: When only the physical type is known (e.g., "length", "time", "angle", "bandwidth"), for example:
	- `Quantity["length"]`

For a comprehensive list of physical types available in Astropy, refer to this [link](https://docs.astropy.org/en/stable/units/ref_api.html#module-astropy.units.physical).

For more details about type annotations with units, see the [Astropy documentation](https://docs.astropy.org/en/stable/units/type_hints.html#).

## TODOs

The next potential steps for this project are outlined in **Section 6** of the project report. You can find the report in the *docs/* folder.

## Authors

### Initial development (fall 2024)

- **Author**: Mathilde Simoni, MSc in Computational Science and Engineering at EPFL (mathilde.simoni@epfl.ch)
- **Supervisors**:
	- Andrew Price (andrew.prince@epfl.ch)
	- Mathieu Udriot (mathieu.udriot@epfl.ch)
- **Referent Professor**:
	- Jean-Paul Kneib

Future authors are encouraged to add their names and details as they contribute to this  project.
