"""Main file to run."""

import sys

from digital_twin import Simulation
from digital_twin.utils import parse_data_file

# Defining data paths
SIMULATION_FOLDER = "data/simulation/"
ORBIT_FOLDER = "data/orbit/"
SPACECRAFT_FOLDER = "data/spacecraft/"
GROUND_STATION_FOLDER = "data/ground_station/"
MISSION_DESIGN_FOLDER = "data/mission_design/"

# Defining result paths
RESULTS_FOLDER = "results/"


def run(args: list[str]) -> None:
    """Main function to instantiate the Simuation class and run the simulation.

    Args:
        args (List[str]): user string arguments
    """
    if len(args) != 5:
        raise TypeError(
            "Please provide 5 data files for the simulaton, orbit, spacecraft, ground station, and mission design parameters!"
        )
    simulation_params = parse_data_file(SIMULATION_FOLDER + args[0])
    orbit_params = parse_data_file(ORBIT_FOLDER + args[1])
    spacecraft_params = parse_data_file(SPACECRAFT_FOLDER + args[2])
    ground_station_params = parse_data_file(GROUND_STATION_FOLDER + args[3])
    mission_design_params = parse_data_file(MISSION_DESIGN_FOLDER + args[4])
    simulation = Simulation(
        simulation_params,
        orbit_params,
        spacecraft_params,
        ground_station_params,
        mission_design_params,
    )
    simulation.run(results_folder=RESULTS_FOLDER)


if __name__ == "__main__":
    # Extract user arguments
    if len(sys.argv) == 2:  # User arguments provided as one string
        user_args = sys.argv[1].split(" ")
        if user_args == [""]:
            run([])
        else:
            run(user_args)
    elif len(sys.argv) > 2:  # User arguments provided as separate strings
        run(sys.argv[1:])
    else:
        run([])  # No user argument
