"""Main file to run.
"""

import sys
from typing import List, Dict

from digital_twin import Simulation
from digital_twin.utils import parse_data_file


def run(args: List[str]) -> None:
    """Main function to instantiate the Simuation class and run the simulation.

    Args:
        args (List[str]): user string arguments
    """
    if len(args) != 4:
        raise TypeError(
            "Please provide 4 data files for the simulaton, orbit, spacecraft and ground station parameters!"
        )
    simulation_params = parse_data_file(args[0])
    orbit_params = parse_data_file(args[1])
    spacecraft_params = parse_data_file(args[2])
    station_params = parse_data_file(args[3])
    simulation = Simulation(
        simulation_params, orbit_params, spacecraft_params, station_params
    )
    simulation.run()


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
