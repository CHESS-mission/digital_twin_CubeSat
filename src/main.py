"""Main file to run the simulation.
"""

import json
import sys
from typing import List, Dict

import numpy as np

from digital_twin import Simulation


def run(args: List[str]) -> None:
    """Main function that runs the simulation.

    Args:
        args (List[str]): user string arguments
    """
    if len(args) != 3:
        raise TypeError(
            "Please provide 3 data files for the simulaton, orbit and spacecraft parameters!"
        )
    simulation_params = parse_data_file(args[0])
    orbit_params = parse_data_file(args[1])
    spacecraft_params = parse_data_file(args[2])
    print("Starting the simulation...")
    simulation = Simulation(simulation_params, orbit_params, spacecraft_params)
    simulation.run()
    print("Simulation ended!")


def parse_data_file(file_path: str) -> Dict:
    """Parse JSON file provided by user

    Args:
        file_path (str): path to JSON file

    Returns:
        Dict: data from file in a dictionary
    """
    f = open(file_path)
    data = json.load(f)
    f.close()
    return data


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
