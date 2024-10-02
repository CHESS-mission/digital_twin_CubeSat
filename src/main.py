"""Main file to run the simulation.
"""

import sys
from typing import List

from digital_twin import Simulation


def run(args: List[str]):
    """Main function that runs the simulation.

    Args:
        args (List[str]): _description_
    """
    print("User arguments: ", args)
    print("Starting the simulation...")
    simulation = Simulation(args)
    simulation.run()
    print("Simulation ended!")


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
