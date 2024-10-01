"""Main file that runs the simulation
"""

from digital_twin.orbit_propagator import OrbitPropagator
from digital_twin.spacecraft import Spacecraft

class Simulation:
    def __init__(self, user_args: str) -> None:
        self.propagator = OrbitPropagator()
        self.spacecraft = Spacecraft()
        
    def run(self) -> None:
        print("Simulation running...")