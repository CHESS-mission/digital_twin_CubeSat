"""Main file for the spacecraft object
"""

from digital_twin.spacecraft.eps import Eps
from digital_twin.spacecraft.telecom import Telecom
from digital_twin.spacecraft.adcs import Adcs

class Spacecraft:
    def __init__(self) -> None:
        print("Initializing the spacecraft...")
        self.eps_subsystem = Eps()
        self.telecom_subsystem = Telecom()
        self.adcs_subsystem = Adcs()