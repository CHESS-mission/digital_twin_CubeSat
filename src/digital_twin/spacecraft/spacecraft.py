"""Main file for the spacecraft object
"""

from typing import Dict

from astropy import units as u
from astropy.units import Quantity

from digital_twin.spacecraft.eps import Eps
from digital_twin.spacecraft.telecom import Telecom
from digital_twin.spacecraft.adcs import Adcs


class Spacecraft:
    def __init__(self, params: Dict) -> None:
        print("Initializing the spacecraft...")
        self.C_D = params["general"]["C_D"] * u.one
        self.cross_section = params["general"]["cross_section_area"] * u.km**2
        self.mass = params["general"]["mass"] * u.kg
        self.A_over_m = ((self.cross_section) / (self.mass)).to_value(
            u.km**2 / u.kg
        ) * (u.km**2 / u.kg)
        self.B = self.C_D * self.A_over_m  # Approximation of the ballistic coefficient

        # SUBSYSTEMS INITIALIZATION
        self.eps_subsystem = Eps(params["eps"])
        self.telecom_subsystem = Telecom(params["telecom"])
        self.adcs_subsystem = Adcs(params["adcs"])

        @property
        def C_D(self) -> Quantity:
            return self.C_D

        @property
        def A_over_m(self) -> Quantity:
            return self.A_over_m

        @property
        def eps(self) -> Eps:
            return self.eps_subsystem

        @property
        def telecom(self) -> Telecom:
            return self.telecom_subsystem
