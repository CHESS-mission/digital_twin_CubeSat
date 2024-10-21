"""Main file for the ADCS subsystem
"""

from typing import Dict

import numpy as np

from astropy import units as u
from astropy.units import Quantity
from astropy.time import TimeDelta

from digital_twin.spacecraft import SubSystem


class Adcs(SubSystem):
    def __init__(self, params: Dict, init_operating_mode: int) -> None:
        print("Initializing ADCS subsystem... ")

        super(Adcs, self).__init__()

        self.consumption_mean = {
            int(k): v * u.W for k, v in params["consumption"].items()
        }

    def update(
        self,
        old_mode: str,
        new_mode: str,
        rv: np.array,
        com_window: bool,
        eclipse_status: bool,
        delta_t: TimeDelta,
    ) -> None:
        pass

    def get_cross_section(self, old_cross_section: Quantity) -> Quantity:
        return old_cross_section  # right now, it is the initial value cubeSat, will need to update this

    def compute_power_consumed(self, mode: int) -> Quantity:
        return self.consumption_mean[mode]

    def __str__(self) -> str:
        return f"ADCS:"
