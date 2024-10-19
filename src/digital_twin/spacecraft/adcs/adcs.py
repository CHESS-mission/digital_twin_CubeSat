"""Main file for the ADCS subsystem
"""

from typing import Dict

import numpy as np

from astropy import units as u
from astropy.units import Quantity

from digital_twin.spacecraft import SubSystem


class Adcs(SubSystem):
    def __init__(self, params: Dict, init_operating_mode: int) -> None:
        print("Initializing ADCS subsystem... ")

        super(Adcs, self).__init__(init_operating_mode)

    def update(
        self,
        old_mode: str,
        new_mode: str,
        rv: np.array,
        com_window: bool,
        eclipse_status: bool,
    ) -> None:
        pass

    def get_cross_section(self) -> Quantity:
        return (
            0.07487e-6 * u.km**2
        )  # right now, it is the initial value cubeSat, will need to update this

    def compute_power_consumed(self) -> float:
        return 0.0

    def __str__(self) -> str:
        return f"ADCS:"
