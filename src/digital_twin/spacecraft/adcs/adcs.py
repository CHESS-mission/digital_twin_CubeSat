"""File for the ADCS subsystem.
"""

from typing import Dict

from astropy import units as u
from astropy.time import TimeDelta
from astropy.units import Quantity
import numpy as np

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
        rv: np.ndarray,
        com_window: bool,
        eclipse_status: bool,
        delta_t: TimeDelta,
    ) -> None:
        pass

    def get_cross_section(self, old_cross_section: Quantity) -> Quantity:
        return old_cross_section  # TODO: update when dynamic cross section feature is implemented

    def compute_power_consumed(self, mode: int) -> Quantity:
        return self.consumption_mean[mode]

    def __str__(self) -> str:
        return f"ADCS:"
