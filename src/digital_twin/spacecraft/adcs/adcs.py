"""File for the ADCS subsystem.
"""

from typing import Dict

from astropy import units as u
from astropy.time import TimeDelta
from astropy.units import Quantity
import numpy as np

from digital_twin.spacecraft import SubSystem
from digital_twin.constants import attitude_mode_dict, attitude_dict


class Adcs(SubSystem):
    def __init__(self, params: Dict, init_operating_mode: int) -> None:
        print("Initializing ADCS subsystem... ")
        self.name = "ADCS"

        super(Adcs, self).__init__()

        self.consumption_mean = {
            int(k): v * u.W for k, v in params["consumption"].items()
        }

        self.safe_flag = False

        self.attitude = attitude_mode_dict[init_operating_mode]
        # print("Initial attitude: ", attitude_dict[self.attitude])

    def update(
        self,
        old_mode: str,
        new_mode: str,
        rv: np.ndarray,
        com_window: bool,
        eclipse_status: bool,
        delta_t: TimeDelta,
    ) -> None:
        # Update the attitude
        self.attitude = attitude_mode_dict[new_mode]

        # SAFE FLAG HANDLING
        # Check safe flag triggers (cannot generate a safe flag if already in safe mode)
        if new_mode != 1 and self.safe_flag == False:
            pass  # not implemented yet for this subsystem
        # Check safe flag resolution
        if self.safe_flag == True:
            pass  # not implemented yet for this subsystem

    def get_cross_section(self, old_cross_section: Quantity) -> Quantity:
        return old_cross_section  # TODO: update when dynamic cross section feature is implemented

    def compute_power_consumed(self, mode: int) -> Quantity:
        return self.consumption_mean[mode]

    def __str__(self) -> str:
        return f"ADCS:"

    def raise_safe_flag(self) -> bool:
        return self.safe_flag

    def get_name(self) -> str:
        return self.name

    def get_attitude(self) -> int:
        return self.attitude
