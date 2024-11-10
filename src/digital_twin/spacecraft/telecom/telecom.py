"""File for the Telecom subsystem.
"""

from typing import Dict

from astropy import units as u
from astropy.time import TimeDelta
from astropy.units import Quantity
import numpy as np

from digital_twin.spacecraft import SubSystem
from digital_twin.spacecraft.payload import Payload
from digital_twin.utils import get_astropy_unit_time


class Telecom(SubSystem):
    def __init__(self, params: Dict, init_operating_mode: int) -> None:
        print("Initializing Telecom subsystem... ")

        super(Telecom, self).__init__()

        self.consumption_mean_uhf = {
            int(k): v * u.W for k, v in params["consumption_uhf"].items()
        }
        self.consumption_mean_x_band = {
            int(k): v * u.W for k, v in params["consumption_x_band"].items()
        }
        com_unit = get_astropy_unit_time(params["communication_max_duration_unit"])
        self.com_max_duration = (params["communication_max_duration"] * com_unit).to(
            u.s
        )

        # variables to track if communication is occuring and how much time (UHF_COM)
        self.com_duration = 0.0 * u.s
        self.is_communicating = False

    def handshake(self) -> bool:
        return True

    def downlink_xband_complete(self, payload: Payload) -> bool:
        # stop downlink when all data transferred
        if payload.data_storage_empty(type="x_band"):
            return True
        return False

    def com_finished(self) -> bool:
        if not self.is_communicating:
            return True
        # stop the communication when maximum time reached
        else:
            if self.com_duration >= self.com_max_duration:
                return True
            return False

    def update(
        self,
        old_mode: str,
        new_mode: str,
        rv: np.ndarray,
        com_window: bool,
        eclipse_status: bool,
        delta_t: TimeDelta,
    ) -> None:
        if new_mode == 3:  # UHF_COM
            self.is_communicating = True
            if old_mode != new_mode:  # just switched to communication mode
                self.com_duration = 0.0 * u.s
            self.com_duration += delta_t  # increase communication time
        else:
            self.is_communicating = False
            self.com_duration = 0.0 * u.s

    def compute_power_consumed(self, mode: int) -> Quantity:
        return self.consumption_mean_uhf[mode] + self.consumption_mean_x_band[mode]

    def __str__(self) -> str:
        return f"Telecom:"
