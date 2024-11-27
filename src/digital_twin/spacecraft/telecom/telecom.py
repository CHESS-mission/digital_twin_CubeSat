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
        self.name = "Telecom"

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
        self.x_band_rate = float(params["x_band_rate"]) * (u.Mbit / u.s)
        self.uhf_rate = float(params["uhf_rate"]) * (u.Mbit / u.s)

        t_com_unit = get_astropy_unit_time(
            params["max_time_without_communication_unit"]
        )
        self.t_max_no_com = (params["max_time_without_communication"] * t_com_unit).to(
            u.s
        )  # maximum time allowed without communicating to ground station
        self.t_no_com = 0.0 * u.s

        # variables to track if communication is occuring and how much time (UHF_COM)
        self.com_duration = 0.0 * u.s
        self.is_communicating = False

        self.alternating = False  # bool which prevents from alternating between x-band-comm and uhf-comm mode during visibility window

        self.safe_flag = False

    def handshake(self) -> bool:
        return True

    def can_downlink(self) -> bool:
        if (
            self.alternating
        ):  # prevents from alternating between uhf-comm and x-band-comm endlessly during visibility window
            return False
        else:
            return True

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

            if old_mode == 4:  # Just passed from x band to uhf
                # We want to prevent the satellite to go back to x-band directly bc of the little GNSS data added
                self.alternating = True

            self.is_communicating = True
            if old_mode != new_mode:  # just switched to communication mode
                self.com_duration = 0.0 * u.s
                self.t_no_com = 0.0 * u.s
            self.com_duration += delta_t  # increase communication time

        else:
            self.is_communicating = False
            if new_mode != 4:
                self.t_no_com += (
                    delta_t  # X_BAND mode (4) is considered as communication
                )
            self.com_duration = 0.0 * u.s
            self.alternating = False  # set back because we re not alternating between uhf-comm and x-band-comm anymore

        # SAFE FLAG HANDLING
        # check safe flag triggers (cannot generate a safe flag if already in safe mode)
        if new_mode != 1 and self.safe_flag == False:
            # there might be other safe flag triggers later
            if self.t_no_com > self.t_max_no_com:
                self.safe_flag = True

        # check safe flag resolution
        if self.safe_flag == True:
            if com_window:
                self.safe_flag = False

    def compute_power_consumed(self, mode: int) -> Quantity:
        return self.consumption_mean_uhf[mode] + self.consumption_mean_x_band[mode]

    def compute_data_update(self, new_mode: int, delta_t: TimeDelta) -> Quantity:
        if new_mode == 4:  # X_BAND
            return -1 * self.x_band_rate * delta_t, 0 * u.Mbit
        elif new_mode == 3:  # UHF
            return 0 * u.Mbit, -1 * self.uhf_rate * delta_t
        else:
            return 0.0 * u.Mbit, 0.0 * u.Mbit

    def __str__(self) -> str:
        return f"Telecom:"

    def raise_safe_flag(self) -> bool:
        return self.safe_flag

    def get_name(self) -> str:
        return self.name
