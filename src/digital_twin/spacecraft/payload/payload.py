"""File for the payload subsystem (gnss and tof together).
"""

import math
from typing import Dict

from astropy import units as u
from astropy.time import TimeDelta
from astropy.units import Quantity
import numpy as np

from digital_twin.spacecraft import SubSystem
from digital_twin.utils import get_astropy_unit_time


class Payload(SubSystem):
    def __init__(self, params: Dict, init_operating_mode: int) -> None:
        print("Initializing Payload subsystem... ")
        self.name = "Payload"

        super(Payload, self).__init__()

        self.consumption_mean_gnss = {
            int(k): v * u.W for k, v in params["consumption_gnss"].items()
        }
        self.consumption_mean_tof = {
            int(k): v * u.W for k, v in params["consumption_tof"].items()
        }

        self.measurement_TOF_rate = float(params["measurement_TOF_rate"]) * (
            u.Mbit / u.s
        )
        self.measurement_GNSS_rate = float(params["measurement_GNSS_rate"]) * (
            u.Mbit / u.s
        )
        measure_unit = get_astropy_unit_time(params["max_duration_unit"])
        self.measurement_max_duration = (params["max_duration"] * measure_unit).to(u.s)
        self.nb_measurement_per_day = int(params["nb_measurements_per_day"])

        self.start_measurement = (
            float(params["measurement_pre_conditioning_time"]) * u.s
        )
        measurement_post_conditioning = (
            float(params["measurement_post_conditioning_time"]) * u.s
        )
        self.stop_measurement = (
            self.measurement_max_duration - measurement_post_conditioning
        )

        self.measurement_duration = float(params["init_measurement_duration"]) * u.s
        self.is_measuring = False if params["init_is_measuring"] == "false" else True
        self.nb_measurement_windows = 0

        self.safe_flag = False if params["init_safe_flag"] == "false" else True

    def can_start_measuring(self, time_elapsed: Quantity) -> bool:
        one_day = 86400 * u.s
        if (
            math.floor(self.nb_measurement_windows / self.nb_measurement_per_day)
            < (time_elapsed / one_day).value
        ):
            return True
        else:
            return False

    def campaign_finished(self) -> bool:
        if not self.is_measuring:
            return True
        else:
            if self.measurement_duration >= self.measurement_max_duration:
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
        if new_mode == 5:  # MEASUREMENT
            if old_mode != new_mode:  # just switched to measurement mode
                self.measurement_duration = 0.0 * u.s
                self.nb_measurement_windows += 1
            self.measurement_duration += delta_t
            self.is_measuring = True

        else:
            self.measurement_duration = 0.0 * u.s
            self.is_measuring = False

        # SAFE FLAG HANDLING
        # check safe flag triggers (cannot generate a safe flag if already in safe mode)
        if new_mode != 1 and self.safe_flag == False:
            pass  # not implemented yet for this subsystem
        # check safe flag resolution
        if self.safe_flag == True:
            pass  # not implemented yet for this subsystem

    def compute_power_consumed(self, mode: int) -> Quantity:
        return self.consumption_mean_gnss[mode] + self.consumption_mean_tof[mode]

    def __str__(self) -> str:
        return f"Payload:"

    def compute_data_update(self, new_mode: int, delta_t: TimeDelta) -> Quantity:
        # only measures during 2 hours after the 25 min of conditioning and before the last 31 min prep for dasta transfer (look at power budget)
        data = 0.0 * u.Mbit
        if (
            new_mode == 5
            and self.measurement_duration >= self.start_measurement
            and self.measurement_duration < self.stop_measurement
        ):
            data += self.measurement_TOF_rate * delta_t

        if new_mode != 1:  # if not in safe mode, GNSS continuously add data
            data += self.measurement_GNSS_rate * delta_t

        return data, 0.0 * u.Mbit  # second term for HK data which is 0 in this case

    def raise_safe_flag(self) -> bool:
        return self.safe_flag

    def get_name(self) -> str:
        return self.name

    def get_measurement_variables(self) -> tuple:
        return self.measurement_duration, self.is_measuring
