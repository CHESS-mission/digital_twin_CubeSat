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
        self.HK_rate = float(params["housekeeping_data_rate"]) * (u.Mbit / u.s)
        self.max_storage = float(params["max_storage"]) * u.Mbit
        measure_unit = get_astropy_unit_time(params["max_duration_unit"])
        self.measurement_max_duration = (params["max_duration"] * measure_unit).to(u.s)
        self.x_band_rate = float(params["x_band_rate"]) * (u.Mbit / u.s)
        self.uhf_rate = float(params["uhf_rate"]) * (u.Mbit / u.s)
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

        self.measurement_duration = 0.0 * u.s

        self.is_measuring = False
        self.nb_measurement_windows = 0

        # data_storage gathers all 3 types of data (TOF, GNSS, HK = Housekeeping)
        self.data_storage = 0.0 * u.Mbit
        self.data_TOF_GNSS = 0.0 * u.Mbit
        self.data_HK = 0.0 * u.Mbit

        self.data_to_downlink = 0.0 * u.Mbit

    def can_start_measuring(self, time_elapsed: Quantity) -> bool:
        one_day = 86400 * u.s
        if (
            math.floor(self.nb_measurement_windows / self.nb_measurement_per_day)
            < (time_elapsed / one_day).value
        ):
            return True
        else:
            return False

    def data_storage_full(self) -> bool:
        if self.data_storage >= self.max_storage:
            return True
        return False

    def data_storage_empty(self, type="all") -> bool:
        if type == "all":
            if self.data_storage <= 0:
                return True
            else:
                return False
        elif type == "x_band":
            if self.data_TOF_GNSS <= 0:
                return True
            else:
                False
        elif type == "uhf":  # housekeeping
            if self.data_HK <= 0:
                return True
            else:
                return False
        else:
            raise NotImplementedError("This type of data does not exist!")

    def x_band_data_empty(self) -> bool:
        if self.data_to_downlink <= 0:
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
        pass
        if new_mode == 4:  # X_BAND
            if new_mode != old_mode:  # just switched to x_band_comm mode
                # update how much data needs to be downlinked
                self.data_to_downlink = (
                    self.data_TOF_GNSS.value * u.Mbit
                )  # not direct = or it does not do a deepcopy

            # cannot remove more than what it has
            if self.data_TOF_GNSS < (self.x_band_rate * delta_t):
                self.data_storage -= self.data_TOF_GNSS
                self.data_to_downlink -= self.data_TOF_GNSS
                self.data_TOF_GNSS = 0.0 * u.Mbit

            else:
                self.data_storage -= self.x_band_rate * delta_t
                self.data_TOF_GNSS -= self.x_band_rate * delta_t
                self.data_to_downlink -= self.x_band_rate * delta_t

            # if self.data_storage.value < 0:
            #     self.data_storage = 0.0 * u.Mbit
            # if self.data_TOF_GNSS.value < 0:
            #     self.data_TOF_GNSS = 0.0 * u.Mbit
            self.measurement_duration = 0.0 * u.s
            self.is_measuring = False

        elif new_mode == 5:  # MEASUREMENT
            if old_mode != new_mode:  # just switched to measurement mode
                self.measurement_duration = 0.0 * u.s
                self.nb_measurement_windows += 1
            # only measures during 2 hours after the 25 min of conditioning and before the last 31 min prep for dasta transfer (look at power budget)
            if (
                self.measurement_duration >= self.start_measurement
                and self.measurement_duration < self.stop_measurement
            ):
                self.data_TOF_GNSS += self.measurement_TOF_rate * delta_t
                self.data_storage += self.measurement_TOF_rate * delta_t

            self.measurement_duration += delta_t
            self.is_measuring = True

        elif new_mode == 3:  # UHF

            if self.data_HK < (self.uhf_rate * delta_t):
                self.data_storage -= self.data_HK
                self.data_HK = 0.0 * u.Mbit
            else:
                self.data_HK -= self.uhf_rate * delta_t
                self.data_storage -= self.uhf_rate * delta_t

            # if self.data_storage.value < 0:
            #     self.data_storage = 0.0 * u.Mbit
            # if self.data_HK.value < 0:
            #     self.data_HK = 0.0 * u.Mbit
            self.measurement_duration = 0.0 * u.s
            self.is_measuring = False

        else:
            self.measurement_duration = 0.0 * u.s
            self.is_measuring = False

        if (
            new_mode != 1  # and new_mode != 4 and new_mode != 3
        ):  # if not in safe mode, GNSS continuously add data
            self.data_TOF_GNSS += self.measurement_GNSS_rate * delta_t
            self.data_storage += self.measurement_GNSS_rate * delta_t
        # Housekeeping (HK) data also added in every mode
        self.data_HK += self.HK_rate * delta_t
        self.data_storage += self.HK_rate * delta_t

    def compute_power_consumed(self, mode: int) -> Quantity:
        return self.consumption_mean_gnss[mode] + self.consumption_mean_tof[mode]

    def __str__(self) -> str:
        return f"Payload:"

    def get_data_storage(self) -> Quantity:
        return self.data_storage, self.data_TOF_GNSS, self.data_HK
