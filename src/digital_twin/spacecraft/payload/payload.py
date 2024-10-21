"""Main file for the payload subsystem (gnss and tof together)
"""

from typing import Dict

import numpy as np

from astropy import units as u
from astropy.units import Quantity
from astropy.time import TimeDelta

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

        self.measurement_rate = float(params["measurement_rate"]) * (u.Mbit / u.s)
        self.max_storage = float(params["max_storage"]) * u.Mbit
        measure_unit = get_astropy_unit_time(params["max_duration_unit"])
        self.measurement_max_duration = (params["max_duration"] * measure_unit).to(u.s)
        self.x_band_rate = float(params["x_band_rate"]) * (u.Mbit / u.s)

        self.measurement_duration = 0.0
        self.data_storage = 0.0
        self.is_measuring = False

    def data_storage_full(self) -> bool:
        if self.data_storage >= self.max_storage:
            return True
        return False

    def data_storage_empty(self) -> bool:
        if self.data_storage <= 0:
            return True
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
        rv: np.array,
        com_window: bool,
        eclipse_status: bool,
        delta_t: TimeDelta,
    ) -> None:
        pass
        if new_mode == 4:  # X_BAND
            self.data_storage -= self.x_band_rate * delta_t
            self.measurement_duration = 0.0 * u.s
            self.is_measuring = False
        if new_mode == 5:  # MEASUREMENT
            self.data_storage += self.measurement_rate * delta_t
            if old_mode != new_mode:  # just switched to measurement mode
                self.measurement_duration = 0.0 * u.s
            self.measurement_duration += delta_t
            self.is_measuring = True
        else:
            self.measurement_duration = 0.0 * u.s
            self.is_measuring = False

    def compute_power_consumed(self, mode: int) -> Quantity:
        return self.consumption_mean_gnss[mode] + self.consumption_mean_tof[mode]

    def __str__(self) -> str:
        return f"Payload:"
