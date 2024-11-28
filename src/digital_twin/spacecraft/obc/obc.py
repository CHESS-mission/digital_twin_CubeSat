"""File for the On-Board Computer subsystem.
"""

from typing import Dict

from astropy import units as u
from astropy.time import TimeDelta
from astropy.units import Quantity
import numpy as np

from digital_twin.spacecraft import SubSystem


class DataStorage:
    def __init__(self, params: Dict):
        self.max_storage = float(params["max_storage"]) * u.Mbit

        # data_storage gathers all 3 types of data (TOF, GNSS, HK = Housekeeping)
        self.data_storage = 0.0 * u.Mbit
        self.data_TOF_GNSS = 0.0 * u.Mbit
        self.data_HK = 0.0 * u.Mbit

        self.data_to_downlink = (
            0.0 * u.Mbit
        )  # Updated at the start of an X-band comm window

    def update_data_storage(
        self, data_update_TOF_GNSS: Quantity, data_update_HK: Quantity
    ) -> None:

        # cannot remove more than what it has
        if (self.data_TOF_GNSS + data_update_TOF_GNSS) < 0:
            self.data_storage -= self.data_TOF_GNSS
            self.data_to_downlink -= self.data_TOF_GNSS
            self.data_TOF_GNSS = 0.0 * u.Mbit
        else:
            self.data_storage += (
                data_update_TOF_GNSS  # data usually has a negative sign
            )
            self.data_TOF_GNSS += data_update_TOF_GNSS
            self.data_to_downlink += data_update_TOF_GNSS

        if (self.data_HK + data_update_HK) < 0:
            self.data_storage -= self.data_HK
            self.data_HK = 0.0 * u.Mbit
        else:
            self.data_storage += data_update_HK
            self.data_HK += data_update_HK

    def x_band_data_empty(self) -> bool:
        if self.data_to_downlink <= 0:
            return True
        else:
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

    def data_storage_full(self) -> bool:
        if self.data_storage >= self.max_storage:
            return True
        return False

    def get_data_storage(self) -> Quantity:
        return self.data_storage, self.data_TOF_GNSS, self.data_HK

    def update(
        self,
        old_mode: str,
        new_mode: str,
        rv: np.ndarray,
        com_window: bool,
        eclipse_status: bool,
        delta_t: TimeDelta,
    ) -> None:
        if new_mode == 4:  # X_BAND
            if new_mode != old_mode:  # just switched to x_band_comm mode
                # update how much data needs to be downlinked
                self.data_to_downlink = (
                    self.data_TOF_GNSS.value * u.Mbit
                )  # not direct = or it does not do a deepcopy


class Obc(SubSystem):
    def __init__(self, params: Dict, init_operating_mode: int) -> None:
        print("Initializing OBC subsystem... ")
        self.name = "OBC"

        super(Obc, self).__init__()

        self.consumption_mean = {
            int(k): v * u.W for k, v in params["consumption"].items()
        }

        # Obc generates housekeeping data
        self.HK_rate = float(params["housekeeping_data_rate"]) * (u.Mbit / u.s)

        self.data_storage = DataStorage(params["data_storage"])

        # Safe flag handling
        self.safe_flag_spacecraft_raised = (
            False  # if a safe flag was raised by ANY subsystem
        )
        self.safe_flag_spacecraft_subsystem = (
            None  # which subsystem raised the same flag
        )

        # Safe flag regarding own subsystem
        self.safe_flag = False

    def update(
        self,
        old_mode: str,
        new_mode: str,
        rv: np.ndarray,
        com_window: bool,
        eclipse_status: bool,
        delta_t: TimeDelta,
    ) -> None:
        self.data_storage.update(
            old_mode, new_mode, rv, com_window, eclipse_status, delta_t
        )

        # SAFE FLAG HANDLING
        # check safe flag triggers (cannot generate a safe flag if already in safe mode)
        if new_mode != 1 and self.safe_flag == False:
            pass  # not implemented yet for this subsystem
        # check safe flag resolution
        if self.safe_flag == True:
            pass  # not implemented yet for this subsystem

    def compute_power_consumed(self, mode: int) -> Quantity:
        return self.consumption_mean[mode]

    def update_data_storage(
        self,
        data_update_TOF_GNSS: Quantity,
        data_update_HK: Quantity,
        delta_t: TimeDelta,
    ) -> None:
        # Housekeeping (HK) data also added in every mode
        data_update_HK += self.HK_rate * delta_t
        self.data_storage.update_data_storage(data_update_TOF_GNSS, data_update_HK)

    def get_data(self) -> Quantity:
        return self.data_storage.get_data_storage()

    def get_data_storage(self) -> DataStorage:
        return self.data_storage

    def __str__(self) -> str:
        return f"Obc:"

    # this is the method for OBC subsystem: was there an issue in the OBC subsystem itself
    def raise_safe_flag(self) -> bool:
        return self.safe_flag

    def get_name(self) -> str:
        return self.name

    def update_safe_flag(
        self, safe_flag_raised: bool, safe_flag_subsystem: str
    ) -> None:
        if safe_flag_raised:
            self.safe_flag_spacecraft_raised = True
            self.safe_flag_spacecraft_subsystem = safe_flag_subsystem
        else:
            self.safe_flag_spacecraft_raised = False
            self.safe_flag_spacecraft_subsystem = None

    def raise_spacecraft_safe_flag(self) -> bool:
        return self.safe_flag_spacecraft_raised
