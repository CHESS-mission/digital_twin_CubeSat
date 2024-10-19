"""Main file for the Telecom subsystem
"""

from typing import Dict

import numpy as np

from astropy import units as u

from digital_twin.spacecraft import SubSystem


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

        self.measurement_duration = 0.0
        self.communication_duration = 0.0
        self.data_storage = 0.0

    def data_storage_full(self) -> bool:
        return False

    def handshake(self) -> bool:
        return True

    def downlink_complete(self) -> bool:
        return True

    # communication finished with ground station
    def com_finished(self) -> bool:
        return False

    # scientific data transfer to ground station session
    def campaign_finished(self) -> bool:
        return False

    def update(
        self,
        old_mode: str,
        new_mode: str,
        rv: np.array,
        com_window: bool,
        eclipse_status: bool,
    ) -> None:
        pass

    def compute_power_consumed(self, mode: int) -> float:
        return self.consumption_mean[mode]

    def __str__(self) -> str:
        return f"Telecom:"
