"""Main file for the Telecom subsystem
"""

from typing import Dict

import numpy as np

from digital_twin.spacecraft import SubSystem


class Telecom(SubSystem):
    def __init__(self, params: Dict, init_operating_mode: int) -> None:
        print("Initializing Telecom subsystem... ")

        super(Telecom, self).__init__(init_operating_mode)

        self.data_storage_full = False
        self.handshake = False
        self.downlink_complete = True
        self.com_finished = True  # communication window
        self.campaign_finished = (
            True  # scientific data transfer to ground station session
        )

    def update(
        self,
        old_mode: str,
        new_mode: str,
        rv: np.array,
        com_window: bool,
        eclipse_status: bool,
    ) -> None:
        pass

    def compute_power_consumed(self) -> float:
        return 0.0

    def __str__(self) -> str:
        return f"Telecom:"
