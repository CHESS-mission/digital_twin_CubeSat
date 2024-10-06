"""Main file for the Telecom subsystem
"""

from typing import Dict


class Telecom:
    def __init__(self, params: Dict) -> None:
        print("Initializing Telecom subsystem... ")

        self.data_storage_full = False
        self.handshake = False
        self.downlink_complete = True
        self.com_finished = True  # communication window
        self.campaign_finished = (
            True  # scientific data transfer to ground station session
        )
