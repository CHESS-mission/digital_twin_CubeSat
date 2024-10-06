"""Main file for the EPS subsystem
"""

from typing import Dict


class Eps:
    def __init__(self, params: Dict) -> None:
        print("Initializing EPS subsystem... ")

        self.min_battery = params["min_battery"]
        self.max_battery = params["max_battery"]
        self.measure_threshold = params["measure_threshold"]
        self.com_threshold = params["com_threshold"]
        self.xb_threshold = params["xb_threshold"]

        # initialize battery level to maximum
        self.battery_level = self.max_battery
