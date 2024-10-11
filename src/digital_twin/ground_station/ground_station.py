"""File to store ground station related functions and objects
"""

from typing import Dict

import numpy as np

from astropy import units as u

from digital_twin.utils import get_astropy_units_angle


class GroundStation:
    def __init__(self, params: Dict) -> None:
        print("Initializing the ground station")

        self.name = params["name"]
        self.coords = np.array([params["latitude"], params["longitude"]])
        self.elev_angle = float(params["elevation_angle"]) * get_astropy_units_angle(
            params["elevation_angle_unit"]
        )

    def __str__(self):
        return f'ground station "{self.name}" located at {self.coords} with an elevation angle of {self.elev_angle}'
