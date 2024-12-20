"""File to store ground station related functions and objects
"""

from typing import Dict, Tuple, Any

from astropy import units as u
from astropy.coordinates import CartesianRepresentation, SphericalRepresentation
from astropy.units import Quantity
import numpy as np

from digital_twin.constants import earth_R
from digital_twin.utils import get_astropy_units_angle


class GroundStation:
    def __init__(self, params: Dict) -> None:
        print("Initializing the ground station")
        self.name = params["name"]
        self.color = params["color"]  # for groundtrack plotting

        # Get coordinates
        self.station_coords = [params["latitude"], params["longitude"]] * u.deg
        spherical_coords = SphericalRepresentation(
            lat=self.station_coords[0],
            lon=self.station_coords[1],
            distance=earth_R,  # Earth's average radius
        )
        self.cartesian_coords = spherical_coords.represent_as(CartesianRepresentation)

        self.elev_angle = float(params["elevation_angle"]) * get_astropy_units_angle(
            params["elevation_angle_unit"]
        )

    @property
    def cartesian_coords(self) -> CartesianRepresentation:
        return self._cartesian_coords

    @cartesian_coords.setter
    def cartesian_coords(self, coords: Any) -> None:
        self._cartesian_coords = coords

    @property
    def elev_angle(self) -> Quantity["angle"]:
        return self._elev_angle

    @elev_angle.setter
    def elev_angle(self, angle: Quantity["angle"]) -> None:
        self._elev_angle = angle

    def __str__(self) -> str:
        return f'ground station "{self.name}" located at {self.station_coords} with an elevation angle of {self.elev_angle}'

    def get_name_pos_color(self) -> Tuple[str, np.ndarray]:
        return self.name, self.station_coords, self.color
