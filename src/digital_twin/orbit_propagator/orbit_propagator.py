"""Wrapper around poliastro orbit propagator
"""

from typing import Dict
import copy

from numba import njit as jit

import numpy as np

from astropy import units as u
from astropy.units import Quantity
from astropy.time import Time
from astropy.coordinates import (
    CartesianRepresentation,
    SphericalRepresentation,
    GCRS,
    ITRS,
    get_body_barycentric_posvel,
)

from poliastro.twobody import Orbit
from poliastro.twobody.propagation import CowellPropagator
from poliastro.core.perturbations import (
    J2_perturbation,
    atmospheric_drag,
)
from poliastro.core.propagation import func_twobody
from poliastro.core.events import eclipse_function

from digital_twin.constants import atmosphere_model, earth_R, J2, sun_R, earth_k
from digital_twin.orbit_propagator.constants import (
    attractor,
    star_string,
    attractor_string,
)
from digital_twin.ground_station import GroundStation
from digital_twin.utils import angle_between_vectors


class OrbitPropagator:
    def __init__(self, orbit_params: Dict, epoch: Time) -> None:
        print("Initializing the propagator...")

        # INITAL CHECKS
        if star_string != "Sun":
            raise ValueError(
                "Simulation only implemented for attractor orbiting about the SUN!"
            )
        if attractor_string != "Earth":
            raise ValueError(
                "Simulation only implemented for satellite orbiting about the EARTH!"
            )

        # INITIALIZE ORBITAL ELEMENTS
        init_SMA = orbit_params["orbital_elements"]["SMA"] * u.km
        init_INC = orbit_params["orbital_elements"]["INC"] * u.deg
        init_ECC = orbit_params["orbital_elements"]["ECC"] * u.one
        init_RAAN = orbit_params["orbital_elements"]["RAAN"] * u.deg
        init_AOP = orbit_params["orbital_elements"]["AOP"] * u.deg
        init_TA = orbit_params["orbital_elements"]["TA"] * u.deg

        # TODO: more general orbit
        self.init_orbit = Orbit.heliosynchronous(
            attractor=attractor,
            a=init_SMA,
            ecc=init_ECC,
            # inc = inc,
            raan=init_RAAN,
            argp=init_AOP,
            nu=init_TA,
            epoch=epoch,
        )

        self.current_orbit = copy.deepcopy(self.init_orbit)

        # Approximation updated every 500 timesteps for now (TODO: update)
        self.rho = atmosphere_model.density(init_SMA - earth_R).to_value(
            u.kg / u.km**3
        ) * (u.kg / u.km**3)
        self.countdown_rho = 500  # rho will be updated when countdown reaches 0

        self.method = CowellPropagator(f=self.f)

        self.drag_params = {
            "C_D": 0.0 * u.one,
            "A_over_m": 0.0 * (u.km**2 / u.kg),
        }  # For drag calculations, updated every time propagate() is called for dynamic drag calculations

    def __str__(self) -> str:
        return f"Initial orbit: {self.init_orbit} with initial period {self.init_orbit.period}"

    @property
    def r(self) -> np.array:
        return self.current_orbit.r

    @property
    def v(self) -> np.array:
        return self.current_orbit.v

    def get_initial_orbit(self) -> Orbit:
        return self.init_orbit

    def propagate(
        self, delta_t: Quantity, C_D: Quantity, A_over_m: Quantity
    ) -> np.array:
        # update grad parameters
        self.drag_params["C_D"] = C_D
        self.drag_params["A_over_m"] = A_over_m

        new_orbit = self.current_orbit.propagate(delta_t, method=self.method)
        self.current_orbit = new_orbit
        rv = np.zeros(6)
        rv[:3] = self.current_orbit.r
        rv[3:] = self.current_orbit.v
        return rv

    def a_d(self, t0, state, k, J2, R, C_D, A_over_m, rho):
        return J2_perturbation(t0, state, k, J2, R) + atmospheric_drag(
            t0, state, k, C_D, A_over_m, rho
        )

    def f(self, t0, state, k):
        du_kep = func_twobody(t0, state, k)

        if self.countdown_rho > 0:
            self.countdown_rho -= 1
        else:
            self.update_rho(state[:3])  # update air density based on current altitude
            self.countdown_rho = 500  # restart countdown

        ax, ay, az = self.a_d(
            t0,
            state,
            k,
            J2=J2.value,
            R=earth_R.value,
            C_D=self.drag_params["C_D"].value,
            A_over_m=self.drag_params["A_over_m"].value,
            rho=self.rho.value,
        )
        du_ad = np.array([0, 0, 0, ax, ay, az])

        return du_kep + du_ad

    def calculate_vis_window(self, ground_stations: np.array) -> np.array:
        visibility = []
        satellite_coords = np.reshape(np.array(self.current_orbit.r), (1, 3))
        obstime = self.current_orbit.epoch

        for gs in ground_stations:
            cartesian_coords = gs.cartesian_coords
            itrs_object = ITRS(cartesian_coords, obstime=obstime)
            gcrs_object = itrs_object.transform_to(GCRS(obstime=obstime))
            gs_coords = np.reshape(np.array(gcrs_object.data.xyz.value), (1, 3))
            distance_vec = satellite_coords - gs_coords
            angle = angle_between_vectors(distance_vec, gs_coords)
            if (
                angle >= (90 * u.deg - gs.elev_angle).value
            ):  # take into account ground station elevation angle
                visibility.append(False)
            else:
                visibility.append(True)

        return np.array(visibility)

    def calculate_eclipse_status(self) -> bool:
        # Taking into account umbra and penumbra for now
        epoch = self.current_orbit.epoch
        r = self.current_orbit.r.value
        v = self.current_orbit.v.value

        # Position vector of Sun wrt Solar System Barycenter
        r_sec_ssb = get_body_barycentric_posvel("Sun", epoch)[0]
        r_pri_ssb = get_body_barycentric_posvel("Earth", epoch)[0]
        r_sec = ((r_sec_ssb - r_pri_ssb).xyz << u.km).value

        eclipse = eclipse_function(
            earth_k.value,
            np.hstack((r, v)),
            r_sec,
            sun_R.value,
            earth_R.value,
            umbra=False,
        )

        return (
            eclipse <= 0.0
        )  # if <= 0, the satellite is in eclipse and doesn't get sunlight

    def update_rho(self, position: np.array) -> None:
        alt = (np.linalg.norm(position) - earth_R.value) * u.km
        self.rho = atmosphere_model.density(alt).to_value(u.kg / u.km**3) * (
            u.kg / u.km**3
        )
