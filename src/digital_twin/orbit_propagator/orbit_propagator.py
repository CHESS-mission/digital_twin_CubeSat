"""Wrapper around poliastro orbit propagator
"""

from typing import Dict
import copy

from numba import njit as jit

import numpy as np

from astropy import units as u
from astropy.units import Quantity
from astropy.time import Time

from poliastro.twobody import Orbit
from poliastro.twobody.propagation import CowellPropagator
from poliastro.core.perturbations import (
    J2_perturbation,
    atmospheric_drag,
)
from poliastro.core.propagation import func_twobody

from digital_twin.constants import atmosphere_model, earth_R, J2
from digital_twin.orbit_propagator.constants import attractor


class OrbitPropagator:
    def __init__(
        self, orbit_params: Dict, epoch: Time, spacecraft_params: Dict
    ) -> None:
        print("Initializing the propagator...")
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

        self.C_D = spacecraft_params["C_D"]
        self.A_over_m = spacecraft_params["A_over_m"]

        self.method = CowellPropagator(f=self.f)

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

    def propagate(self, delta_t: Quantity) -> np.array:
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
            C_D=self.C_D.value,
            A_over_m=self.A_over_m.value,
            rho=self.rho.value,
        )
        du_ad = np.array([0, 0, 0, ax, ay, az])

        return du_kep + du_ad

    def calculate_com_window(self):
        pass

    def calculate_eclipse_status(self):
        pass

    def update_cross_section_mass_ratio(self, new_A_over_m: Quantity) -> None:
        self.A_over_m = new_A_over_m

    def update_rho(self, position: np.array) -> None:
        alt = (np.linalg.norm(position) - earth_R.value) * u.km
        self.rho = atmosphere_model.density(alt).to_value(u.kg / u.km**3) * (
            u.kg / u.km**3
        )
