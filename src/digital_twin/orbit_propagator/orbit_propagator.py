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
        self.init_SMA = orbit_params["orbital_elements"]["SMA"] * u.km
        self.init_INC = orbit_params["orbital_elements"]["INC"] * u.deg
        self.init_ECC = orbit_params["orbital_elements"]["ECC"] * u.one
        self.init_RAAN = orbit_params["orbital_elements"]["RAAN"] * u.deg
        self.init_AOP = orbit_params["orbital_elements"]["AOP"] * u.deg
        self.init_TA = orbit_params["orbital_elements"]["TA"] * u.deg
        self.init_epoch = epoch
        self.attractor = attractor
        self.init_alt = self.init_SMA - earth_R

        # Define orbit, time units and number of timesteps
        # TODO: more general orbit
        self.init_orbit = Orbit.heliosynchronous(
            attractor=attractor,
            a=self.init_SMA,
            ecc=self.init_ECC,
            # inc = inc,
            raan=self.init_RAAN,
            argp=self.init_AOP,
            nu=self.init_TA,
            epoch=epoch,
        )

        self.current_orbit = copy.deepcopy(self.init_orbit)

        # Constant approximation for now
        # TODO: update
        self.rho = atmosphere_model.density(self.init_alt).to_value(
            u.kg / u.km**3
        ) * (u.kg / u.km**3)

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

    def update_cross_section_mass_ratio(self, new_A_over_m: Quantity):
        self.A_over_m = new_A_over_m
