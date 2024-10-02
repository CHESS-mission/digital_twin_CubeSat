"""Main file that runs the simulation
"""

from typing import Dict

import numpy as np

from astropy.time import Time, TimeDelta
from digital_twin.orbit_propagator import OrbitPropagator
from digital_twin.spacecraft import Spacecraft
from digital_twin.utils import get_astropy_unit


class Simulation:
    def __init__(
        self, simulation_params: Dict, orbit_params: Dict, spacecraft_params: Dict
    ) -> None:

        # TIME INITIALIZATION
        self.sim_unit = get_astropy_unit(simulation_params["units_delta_t"])
        delta_t_val = simulation_params["delta_t"] * self.sim_unit
        self.delta_t = TimeDelta(delta_t_val)
        self.init_time_local = 0 * self.sim_unit

        self.duration_sim = simulation_params["duration_sim"] * get_astropy_unit(
            simulation_params["units_duration_sim"]
        )
        self.n_timesteps = int(
            self.duration_sim.to_value(self.sim_unit) / delta_t_val.value
        )
        self.end_time_local = (self.n_timesteps * delta_t_val.value) * self.sim_unit

        self.epoch = Time(orbit_params["epoch"], format="iso", scale="utc")
        self.init_time = self.init_time_local + self.epoch
        self.end_time = self.end_time_local + self.epoch

        print("Init time (local): ", self.init_time_local)
        print("End time (local): ", self.end_time_local)
        print("Init time (global): ", self.init_time.to_value("iso"))
        print("End time (global): ", self.end_time.to_value("iso"))

        # SPACECRAFT INITIALIZATION
        self.spacecraft = Spacecraft(spacecraft_params)

        spacecraft_args = {
            "C_D": self.spacecraft.C_D,
            "A_over_m": self.spacecraft.A_over_m,
        }
        self.propagator = OrbitPropagator(orbit_params, self.epoch, spacecraft_args)

        self.tofs = TimeDelta(
            np.linspace(
                self.init_time_local, self.end_time_local, num=self.n_timesteps + 1
            )
        )  # gives the results in days

    def run(self) -> None:
        print("Simulation running...")
        eph = np.zeros((self.n_timesteps + 1, 6))
        eph[0, :3] = self.propagator.r
        eph[0, 3:] = self.propagator.v

        for t in range(0, self.n_timesteps):
            rv = self.propagator.propagate(self.delta_t)
            eph[t + 1, :3] = rv[:3]
            eph[t + 1, 3:] = rv[3:]

        # TODO: need to extract the data and plot afterwards
