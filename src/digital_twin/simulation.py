"""Main file that runs the simulation
"""

from typing import Dict
import time

import numpy as np

from astropy.time import Time, TimeDelta

from poliastro.core.elements import rv2coe

from digital_twin.orbit_propagator import OrbitPropagator
from digital_twin.spacecraft import Spacecraft
from digital_twin.utils import get_astropy_unit
from digital_twin.constants import earth_R, earth_k
from digital_twin.plotting import plot_1d, seconds_to_days


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

        print("Number of timesteps:", self.n_timesteps)
        start_for_loop = time.time()

        for t in range(0, self.n_timesteps):
            if t % 500 == 0:
                print(t)
            rv = self.propagator.propagate(self.delta_t)
            eph[t + 1, :3] = rv[:3]
            eph[t + 1, 3:] = rv[3:]

        end_for_loop = time.time()
        duration_for_loop = end_for_loop - start_for_loop
        print("Time: ", duration_for_loop)

        # Extract the data
        rr, vv, SMAs, ECCs, INCs, RAANs, AOPs, TAs, altitudes = (
            self.extract_propagation_data(eph)
        )
        print("Max change in altitude:", (np.max(altitudes) - np.min(altitudes)))

        # Plot the data
        self.plot_propagation_data(self.tofs, RAANs, AOPs, ECCs, INCs, altitudes)

    def extract_propagation_data(self, eph: np.array):
        rr = eph[:, :3]
        vv = eph[:, 3:]
        orbital_params = np.array([rv2coe(earth_k, r, v) for r, v in zip(rr, vv)])
        ps = orbital_params[:, 0]
        ECCs = orbital_params[:, 1]
        INCs = orbital_params[:, 2]
        RAANs = orbital_params[:, 3]
        AOPs = orbital_params[:, 4]
        TAs = orbital_params[:, 5]
        SMAs = np.divide(
            ps, 1 - np.multiply(ECCs, ECCs)
        )  # Formula linking semi-latus rectum to semi-major axis
        altitudes = np.linalg.norm(eph, axis=1) - earth_R.value

        return rr, vv, SMAs, ECCs, INCs, RAANs, AOPs, TAs, altitudes

    def plot_propagation_data(
        self,
        tofs: np.array,
        RAANs: np.array,
        AOPs: np.array,
        ECCs: np.array,
        INCs: np.array,
        altitudes: np.array,
    ):
        ticks_angle = np.array([0, 1 / 2 * np.pi, np.pi, 3 / 2 * np.pi, 2 * np.pi])
        tick_labels_angle = np.array(
            [r"$0$", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"]
        )
        plot_1d(
            tofs.to_value("second"),
            RAANs,
            "RAAN Over Time",
            r"Time ($day$)",
            r"RAAN ($rad$)",
            step=100,
            fill_under=False,
            remove_box=True,
            y_range=(min(RAANs) - 0.2, max(RAANs) + 0.2),
            scatter=True,
            custom_y_ticks=True,
            y_ticks=ticks_angle,
            y_tick_labels=tick_labels_angle,
            x_label_f=seconds_to_days,
        )

        plot_1d(
            tofs.to_value("second"),
            AOPs,
            "Argument of Periapsis Over Time",
            r"Time ($day$)",
            r"Argps ($rad$)",
            step=100,
            fill_under=False,
            remove_box=True,
            y_range=(min(AOPs) - 0.2, max(AOPs) + 0.2),
            custom_y_ticks=True,
            y_ticks=ticks_angle,
            y_tick_labels=tick_labels_angle,
            x_label_f=seconds_to_days,
        )

        plot_1d(
            tofs.to_value("second"),
            ECCs,
            "Eccentricity Over Time",
            r"Time ($day$)",
            r"Ecc",
            step=100,
            fill_under=False,
            remove_box=True,
            y_range=(min(ECCs) - 0.1, max(ECCs) + 0.1),
            scatter=True,
            x_label_f=seconds_to_days,
        )

        plot_1d(
            tofs.to_value("second"),
            INCs,
            "Inclination over time",
            r"Time ($day$)",
            r"Inc ($rad$)",
            step=100,
            fill_under=False,
            remove_box=True,
            y_range=(min(INCs) - 0.2, max(INCs) + 0.2),
            custom_y_ticks=True,
            y_ticks=ticks_angle,
            y_tick_labels=tick_labels_angle,
            x_label_f=seconds_to_days,
        )

        plot_1d(
            tofs.to_value("second"),
            altitudes,
            "Spacecraft Altitude Over Time",
            r"Time ($day$)",
            r"Altitude ($km$)",
            step=100,
            fill_under=False,
            remove_box=True,
            y_range=(min(altitudes) - 20, max(altitudes) + 20),
            x_label_f=seconds_to_days,
        )
