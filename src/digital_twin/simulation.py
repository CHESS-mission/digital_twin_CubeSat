"""Main file that runs the simulation
"""

from typing import Dict
import time

import numpy as np

from astropy.time import Time, TimeDelta

from digital_twin.orbit_propagator import OrbitPropagator
from digital_twin.orbit_propagator.constants import attractor_string
from digital_twin.spacecraft import Spacecraft
from digital_twin.mode_switch import ModeSwitch
from digital_twin.utils import (
    get_astropy_unit,
    check_and_empty_folder,
    extract_propagation_data_from_ephemeris,
)
from digital_twin.plotting import (
    plot_1d,
    seconds_to_days,
    plot_orbit_trajectory_3d,
    plot_orbit_2d,
    plot_groundtrack,
)


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
        self.epochs_array = np.array(self.epoch + self.tofs)

        # INITIALIZE MODE SWITCH ALGORITHM
        self.switch_algo = ModeSwitch(
            init_mode=simulation_params["init_operating_mode"]
        )
        print(f"Init operating mode: {self.switch_algo.operating_mode}")

        # report parameters
        self.report_params = simulation_params["report"]

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

            # 1. propagate to next position and store the results
            rv = self.propagator.propagate(self.delta_t)
            eph[t + 1, :3] = rv[:3]
            eph[t + 1, 3:] = rv[3:]

            # 2. Calculate position based params: communication window, eclipse
            # com_window = self.propagator.calculate_com_window()
            # eclipse_status = self.propagator.calculate_eclipse_status()
            # measurement_session = (
            # True  # initiate this way, always want to measure if possible
            # )

            # 3. Check for potential flags raised by OBS
            # safe_flag = False

            # 4. Switch mode based on location, Eps and Telecom states
            # self.switch_algo.switch_mode(
            #     self.spacecraft.eps,
            #     self.spacecraft.telecom,
            #     com_window,
            #    eclipse_status,
            #    measurement_session,
            #     safe_flag,
            # )

        end_for_loop = time.time()
        duration_for_loop = end_for_loop - start_for_loop
        print("Time: ", duration_for_loop)

        # Extract the data
        rr, vv, SMAs, ECCs, INCs, RAANs, AOPs, TAs, altitudes = (
            extract_propagation_data_from_ephemeris(eph)
        )
        data_results = {
            "tofs": self.tofs,
            "rr": rr,
            "vv": vv,
            "SMAs": SMAs,
            "ECCs": ECCs,
            "INCs": INCs,
            "RAANs": RAANs,
            "AOPs": AOPs,
            "TAs": TAs,
            "altitudes": altitudes,
        }
        # Produce report
        self.produce_report(data_results)

    def produce_report(self, data: str):
        check_and_empty_folder(self.report_params["folder"])
        if self.report_params["orbital_elem_evolution"] == "yes":
            plot_orbital_elem_evolution(
                data["tofs"],
                data["RAANs"],
                data["AOPs"],
                data["ECCs"],
                data["INCs"],
                data["altitudes"],
                self.report_params["folder"],
            )
        if self.report_params["trajectory_2d"] == "yes":
            plot_orbit_2d(
                self.report_params["title_figures"],
                attractor_string,
                self.report_params["folder"],
                "initial orbit",
                orbit=self.propagator.get_initial_orbit(),
            )
        if self.report_params["trajectory_3d"] == "yes":
            plot_orbit_trajectory_3d(
                self.report_params["title_figures"],
                attractor_string,
                self.report_params["folder"],
                orbit=self.propagator.get_initial_orbit(),
                label_orbit="initial orbit",
                label_traj="trajectory",
                traj=data["rr"],
            )

        if self.report_params["groundtrack"] == "yes":
            plot_groundtrack(
                self.report_params["title_figures"],
                data["rr"],
                self.epochs_array,
                "CHESS_1 cubeSat",
                self.report_params["folder"],
                station_coords=np.array([46.31, 6.38]),
                station_name="Lausanne",
            )


# Methods not part of the class not not general enough to be in utils file
def plot_orbital_elem_evolution(
    tofs: np.array,
    RAANs: np.array,
    AOPs: np.array,
    ECCs: np.array,
    INCs: np.array,
    altitudes: np.array,
    folder: str,
):
    ticks_angle = np.array([0, 1 / 2 * np.pi, np.pi, 3 / 2 * np.pi, 2 * np.pi])
    tick_labels_angle = np.array([r"$0$", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
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
        show=False,
        save_filename=folder + "RAAN.pdf",
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
        show=False,
        save_filename=folder + "AOP.pdf",
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
        show=False,
        save_filename=folder + "ECC.pdf",
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
        show=False,
        save_filename=folder + "INC.pdf",
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
        show=False,
        save_filename=folder + "altitude.pdf",
    )
