"""Main file to run the simulation
"""

import time
from typing import Dict

from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.units import Quantity
import numpy as np

from digital_twin.ground_station import GroundStation
from digital_twin.mode_switch import ModeSwitch
from digital_twin.orbit_propagator import OrbitPropagator
from digital_twin.orbit_propagator.constants import attractor_string
from digital_twin.plotting import (
    plot_1d,
    seconds_to_days,
    seconds_to_minutes,
    seconds_to_hours,
    plot_orbit_trajectory_3d,
    plot_orbit_2d,
    plot_groundtrack,
    plot_operating_modes,
)
from digital_twin.spacecraft import Spacecraft
from digital_twin.utils import (
    get_astropy_unit_time,
    check_and_empty_folder,
    extract_propagation_data_from_ephemeris,
)


class Simulation:
    """Manager class which gathers simulation objects and runs the main simulation loop."""

    def __init__(
        self,
        simulation_params: Dict,
        orbit_params: Dict,
        spacecraft_params: Dict,
        station_params: Dict,
    ) -> None:

        # TIME INITIALIZATION
        if simulation_params["units_delta_t"] != "second":
            raise ValueError(
                "Simulation no implemented for timesteps not expressed in seconds"
            )
        self.sim_unit = get_astropy_unit_time(simulation_params["units_delta_t"])
        self.delta_t = simulation_params["delta_t"] * self.sim_unit
        init_time_local = 0 * self.sim_unit

        self.duration_sim = simulation_params["duration_sim"] * get_astropy_unit_time(
            simulation_params["units_duration_sim"]
        )
        self.n_timesteps = int(
            self.duration_sim.to_value(self.sim_unit) / self.delta_t.value
        )
        end_time_local = (self.n_timesteps * self.delta_t.value) * self.sim_unit

        epoch = Time(orbit_params["epoch"], format="iso", scale="utc")
        self.init_time = init_time_local + epoch
        self.end_time = end_time_local + epoch

        self.tofs = TimeDelta(
            np.linspace(init_time_local, end_time_local, num=self.n_timesteps + 1)
        )  # gives the results in days
        self.epochs_array = np.array(epoch + self.tofs)

        # MODE SWITCH ALGORITHM INITIALIZATION
        self.switch_algo = ModeSwitch(
            init_mode=simulation_params["init_operating_mode"]
        )

        # SPACECRAFT INITIALIZATION
        self.spacecraft = Spacecraft(spacecraft_params, self.switch_algo.operating_mode)

        # GROUND STATION INITIALIZATION
        self.ground_stations = []
        for param in station_params["stations"]:
            self.ground_stations.append(GroundStation(param))
        self.ground_stations = np.array(self.ground_stations)

        # PROPAGATOR INITIALIZATION
        self.propagator = OrbitPropagator(orbit_params, epoch)

        # data
        self.report_params = simulation_params["report"]
        self.print_parameters()

    def run(self) -> None:
        """Main simulation loop."""
        print("Simulation running...")
        eph = np.zeros((self.n_timesteps + 1, 6))
        eph[0, :3] = self.propagator.r
        eph[0, 3:] = self.propagator.v
        modes = np.zeros(self.n_timesteps + 1)
        modes[0] = self.switch_algo.operating_mode

        print("Number of timesteps:", self.n_timesteps)
        start_for_loop = time.time()

        for t in range(0, self.n_timesteps):
            if t % 500 == 0:
                print(t)  # TODO: remove

            # 1. propagate to next position and store the results
            rv = self.propagator.propagate(
                self.delta_t, self.spacecraft.C_D, self.spacecraft.A_over_m
            )
            eph[t + 1, :3] = rv[:3]
            eph[t + 1, 3:] = rv[3:]

            # 2. Calculate position based params: communication window, eclipse
            visibility = self.propagator.calculate_vis_window(self.ground_stations)
            eclipse_status = self.propagator.calculate_eclipse_status()

            # 3. Calculate user-scheduled params
            measurement_session = (
                True  # Assumption for now, later will depend on user-defined scheduler
            )
            com_window = visibility  # Assumption for now, later might depend on user-defined scheduler

            # 3. Check for potential flags raised by OBS
            safe_flag = False

            # 4. Switch mode based on location, Eps and Telecom states
            old_mode = self.switch_algo.operating_mode
            self.switch_algo.switch_mode(
                self.spacecraft.get_eps(),
                self.spacecraft.get_telecom(),
                self.spacecraft.get_payload(),
                com_window,
                eclipse_status,
                measurement_session,
                safe_flag,
            )
            new_mode = self.switch_algo.operating_mode
            modes[t + 1] = new_mode

            # 5. Ask spacecraft to update all subsystems
            # self.spacecraft.update_subsystems(
            #     old_mode, new_mode, rv, com_window, eclipse_status, self.delta_t

        end_for_loop = time.time()
        duration = end_for_loop - start_for_loop
        print("Simulation ended!")
        print(f"Time: {duration} s")

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
            "modes": modes,
        }
        # Produce report
        self.produce_report(data_results)

    def produce_report(self, data: Dict) -> None:
        """Produce plots and other report data."""
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
                self.duration_sim,
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

        stations_coords = []
        stations_names = []
        for station in self.ground_stations:
            name, pos = station.get_name_pos()
            stations_coords.append(pos)
            stations_names.append(name)
        if self.report_params["groundtrack"] == "yes":
            plot_groundtrack(
                self.report_params["title_figures"],
                data["rr"],
                self.epochs_array,
                "CHESS_1 cubeSat",
                self.report_params["folder"],
                stations_coords=np.array(stations_coords),
                stations_name=np.array(stations_names),
            )

        if self.report_params["modes"] == "yes":
            save_filename = self.report_params["folder"] + "modes.pdf"
            plot_operating_modes(
                data["modes"],
                data["tofs"].to_value("second"),
                self.duration_sim,
                save_filename=save_filename,
                show=False,
            )

    def print_parameters(self) -> None:
        """Print a summary of the the simulation objects."""
        print("")
        print("*******************")
        print("INITIAL PARAMETERS:\n")
        print("Simulation:")
        print(f"- init time: {self.init_time}")
        print(f"- end time: {self.end_time}")
        print(f"- init operating mode: {self.switch_algo.print_operating_mode()}")
        print("")
        print(str(self.spacecraft))
        print("\nGround Stations:")
        for ground_station in self.ground_stations:
            print("- " + str(ground_station))
        print("")
        print(str(self.propagator))
        print("*******************")
        print("")


# Method not part of the class but not general enough to be in utils file
def plot_orbital_elem_evolution(
    tofs: np.ndarray,
    RAANs: np.ndarray,
    AOPs: np.ndarray,
    ECCs: np.ndarray,
    INCs: np.ndarray,
    altitudes: np.ndarray,
    folder: str,
    duration_sim: Quantity["time"],
):
    """Plot the evolution of RAAN, AOP, ECC, INC, and altitude during the simulation."""
    ticks_angle = np.array([0, 1 / 2 * np.pi, np.pi, 3 / 2 * np.pi, 2 * np.pi])
    tick_labels_angle = np.array([r"$0$", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])

    # update what the x_axis scale should be
    if duration_sim <= 3 * u.h:
        x_label_f = seconds_to_minutes
        x_label = r"Time ($min$)"
    elif duration_sim <= 3 * u.day:
        x_label_f = seconds_to_hours
        x_label = r"Time ($hour$)"
    else:
        x_label_f = seconds_to_days
        x_label = r"Time ($day$)"

    # update how many values to step to obtain clear plotting (approximately 50 points)
    step = int(len(tofs) / 50)

    plot_1d(
        tofs.to_value("second"),
        RAANs,
        "RAAN Over Time",
        x_label,
        r"RAAN ($rad$)",
        step=step,
        fill_under=False,
        remove_box=True,
        y_range=(min(RAANs) - 0.2, max(RAANs) + 0.2),
        scatter=True,
        custom_y_ticks=True,
        y_ticks=ticks_angle,
        y_tick_labels=tick_labels_angle,
        x_label_f=x_label_f,
        show=False,
        save_filename=folder + "RAAN.pdf",
    )

    plot_1d(
        tofs.to_value("second"),
        AOPs,
        "Argument of Periapsis Over Time",
        x_label,
        r"Argps ($rad$)",
        step=step,
        fill_under=False,
        remove_box=True,
        y_range=(min(AOPs) - 0.2, max(AOPs) + 0.2),
        custom_y_ticks=True,
        y_ticks=ticks_angle,
        y_tick_labels=tick_labels_angle,
        x_label_f=x_label_f,
        show=False,
        save_filename=folder + "AOP.pdf",
    )

    plot_1d(
        tofs.to_value("second"),
        ECCs,
        "Eccentricity Over Time",
        x_label,
        r"Ecc",
        step=step,
        fill_under=False,
        remove_box=True,
        y_range=(min(ECCs) - 0.1, max(ECCs) + 0.1),
        scatter=True,
        x_label_f=x_label_f,
        show=False,
        save_filename=folder + "ECC.pdf",
    )

    plot_1d(
        tofs.to_value("second"),
        INCs,
        "Inclination over time",
        x_label,
        r"Inc ($rad$)",
        step=step,
        fill_under=False,
        remove_box=True,
        y_range=(min(INCs) - 0.2, max(INCs) + 0.2),
        custom_y_ticks=True,
        y_ticks=ticks_angle,
        y_tick_labels=tick_labels_angle,
        x_label_f=x_label_f,
        show=False,
        save_filename=folder + "INC.pdf",
    )

    plot_1d(
        tofs.to_value("second"),
        altitudes,
        "Spacecraft Altitude Over Time",
        x_label,
        r"Altitude ($km$)",
        step=step,
        fill_under=False,
        remove_box=True,
        y_range=(min(altitudes) - 20, max(altitudes) + 20),
        x_label_f=x_label_f,
        show=False,
        save_filename=folder + "altitude.pdf",
    )
