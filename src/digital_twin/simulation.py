"""Main file to run the simulation
"""

import time
from typing import Dict

from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.units import Quantity
import numpy as np

from digital_twin.constants import earth_R
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
    find_x_scale,
    plot_boolean_bars,
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

        # data to store
        modes = np.zeros(self.n_timesteps + 1)
        modes[0] = self.switch_algo.operating_mode

        battery_energies = np.zeros(self.n_timesteps + 1)
        battery_energies[0] = self.spacecraft.get_eps().get_battery_energy().value
        power_consumption = np.zeros(self.n_timesteps + 1)
        power_consumption[0] = self.spacecraft.get_eps().get_power_consumption().value
        power_generation = np.zeros(self.n_timesteps + 1)
        power_generation[0] = self.spacecraft.get_eps().get_power_generation().value

        data_storage = np.zeros(self.n_timesteps + 1)
        data_storage_GNSS_TOF = np.zeros(self.n_timesteps + 1)
        data_storage_HK = np.zeros(self.n_timesteps + 1)
        all, GNSS_TOF, HK = self.spacecraft.get_payload().get_data_storage()
        data_storage[0] = all.value
        data_storage_GNSS_TOF[0] = GNSS_TOF.value
        data_storage_HK[0] = HK.value

        vis_windows = np.zeros(
            (self.n_timesteps + 1, len(self.ground_stations))
        )  # for each ground station
        vis_windows[0] = self.propagator.calculate_vis_window(self.ground_stations)

        eclipse_windows = np.zeros(self.n_timesteps + 1)
        eclipse_windows[0] = self.propagator.calculate_eclipse_status()

        print("Number of timesteps:", self.n_timesteps)
        start_for_loop = time.time()

        for t in range(0, self.n_timesteps):

            # 1. propagate to next position and store the results
            rv = self.propagator.propagate(
                self.delta_t, self.spacecraft.C_D, self.spacecraft.A_over_m
            )
            eph[t + 1, :3] = rv[:3]
            eph[t + 1, 3:] = rv[3:]

            if t % 500 == 0:
                print(t)  # TODO: remove
                if np.linalg.norm(rv[:3] - earth_R.value) < 180:
                    print("Altitude decay: deorbiting")
                    break

            # 2. Calculate position based params: communication window, eclipse
            visibility = self.propagator.calculate_vis_window(self.ground_stations)
            eclipse_status = self.propagator.calculate_eclipse_status()

            # 3. Calculate user-scheduled params
            measurement_session = (
                True  # Assumption for now, later will depend on user-defined scheduler
            )
            measurement_session = self.spacecraft.get_payload().can_start_measuring(
                self.tofs[t + 1].to("second")
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

            # # 5. Ask spacecraft to update all subsystems
            self.spacecraft.update_subsystems(
                old_mode, new_mode, rv, com_window, eclipse_status, self.delta_t
            )

            # 6. Save data
            vis_windows[t + 1] = np.array([int(vis) for vis in visibility])
            eclipse_windows[t + 1] = int(eclipse_status)
            battery_energies[t + 1] = (
                self.spacecraft.get_eps().get_battery_energy().value
            )
            power_consumption[t + 1] = (
                self.spacecraft.get_eps().get_power_consumption().value
            )
            power_generation[t + 1] = (
                self.spacecraft.get_eps().get_power_generation().value
            )
            # data_storage[t + 1] = self.spacecraft.get_payload().get_data_storage().value
            all, GNSS_TOF, HK = self.spacecraft.get_payload().get_data_storage()
            data_storage[t + 1] = all.value
            data_storage_GNSS_TOF[t + 1] = GNSS_TOF.value
            data_storage_HK[t + 1] = HK.value

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
            "vis": vis_windows,
            "eclipse": eclipse_windows,
            "battery": battery_energies,
            "consumption": power_consumption,
            "generation": power_generation,
            "storage": data_storage,
            "storage_GNSS_TOF": data_storage_GNSS_TOF,
            "storage_HK": data_storage_HK,
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

        x_label, x_label_f = find_x_scale(self.duration_sim)
        step = int(len(data["tofs"]) / 100)

        if self.report_params["battery_energy"] == "yes":
            plot_1d(
                data["tofs"].to_value("second"),
                (data["battery"] * (u.W * u.s)).to(u.W * u.h),
                "Battery Energy Over Time",
                x_label,
                r"Battery Energy ($Wh$)",
                step=1,
                fill_under=False,
                remove_box=True,
                scatter=False,
                x_label_f=x_label_f,
                show=False,
                save_filename=self.report_params["folder"] + "battery_energy.pdf",
                markersize_plot=0,
            )

        if self.report_params["power_consumption"] == "yes":
            plot_1d(
                data["tofs"].to_value("second")[1:],
                (data["consumption"][1:] * (u.W * u.s)).to(u.W * u.h),
                "Power Consumption Over Time",
                x_label,
                r"Power Consumption ($Wh$)",
                step=1,  # no step because we want to capture small intervals of time
                fill_under=False,
                remove_box=True,
                scatter=False,
                x_label_f=x_label_f,
                show=False,
                save_filename=self.report_params["folder"] + "power_consumption.pdf",
                markersize_plot=0,
            )

        if self.report_params["power_generation"] == "yes":
            plot_1d(
                data["tofs"].to_value("second")[1:],
                (data["generation"][1:] * (u.W * u.s)).to(u.W * u.h),
                "Power Generation Over Time",
                x_label,
                r"Power Generation ($Wh$)",
                step=1,
                fill_under=False,
                remove_box=True,
                scatter=False,
                x_label_f=x_label_f,
                show=False,
                save_filename=self.report_params["folder"] + "power_generation.pdf",
                markersize_plot=0,
            )

        if self.report_params["data_storage"] == "yes":
            plot_1d(
                data["tofs"].to_value("second"),
                data["storage"],
                "Data Storage Over Time",
                x_label,
                r"Data Storage ($Mbit$)",
                step=step,
                fill_under=False,
                remove_box=True,
                scatter=False,
                x_label_f=x_label_f,
                show=False,
                save_filename=self.report_params["folder"] + "data_storage.pdf",
            )
            plot_1d(
                data["tofs"].to_value("second"),
                data["storage_GNSS_TOF"],
                "Data Storage Over Time (GNSS and TOF)",
                x_label,
                r"GNSS/TOF Data Storage ($Mbit$)",
                step=step,
                fill_under=False,
                remove_box=True,
                scatter=False,
                x_label_f=x_label_f,
                show=False,
                save_filename=self.report_params["folder"]
                + "data_storage_GNSS_TOF.pdf",
            )
            plot_1d(
                data["tofs"].to_value("second"),
                data["storage_HK"],
                "Data Storage Over Time (Housekeeping Data)",
                x_label,
                r"HK Data Storage ($Mbit$)",
                step=step,
                fill_under=False,
                remove_box=True,
                scatter=False,
                x_label_f=x_label_f,
                show=False,
                save_filename=self.report_params["folder"] + "data_storage_HK.pdf",
            )

        if self.report_params["visibility_windows"] == "yes":
            save_filename = self.report_params["folder"] + "visibility_windows.pdf"
            plot_boolean_bars(
                data["vis"],
                data["tofs"].to_value("second"),
                self.duration_sim,
                "Visibility Windows",
                save_filename=save_filename,
                show=False,
            )

        if self.report_params["eclipse_windows"] == "yes":
            save_filename = self.report_params["folder"] + "eclipse_windows.pdf"
            plot_boolean_bars(
                data["eclipse"],
                data["tofs"].to_value("second"),
                self.duration_sim,
                "Eclipse Windows",
                save_filename=save_filename,
                show=False,
            )
        if self.report_params["save_np_arrays_telecom_data"] == "yes":
            folder = self.report_params["folder"]
            with open(folder + "times.npy", "wb") as f:
                np.save(f, data["tofs"].to_value("second"))
            with open(folder + "visibility.npy", "wb") as f:
                np.save(f, data["vis"])
            with open(folder + "modes.npy", "wb") as f:
                np.save(f, data["modes"])

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

    x_label, x_label_f = find_x_scale(duration_sim)

    # update how many values to step to obtain clear plotting (approximately 100 points)
    step = int(len(tofs) / 100)

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
