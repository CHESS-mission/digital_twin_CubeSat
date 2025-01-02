"""Definition of the Simulation class which acts as the central manager for orchestrating the various components
of a satellite mission simulation. It brings together the orbit propagator, spacecraft, ground station, and other
subsystems to run a cohesive simulation based on user-defined parameters.
"""

import time
from typing import Any

from astropy import units as u
from astropy.time import Time, TimeDelta
import numpy as np

from digital_twin.constants import earth_R, simulation_unit, simulation_unit_string
from digital_twin.ground_station import GroundStation
from digital_twin.mode_switch import ModeSwitch
from digital_twin.orbit_propagator import OrbitPropagator
from digital_twin.report import produce_report
from digital_twin.spacecraft import Spacecraft
from digital_twin.utils import (
    get_astropy_unit_time,
    extract_propagation_data_from_ephemeris,
)


class Simulation:
    """Manager class which gathers simulation objects and runs the main simulation loop."""

    def __init__(
        self,
        simulation_params: dict,
        orbit_params: dict,
        spacecraft_params: dict,
        station_params: dict,
        mission_design_params: dict,
    ) -> None:

        # TIME INITIALIZATION
        self.sim_unit_string = simulation_unit_string
        self.sim_unit = simulation_unit
        self.delta_t = (
            simulation_params["delta_t"]
            * get_astropy_unit_time(simulation_params["units_delta_t"])
        ).to(self.sim_unit)
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
        )  # Gives the results in days
        self.epochs_array = np.array(epoch + self.tofs)

        # MODE SWITCH ALGORITHM INITIALIZATION
        self.switch_algo = ModeSwitch(
            init_mode=int(spacecraft_params["general"]["init_operating_mode"])
        )

        # SPACECRAFT INITIALIZATION
        self.spacecraft = Spacecraft(spacecraft_params, self.switch_algo.operating_mode)

        # GROUND STATION INITIALIZATION
        self.ground_stations = []
        for param in station_params["stations"]:
            self.ground_stations.append(GroundStation(param))
        self.ground_stations = np.array(self.ground_stations)

        # PROPAGATOR INITIALIZATION
        update_air_density_timestep = simulation_params[
            "update_air_density_timestep"
        ] * get_astropy_unit_time(simulation_params["update_air_density_timestep_unit"])
        atmosphere_model = simulation_params["atmosphere_model"]
        self.propagator = OrbitPropagator(
            orbit_params, epoch, atmosphere_model, update_air_density_timestep
        )

        # Additional user input
        user_input = mission_design_params["user_input"]
        self.handle_user_input(user_input)

        # Print user parameters
        self.report_params = mission_design_params["report"]
        self.print_parameters()

    def run(self) -> None:
        """Function to run the simulation, which contains the main simulation loop."""
        print("Simulation running...")

        # INITIALIZATION
        eph = np.zeros((self.n_timesteps + 1, 6))
        eph[0, :3] = self.propagator.r
        eph[0, 3:] = self.propagator.v

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
        all, GNSS_TOF, HK = self.spacecraft.get_obc().get_data()
        data_storage[0] = all.value
        data_storage_GNSS_TOF[0] = GNSS_TOF.value
        data_storage_HK[0] = HK.value

        vis_windows = np.zeros(
            (self.n_timesteps + 1, len(self.ground_stations))
        )  # for each ground station
        vis, _ = self.propagator.calculate_vis_window(self.ground_stations)
        vis_windows[0] = vis

        eclipse_windows = np.zeros(self.n_timesteps + 1)
        eclipse, _ = self.propagator.calculate_eclipse_status()
        eclipse_windows[0] = eclipse

        density_array = np.zeros(self.n_timesteps + 1)
        density_array[0] = self.propagator.get_density().value

        print("Number of timesteps:", self.n_timesteps)

        # MAIN SIMULATION LOOP
        start_for_loop = time.time()
        for t in range(0, self.n_timesteps):

            # 1. propagate to next position and store the results
            try:
                rv = self.propagator.propagate(
                    self.delta_t, self.spacecraft.C_D, self.spacecraft.A_over_m
                )
            except (
                RuntimeError,
                ValueError,
            ):  # Spacecraft cannot be propagated anymore (usually because altitude is too low)
                break
            except ZeroDivisionError:
                break

            eph[t + 1, :3] = rv[:3]
            eph[t + 1, 3:] = rv[3:]

            if t % 1000 == 0:  # For debugging
                print("iter ", t)
                if np.linalg.norm(rv[:3] - earth_R.value) < 180:
                    print("Altitude decay: deorbiting")
                    break

            # # 2. Calculate position based params: communication window, eclipse status
            # visibility, gs_coords_array = self.propagator.calculate_vis_window(
            #     self.ground_stations
            # )
            # eclipse_status, r_earth_sun = self.propagator.calculate_eclipse_status()
            # r_earth_sun = r_earth_sun

            # # 3. Calculate user-scheduled params
            # measurement_session = self.spacecraft.get_payload().can_start_measuring(
            #     self.tofs[t + 1].to("second")
            # )  # Currently implemented with a user parameter deciding the maximum number of measurement sessions per day
            # com_window = False
            # gs_coords = None
            # for i, vis in enumerate(visibility):
            #     if vis:
            #         com_window = True
            #         gs_coords = (gs_coords_array[i] * u.km).flatten()
            #         break  # Right now, we only consider the first ground station which is visible from satellite

            # # 4. Check for potential flags raised by OBC
            # safe_flag = self.spacecraft.get_obc().raise_spacecraft_safe_flag()

            # # 5. Switch mode based on current state
            # old_mode = self.switch_algo.operating_mode
            # self.switch_algo.switch_mode(
            #     self.spacecraft.get_eps(),
            #     self.spacecraft.get_telecom(),
            #     self.spacecraft.get_payload(),
            #     self.spacecraft.get_data_storage(),
            #     com_window,
            #     eclipse_status,
            #     measurement_session,
            #     safe_flag,
            # )
            # new_mode = self.switch_algo.operating_mode
            # modes[t + 1] = new_mode

            # # 6. Ask spacecraft to update all subsystems
            # self.spacecraft.update_subsystems(
            #     old_mode,
            #     new_mode,
            #     rv,
            #     com_window,
            #     eclipse_status,
            #     self.delta_t,
            #     r_earth_sun,
            #     gs_coords,
            # )

            # # 7. Save data at current timestep
            # vis_windows[t + 1] = np.array([int(vis) for vis in visibility])
            # eclipse_windows[t + 1] = int(eclipse_status)
            # battery_energies[t + 1] = (
            #     self.spacecraft.get_eps().get_battery_energy().value
            # )
            # power_consumption[t + 1] = (
            #     self.spacecraft.get_eps().get_power_consumption().value
            # )
            # power_generation[t + 1] = (
            #     self.spacecraft.get_eps().get_power_generation().value
            # )
            # all, GNSS_TOF, HK = self.spacecraft.get_obc().get_data()
            # data_storage[t + 1] = all.value
            # data_storage_GNSS_TOF[t + 1] = GNSS_TOF.value
            # data_storage_HK[t + 1] = HK.value

            density_array[t + 1] = self.propagator.get_density().value

        end_for_loop = time.time()
        duration = end_for_loop - start_for_loop
        print("Simulation ended!")
        print(f"Time: {duration} s")

        # If the simulation finished earlier than planned
        if t < self.n_timesteps - 1:
            print("WARNING: simulation couldn't end because of propagation error")
            eph = eph[: t + 1]

        # Extract the orbital elements data
        rr, vv, SMAs, ECCs, INCs, RAANs, AOPs, TAs, altitudes = (
            extract_propagation_data_from_ephemeris(eph)
        )

        # Remove data if propagation continued below altitude 0 or stopped earlier than planned
        neg_index = np.argmax(altitudes < 0)
        if neg_index != 0 or t < self.n_timesteps - 1:
            last_ind = neg_index if neg_index != 0 else t

            altitudes = altitudes[:last_ind]
            rr = rr[:last_ind]
            vv = vv[:last_ind]
            SMAs = SMAs[:last_ind]
            ECCs = ECCs[:last_ind]
            INCs = INCs[:last_ind]
            RAANs = RAANs[:last_ind]
            AOPs = AOPs[:last_ind]
            TAs = TAs[:last_ind]

            eph = eph[:last_ind]
            vis_windows = vis_windows[:last_ind]
            modes = modes[:last_ind]
            battery_energies = battery_energies[:last_ind]
            power_consumption = power_consumption[:last_ind]
            power_generation = power_generation[:last_ind]
            data_storage = data_storage[:last_ind]
            data_storage_GNSS_TOF = data_storage_GNSS_TOF[:last_ind]
            data_storage_HK = data_storage_HK[:last_ind]
            eclipse_windows = eclipse_windows[:last_ind]
            density_array = density_array[:last_ind]

            self.duration_sim = last_ind * self.delta_t
            self.tofs = self.tofs[:last_ind]
            self.epochs_array = self.epochs_array[:last_ind]

        # Extract orbit and spacecraft states
        orbit_state = self.propagator.save_state()
        spacecraft_state = self.spacecraft.save_state()
        spacecraft_state["general"][
            "init_operating_mode"
        ] = self.switch_algo.operating_mode

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
            "duration_sim": self.duration_sim,
            "epochs_array": self.epochs_array,
            "initial_orbit": self.propagator.get_initial_orbit(),
            "ground_stations": self.ground_stations,
            "orbit_state": orbit_state,
            "spacecraft_state": spacecraft_state,
            "density_array": density_array,
        }

        # Produce report
        produce_report(data_results, self.report_params)

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

    def handle_user_input(self, user_input: dict) -> None:
        """Handle additional user input. Currently only implemented for uplink safe mode trigger.

        Args:
            user_input (dict): Dictionary containing the user input to consider.
        """
        if user_input["uplink_safe_mode"]:  # If dic is not empty
            self.spacecraft.get_telecom().add_uplink_safe_mode(
                user_input["uplink_safe_mode"]
            )
