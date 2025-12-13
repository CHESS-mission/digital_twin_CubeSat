"""Definition of the Simulation class which acts as the central manager for orchestrating the various components
of a satellite mission simulation. It brings together the orbit propagator, spacecraft, ground station, and other
subsystems to run a cohesive simulation based on user-defined parameters.
"""

import time
from typing import Any
from queue import Queue

from astropy import units as u
from astropy.time import Time, TimeDelta
import numpy as np

from digital_twin.constants import earth_R, simulation_unit, simulation_unit_string
from digital_twin.ground_station import GroundStation
from digital_twin.mode_switch import ModeSwitch
from digital_twin.commands import CommandProcessor

from digital_twin.orbit_propagator import OrbitPropagator
from digital_twin.report import produce_report
from digital_twin.spacecraft import Spacecraft
from digital_twin.utils import (
    get_astropy_unit_time,
    extract_propagation_data_from_ephemeris,
    convert_cartesian_to_spherical,
)
import influxdb_client, os
from influxdb_client.client.write_api import SYNCHRONOUS
from datetime import datetime
from dotenv import load_dotenv

class Simulation:
    """Manager class which gathers simulation objects and runs the main simulation loop."""

    def __init__(
        self,
        simulation_params: dict,
        orbit_params: dict,
        spacecraft_params: dict,
        station_params: dict,
        mission_design_params: dict,
        env_file: str,
    ) -> None:
        self.verbose = True if simulation_params["verbose"] == "yes" else False

        # TIME INITIALIZATION
        self.sim_unit_string = simulation_unit_string
        self.sim_unit = simulation_unit
        self.delta_t = (
            simulation_params["delta_t"]
            * get_astropy_unit_time(simulation_params["delta_t_unit"])
        ).to(self.sim_unit)
        init_time_local = 0 * self.sim_unit

        self.duration_sim = simulation_params["duration_sim"] * get_astropy_unit_time(
            simulation_params["duration_sim_unit"]
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
        self.epochs_array = epoch + self.tofs

        self.run_real_time = simulation_params.get("run_real_time", False)
        print(f"ATTENTION: Running simulation in real-time mode (this is very inefficient!)") if self.run_real_time and self.verbose else None

        self.command_queue = Queue()
        self.command_processor = CommandProcessor(simulation=self)

        # MODE SWITCH ALGORITHM INITIALIZATION
        self.switch_algo = ModeSwitch(
            init_mode=int(spacecraft_params["general"]["init_operating_mode"]),
            verbose=self.verbose,
        )

        # SPACECRAFT INITIALIZATION
        self.spacecraft = Spacecraft(
            spacecraft_params, self.switch_algo.operating_mode, self.verbose
        )

        # GROUND STATION INITIALIZATION
        self.ground_stations = []
        for param in station_params["stations"]:
            self.ground_stations.append(GroundStation(param, self.verbose))
        self.ground_stations = np.array(self.ground_stations)

        # PROPAGATOR INITIALIZATION
        update_air_density_timestep = simulation_params[
            "update_air_density_timestep"
        ] * get_astropy_unit_time(simulation_params["update_air_density_timestep_unit"])
        atmosphere_model = simulation_params["atmosphere_model"]
        self.propagator = OrbitPropagator(
            orbit_params,
            epoch,
            atmosphere_model,
            update_air_density_timestep,
            self.verbose,
        )
        self.propagation_only = (
            True if simulation_params["propagation_only"] == "yes" else False
        )

        # Additional user input
        user_input = mission_design_params["user_input"]
        if user_input.get("uplink_safe_mode"):
            self.send_command("uplink_safe_mode", user_input["uplink_safe_mode"])

        # Report and printing
        self.report_params = mission_design_params["report"]
        if simulation_params["print_initial_parameters"] == "yes":
            self.print_parameters()  # Print user parameters

        # InfluxDB initialization
        if simulation_params.get("influxdb_delta_t", -1) > 0:
            self.influxdb_delta_t = (simulation_params["influxdb_delta_t"] * get_astropy_unit_time(simulation_params.get("influxdb_delta_t_unit", simulation_params["delta_t_unit"]))).to(self.sim_unit) 
            self.influxdb_only_visible = simulation_params.get("influxdb_only_visible", False)

            #Loading InfluxDB credentials from .env file
            load_dotenv(env_file)
            token = os.environ.get("INFLUXDB_TOKEN")
            url = os.environ.get("INFLUXDB_URL", "http://localhost:8086")  # Default URL
            self.influxdborg = os.environ.get("INFLUXDB_ORG", "EST")  # Default to "EST" if not set
            self.influxdbbucket = os.environ.get("INFLUXDB_BUCKET", "NICE")         
            if not token:
                raise ValueError("InfluxDB token not found in environment variables.")

            client = influxdb_client.InfluxDBClient(url=url, token=token, org=self.influxdborg)

            # Delete all previous data from bucket
            start = "1970-01-01T00:00:00Z"
            stop =  datetime(2070, 1, 1, 0, 0, 0)
            client.delete_api().delete(start, stop, '', bucket=self.influxdbbucket, org=self.influxdborg)

            self.write_api = client.write_api(write_options=SYNCHRONOUS)
            self.influxdb_last_upload_t = -1
        else:
            self.influxdb_delta_t = -1 * self.sim_unit

    def run(self, results_folder: str = "results/") -> None:
        """Function to run the simulation, which contains the main simulation loop."""
        print("Simulation running...") if self.verbose else None

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
        solar_cells_efficiency = np.zeros(self.n_timesteps + 1)
        solar_cells_efficiency[0] = self.spacecraft.get_eps().get_solar_cells_efficiency().value
        
        data_storage = np.zeros(self.n_timesteps + 1)
        data_storage_payload = np.zeros(self.n_timesteps + 1)
        data_storage_HK = np.zeros(self.n_timesteps + 1)
        all, payload, HK = self.spacecraft.get_obc().get_data()
        data_storage[0] = all.value
        data_storage_payload[0] = payload.value
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

        print("Number of timesteps:", self.n_timesteps) if self.verbose else None

        # MAIN SIMULATION LOOP
        start_for_loop = time.time()
        step_10_percent = max(1, self.n_timesteps // 10)
        for t in range(0, self.n_timesteps):

            # Process commands            
            self._process_commands()

            if self.run_real_time:
                # Calculate the target time for the current timestep
                target_time = start_for_loop + (t + 1) * self.delta_t.to_value(u.second)
                # Calculate the time to wait until the target time
                now = time.time()
                time_to_wait = (target_time - now)
                if time_to_wait > 0:
                    time.sleep(time_to_wait)


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

            if t % step_10_percent == 0 and self.verbose:
                percent = int((t / self.n_timesteps) * 100)
                print(f"> iter {t} ({percent}%)")

            # For simulations where only propagation matters, skip the next steps
            if not self.propagation_only:

                # 2. Calculate position based params: communication window, eclipse status
                visibility, gs_coords_array = self.propagator.calculate_vis_window(
                    self.ground_stations
                )
                eclipse_status, r_earth_sun = self.propagator.calculate_eclipse_status()

                # 3. Calculate user-scheduled params
                measurement_session = self.spacecraft.get_payload().can_start_measuring(
                    self.tofs[t + 1].to("second")
                )  # Currently implemented with a user parameter deciding the maximum number of measurement sessions per day
                com_window = False
                gs_coords = None
                for i, vis in enumerate(visibility):
                    if vis:
                        com_window = True
                        gs_coords = (gs_coords_array[i] * u.km).flatten()
                        break  # Right now, we only consider the first ground station which is visible from satellite

                # 4. Check for potential flags raised by OBC
                safe_flag = self.spacecraft.get_obc().raise_spacecraft_safe_flag()

                # 5. Switch mode based on current state
                old_mode = self.switch_algo.operating_mode
                self.switch_algo.switch_mode(
                    self.spacecraft.get_eps(),
                    self.spacecraft.get_telecom(),
                    self.spacecraft.get_payload(),
                    self.spacecraft.get_data_storage(),
                    com_window,
                    eclipse_status,
                    measurement_session,
                    safe_flag,
                )
                new_mode = self.switch_algo.operating_mode
                modes[t + 1] = new_mode

                # 6. Ask spacecraft to update all subsystems
                self.spacecraft.update_subsystems(
                    old_mode,
                    new_mode,
                    rv,
                    com_window,
                    eclipse_status,
                    self.delta_t,
                    r_earth_sun,
                    gs_coords,
                    t
                )

                # 7. Save data at current timestep
                vis_windows[t + 1] = visibility.astype(int)
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
                solar_cells_efficiency[t + 1] = (
                    self.spacecraft.get_eps().get_solar_cells_efficiency().value
                )
                
                all, payload, HK = self.spacecraft.get_obc().get_data()
                data_storage[t + 1] = all.value
                data_storage_payload[t + 1] = payload.value
                data_storage_HK[t + 1] = HK.value

                density_array[t + 1] = self.propagator.get_density().value


                # 8. If live InfluxDB upload is activated, send data to database
                if self.influxdb_delta_t > 0:
                    # Check if enough simulation time has passed since the last upload or we're at the final timestep
                    last_tof = self.tofs[self.influxdb_last_upload_t] if self.influxdb_last_upload_t >=0 else 0
                    elapsed_since_last = (self.tofs[t] - last_tof).to(self.sim_unit)
                    enough_time_elapsed = elapsed_since_last >= self.influxdb_delta_t
                    is_final_timestep = (t == self.n_timesteps - 1)

                    is_visible_check = np.any(vis_windows[t + 1]) if self.influxdb_only_visible else True
                    if (enough_time_elapsed and is_visible_check) or is_final_timestep:

                        print("Uploading data to InfluxDB...") if self.verbose else None

                        # Only compute data for the new points since last upload
                        start_idx = self.influxdb_last_upload_t + 1
                        end_idx = t + 1
                        n_new_points = end_idx - start_idx

                        # Extract orbital elements
                        rr_new, vv_new, SMAs_new, ECCs_new, INCs_new, RAANs_new, AOPs_new, TAs_new, altitudes_new = (
                            extract_propagation_data_from_ephemeris(eph[start_idx:end_idx])
                        )
                        
                        # Convert coordinates
                        latitude_deg_new, longitude_deg_new = convert_cartesian_to_spherical(rr_new, self.epochs_array[start_idx:end_idx])
                        
                        # Handle visibility data
                        vis_slice = vis_windows[start_idx:end_idx]
                        if vis_slice.ndim > 1:  # Multiple ground stations
                            vis_sums = np.sum(vis_slice, axis=1).astype(int)
                            is_visible = (vis_sums > 0).astype(int)
                        else:  # Single ground station or already summed
                            vis_sums = vis_slice.astype(int)
                            is_visible = vis_sums

                        timestamps = self.epochs_array[start_idx:end_idx].to_datetime()

                        # Extract all data slices as numpy arrays
                        data_arrays = {
                            'modes': modes[start_idx:end_idx].astype(int),
                            'storage': data_storage[start_idx:end_idx],
                            'storage_payload': data_storage_payload[start_idx:end_idx],
                            'storage_hk': data_storage_HK[start_idx:end_idx],
                            'battery': battery_energies[start_idx:end_idx],
                            'consumption': power_consumption[start_idx:end_idx],
                            'generation': power_generation[start_idx:end_idx],
                            'eclipse': eclipse_windows[start_idx:end_idx].astype(int),
                            'density': density_array[start_idx:end_idx],
                            'solar_eff': solar_cells_efficiency[start_idx:end_idx],
                            'altitudes': altitudes_new,
                            'RAANs': RAANs_new,
                            'AOPs': AOPs_new,
                            'ECCs': ECCs_new,
                            'INCs': INCs_new,
                            'latitude': latitude_deg_new,
                            'longitude': longitude_deg_new,
                            'vis_sums': vis_sums,
                            'is_visible': is_visible
                        }
                        # Use dictionary comprehension for batch creation (more efficient than explicit loop)
                        batch_data = [
                            {
                                "measurement": "satellite_data",
                                "tags": {"mode": int(data_arrays['modes'][i]), "visible": int(data_arrays['is_visible'][i])},
                                "fields": {
                                    "visibility": float(data_arrays['vis_sums'][i]),
                                    "data": float(data_arrays['storage'][i]),
                                    "data_payload": float(data_arrays['storage_payload'][i]),
                                    "data_HK": float(data_arrays['storage_hk'][i]),
                                    "battery": float(data_arrays['battery'][i]),
                                    "consumption": float(data_arrays['consumption'][i]),
                                    "generation": float(data_arrays['generation'][i]),
                                    "eclipse": float(data_arrays['eclipse'][i]),
                                    "modes": float(data_arrays['modes'][i]),
                                    "altitude": float(data_arrays['altitudes'][i]),
                                    "RAAN": float(data_arrays['RAANs'][i]),
                                    "AOP": float(data_arrays['AOPs'][i]),
                                    "ECC": float(data_arrays['ECCs'][i]),
                                    "INC": float(data_arrays['INCs'][i]),
                                    "density": float(data_arrays['density'][i]),
                                    "Lat": float(data_arrays['latitude'][i]),
                                    "Lng": float(data_arrays['longitude'][i]),
                                    "solar_cells_efficiency": float(data_arrays['solar_eff'][i])
                                },
                                "time": timestamps[i]
                            }
                            for i in range(n_new_points)
                        ]
                        
                        # Send batch data - InfluxDB client optimizes this internally
                        self.write_api.write(bucket=self.influxdbbucket, org=self.influxdborg, record=batch_data)
                        self.influxdb_last_upload_t = t

        end_for_loop = time.time()
        duration = end_for_loop - start_for_loop
        print("Simulation ended!") if self.verbose else None
        print(f"Duration: {duration} s") if self.verbose else None

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
            solar_cells_efficiency[:last_ind]
            data_storage = data_storage[:last_ind]
            data_storage_payload = data_storage_payload[:last_ind]
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
            "storage_payload": data_storage_payload,
            "storage_HK": data_storage_HK,
            "duration_sim": self.duration_sim,
            "epochs_array": self.epochs_array,
            "initial_orbit": self.propagator.get_initial_orbit(),
            "ground_stations": self.ground_stations,
            "orbit_state": orbit_state,
            "spacecraft_state": spacecraft_state,
            "density_array": density_array,
            "solar_cells_efficiency": solar_cells_efficiency
        }

        # Produce report
        produce_report(data_results, self.report_params, results_folder, self.verbose)
        print("Results saved!") if self.verbose else None

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

    def _process_commands(self) -> None:
        """Process all pending commands in the queue."""
        
        while not self.command_queue.empty():
            try:
                command = self.command_queue.get_nowait()
                success, message = self.command_processor.execute_command(command)
                if self.verbose:
                    print(f"Command {command['command']}: {'✓' if success else '✗'} {message}")
            except Exception as e:
                print(f"Error processing command: {e}")

    def send_command(self, command: str, params: dict = {}) -> None:
        """Send a command to the simulation.
        
        Args:
            command: Command name (e.g., 'set_mode', 'uplink_safe_mode')
            params: Command parameters
        """
        command_dict = {"command": command, "params": params}
        self.command_queue.put(command_dict)
