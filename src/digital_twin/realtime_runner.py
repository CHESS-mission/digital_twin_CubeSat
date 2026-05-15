"""Controlled real-time stepping support for the Digital Twin.

This module intentionally does not modify the existing batch ``Simulation.run``
path. It mirrors the batch loop one step at a time so an external bridge can
observe snapshots and inject commands between steps.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any

import numpy as np
from astropy import units as u

from digital_twin import Simulation
from digital_twin.utils import (
    convert_cartesian_to_spherical,
    parse_data_file,
)
from digital_twin.constants import earth_R


SIMULATION_FOLDER = "data/simulation"
ORBIT_FOLDER = "data/orbit"
SPACECRAFT_FOLDER = "data/spacecraft"
GROUND_STATION_FOLDER = "data/ground_station"
MISSION_DESIGN_FOLDER = "data/mission_design"
ENV_FILE = ".env"


@dataclass(frozen=True)
class SimulationFiles:
    """Input files used to create a Digital Twin run."""

    simulation: str = "simulation_template.json"
    orbit: str = "orbit_template.json"
    spacecraft: str = "spacecraft_template.json"
    ground_station: str = "ground_station_template.json"
    mission_design: str = "mission_design_template.json"


@dataclass(frozen=True)
class SimulationSnapshot:
    """A normalized bridge-facing Digital Twin state snapshot."""

    step: int
    simulation_time_s: float
    epoch: str
    mode: int
    eclipse: bool
    visibility: list[int]
    visible: bool
    latitude_deg: float
    longitude_deg: float
    altitude_km: float
    position_km: list[float]
    velocity_km_s: list[float]
    battery_ws: float
    battery_wh: float
    power_consumption_w: float
    power_generation_w: float
    data_storage_mbit: float
    data_storage_payload_mbit: float
    data_storage_hk_mbit: float
    density: float
    solar_cells_efficiency: float

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def build_simulation(
    project_root: Path,
    files: SimulationFiles | None = None,
    *,
    simulation_overrides: dict[str, Any] | None = None,
    disable_influxdb: bool = True,
    quiet: bool = True,
) -> Simulation:
    """Build a ``Simulation`` from standard Digital Twin config files."""

    files = files or SimulationFiles()

    simulation_params = parse_data_file(
        str(project_root / SIMULATION_FOLDER / files.simulation)
    )
    orbit_params = parse_data_file(str(project_root / ORBIT_FOLDER / files.orbit))
    spacecraft_params = parse_data_file(
        str(project_root / SPACECRAFT_FOLDER / files.spacecraft)
    )
    station_params = parse_data_file(
        str(project_root / GROUND_STATION_FOLDER / files.ground_station)
    )
    mission_design_params = parse_data_file(
        str(project_root / MISSION_DESIGN_FOLDER / files.mission_design)
    )

    if simulation_overrides:
        simulation_params.update(simulation_overrides)

    # The wrapper owns pacing. The batch runner's real-time sleep and InfluxDB
    # upload path are not used for bridge-facing stepping.
    simulation_params["run_real_time"] = False
    if disable_influxdb:
        simulation_params["influxdb_delta_t"] = -1
    if quiet:
        simulation_params["verbose"] = "no"
        simulation_params["print_initial_parameters"] = "no"

    return Simulation(
        simulation_params,
        orbit_params,
        spacecraft_params,
        station_params,
        mission_design_params,
        str(project_root / ENV_FILE),
    )


class RealtimeSimulationRunner:
    """Step a Digital Twin simulation and expose bridge-facing snapshots."""

    def __init__(self, simulation: Simulation) -> None:
        self.simulation = simulation
        self.current_step = 0
        self.status = "ready"
        self.error: str | None = None

        n = self.simulation.n_timesteps + 1
        self.eph = np.zeros((n, 6))
        self.eph[0, :3] = self.simulation.propagator.r
        self.eph[0, 3:] = self.simulation.propagator.v

        self.modes = np.zeros(n)
        self.modes[0] = self.simulation.switch_algo.operating_mode

        self.battery_energies = np.zeros(n)
        self.battery_energies[0] = (
            self.simulation.spacecraft.get_eps().get_battery_energy().value
        )

        self.power_consumption = np.zeros(n)
        self.power_consumption[0] = (
            self.simulation.spacecraft.get_eps().get_power_consumption().value
        )

        self.power_generation = np.zeros(n)
        self.power_generation[0] = (
            self.simulation.spacecraft.get_eps().get_power_generation().value
        )

        self.solar_cells_efficiency = np.zeros(n)
        self.solar_cells_efficiency[0] = (
            self.simulation.spacecraft.get_eps().get_solar_cells_efficiency().value
        )

        self.data_storage = np.zeros(n)
        self.data_storage_payload = np.zeros(n)
        self.data_storage_hk = np.zeros(n)
        stored_all, stored_payload, stored_hk = (
            self.simulation.spacecraft.get_obc().get_data()
        )
        self.data_storage[0] = stored_all.value
        self.data_storage_payload[0] = stored_payload.value
        self.data_storage_hk[0] = stored_hk.value

        self.vis_windows = np.zeros((n, len(self.simulation.ground_stations)))
        visibility, _ = self.simulation.propagator.calculate_vis_window(
            self.simulation.ground_stations
        )
        self.vis_windows[0] = visibility.astype(int)

        self.eclipse_windows = np.zeros(n)
        eclipse, _ = self.simulation.propagator.calculate_eclipse_status()
        self.eclipse_windows[0] = int(eclipse)

        self.density_array = np.zeros(n)
        self.density_array[0] = self.simulation.propagator.get_density().value

        self.latest_snapshot = self._snapshot_for_step(0)

    @property
    def delta_t_s(self) -> float:
        return float(self.simulation.delta_t.to_value(u.second))

    @property
    def total_steps(self) -> int:
        return int(self.simulation.n_timesteps)

    def send_command(self, command: str, params: dict[str, Any] | None = None) -> None:
        self.simulation.send_command(command, params or {})

    def step(self) -> SimulationSnapshot:
        """Advance the simulation by one timestep and return the new snapshot."""

        if self.status == "completed":
            return self.latest_snapshot
        if self.current_step >= self.total_steps:
            self.status = "completed"
            return self.latest_snapshot

        step_index = self.current_step
        self.simulation._process_commands()

        try:
            rv = self.simulation.propagator.propagate(
                self.simulation.delta_t,
                self.simulation.spacecraft.C_D,
                self.simulation.spacecraft.A_over_m,
            )
        except (RuntimeError, ValueError, ZeroDivisionError) as exc:
            self.status = "failed"
            self.error = str(exc)
            return self.latest_snapshot

        next_index = step_index + 1
        self.eph[next_index, :3] = rv[:3]
        self.eph[next_index, 3:] = rv[3:]

        if self.simulation.propagation_only:
            self._copy_previous_state(next_index)
        else:
            self._update_spacecraft_state(step_index, next_index, rv)

        self.current_step = next_index
        self.latest_snapshot = self._snapshot_for_step(next_index)
        if self.current_step >= self.total_steps:
            self.status = "completed"
        else:
            self.status = "running"
        return self.latest_snapshot

    def _copy_previous_state(self, next_index: int) -> None:
        prev = next_index - 1
        self.modes[next_index] = self.modes[prev]
        self.vis_windows[next_index] = self.vis_windows[prev]
        self.eclipse_windows[next_index] = self.eclipse_windows[prev]
        self.battery_energies[next_index] = self.battery_energies[prev]
        self.power_consumption[next_index] = self.power_consumption[prev]
        self.power_generation[next_index] = self.power_generation[prev]
        self.solar_cells_efficiency[next_index] = self.solar_cells_efficiency[prev]
        self.data_storage[next_index] = self.data_storage[prev]
        self.data_storage_payload[next_index] = self.data_storage_payload[prev]
        self.data_storage_hk[next_index] = self.data_storage_hk[prev]
        self.density_array[next_index] = self.simulation.propagator.get_density().value

    def _update_spacecraft_state(
        self,
        step_index: int,
        next_index: int,
        rv: np.ndarray,
    ) -> None:
        visibility, gs_coords_array = self.simulation.propagator.calculate_vis_window(
            self.simulation.ground_stations
        )
        eclipse_status, r_earth_sun = (
            self.simulation.propagator.calculate_eclipse_status()
        )

        measurement_session = (
            self.simulation.spacecraft.get_payload().can_start_measuring(
                self.simulation.tofs[next_index].to("second")
            )
        )

        com_window = False
        gs_coords = None
        for idx, visible in enumerate(visibility):
            if visible:
                com_window = True
                gs_coords = (gs_coords_array[idx] * u.km).flatten()
                break

        safe_flag = self.simulation.spacecraft.get_obc().raise_spacecraft_safe_flag()

        old_mode = self.simulation.switch_algo.operating_mode
        self.simulation.switch_algo.switch_mode(
            self.simulation.spacecraft.get_eps(),
            self.simulation.spacecraft.get_telecom(),
            self.simulation.spacecraft.get_payload(),
            self.simulation.spacecraft.get_data_storage(),
            com_window,
            eclipse_status,
            measurement_session,
            safe_flag,
        )
        new_mode = self.simulation.switch_algo.operating_mode

        self.simulation.spacecraft.update_subsystems(
            old_mode,
            new_mode,
            rv,
            com_window,
            eclipse_status,
            self.simulation.delta_t,
            r_earth_sun,
            gs_coords,
            step_index,
        )

        self.modes[next_index] = new_mode
        self.vis_windows[next_index] = visibility.astype(int)
        self.eclipse_windows[next_index] = int(eclipse_status)
        self.battery_energies[next_index] = (
            self.simulation.spacecraft.get_eps().get_battery_energy().value
        )
        self.power_consumption[next_index] = (
            self.simulation.spacecraft.get_eps().get_power_consumption().value
        )
        self.power_generation[next_index] = (
            self.simulation.spacecraft.get_eps().get_power_generation().value
        )
        self.solar_cells_efficiency[next_index] = (
            self.simulation.spacecraft.get_eps().get_solar_cells_efficiency().value
        )

        stored_all, stored_payload, stored_hk = (
            self.simulation.spacecraft.get_obc().get_data()
        )
        self.data_storage[next_index] = stored_all.value
        self.data_storage_payload[next_index] = stored_payload.value
        self.data_storage_hk[next_index] = stored_hk.value
        self.density_array[next_index] = self.simulation.propagator.get_density().value

    def _snapshot_for_step(self, step: int) -> SimulationSnapshot:
        state_vector = self.eph[step]
        position = state_vector[:3]
        velocity = state_vector[3:]

        latitude, longitude = convert_cartesian_to_spherical(
            position.reshape(1, 3),
            self.simulation.epochs_array[step : step + 1],
        )
        altitude_km = float(np.linalg.norm(position) - earth_R.value)

        battery_ws = float(self.battery_energies[step])
        visibility = [int(value) for value in np.asarray(self.vis_windows[step]).flat]

        return SimulationSnapshot(
            step=int(step),
            simulation_time_s=float(self.simulation.tofs[step].to_value("second")),
            epoch=self.simulation.epochs_array[step].isot,
            mode=int(self.modes[step]),
            eclipse=bool(self.eclipse_windows[step]),
            visibility=visibility,
            visible=any(bool(value) for value in visibility),
            latitude_deg=float(latitude[0]),
            longitude_deg=float(longitude[0]),
            altitude_km=altitude_km,
            position_km=[float(value) for value in position],
            velocity_km_s=[float(value) for value in velocity],
            battery_ws=battery_ws,
            battery_wh=battery_ws / 3600.0,
            power_consumption_w=float(self.power_consumption[step]),
            power_generation_w=float(self.power_generation[step]),
            data_storage_mbit=float(self.data_storage[step]),
            data_storage_payload_mbit=float(self.data_storage_payload[step]),
            data_storage_hk_mbit=float(self.data_storage_hk[step]),
            density=float(self.density_array[step]),
            solar_cells_efficiency=float(self.solar_cells_efficiency[step]),
        )
