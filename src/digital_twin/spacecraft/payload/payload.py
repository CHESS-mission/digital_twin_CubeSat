"""
File for the Payload subsystem, which handles GNSS and TOF operations.

This module defines the Payload class responsible for managing the spacecraft's 
payload operations, including GNSS and TOF measurements.

Classes:
    Payload: Represents the Payload subsystem managing GNSS and TOF systems.
"""

import math

from astropy import units as u
from astropy.time import TimeDelta
from astropy.units import Quantity
import numpy as np

from digital_twin.spacecraft import SubSystem
from digital_twin.utils import get_astropy_unit_time


class Payload(SubSystem):
    """
    Represent the Payload subsystem, which includes the GNSS and TOF systems for the spacecraft.

    Attributes:
        name (str): The name of the subsystem.
        consumption_mean_gnss (dict[int, Quantity["power"]]): Mean power consumption for the GNSS system at different operating modes.
        consumption_mean_tof (dict[int, Quantity["power"]]): Mean power consumption for the TOF system at different operating modes.
        measurement_TOF_rate (Quantity["bandwidth"]): Rate at which TOF data is generated.
        measurement_GNSS_rate (Quantity["bandwidth"]): Rate at which GNSS data is generated.
        measurement_max_duration (Quantity["time"]): Maximum duration for a single measurement session.
        nb_measurement_per_day (int): Number of measurements allowed per day.
        start_measurement (Quantity["time"]): Duration of the pre-conditioning period before measurements can begin.
        stop_measurement (Quantity["time"]): Duration after which the measurement session stops due to post-conditioning.
        measurement_duration (Quantity["time"]): Duration of the current measurement session.
        is_measuring (bool): Indicates whether the subsystem is currently measuring.
        nb_measurement_windows (int): Counter for the number of measurement windows during the simulation.
        safe_flag (bool): Indicates whether the Payload subsystem is in safe mode.
    """

    def __init__(self, params: dict, init_operating_mode: int) -> None:
        print("Initializing Payload subsystem... ")
        self.name = "Payload"

        super(Payload, self).__init__()

        self.consumption_mean_gnss = {
            int(k): v * u.W for k, v in params["consumption_gnss"].items()
        }
        self.consumption_mean_tof = {
            int(k): v * u.W for k, v in params["consumption_tof"].items()
        }

        self.measurement_TOF_rate = float(params["measurement_TOF_rate"]) * (
            u.Mbit / u.s
        )
        self.measurement_GNSS_rate = float(params["measurement_GNSS_rate"]) * (
            u.Mbit / u.s
        )
        measure_unit = get_astropy_unit_time(params["max_duration_unit"])
        self.measurement_max_duration = (params["max_duration"] * measure_unit).to(u.s)
        self.nb_measurement_per_day = int(params["nb_measurements_per_day"])

        # Duration of the pre conditioning before the TOF can actually generate data (at the beginning of the measurement session)
        self.start_measurement = (
            float(params["measurement_pre_conditioning_time"]) * u.s
        )
        # Duration of post conditioning at the end of a measurement session
        measurement_post_conditioning = (
            float(params["measurement_post_conditioning_time"]) * u.s
        )
        self.stop_measurement = (
            self.measurement_max_duration - measurement_post_conditioning
        )  # Stop the measurement earlier than measurement_max_duration to have enough time to do the post conditioning before the measurement session stops

        self.measurement_duration = float(params["init_measurement_duration"]) * u.s
        self.is_measuring = False if params["init_is_measuring"] == "false" else True
        self.nb_measurement_windows = (
            0  # Count the number of measurement windows during the simulation
        )

        self.safe_flag = False if params["init_safe_flag"] == "false" else True

    def can_start_measuring(self, time_elapsed: Quantity["time"]) -> bool:
        """
        Determine if a new measurement session can start based on elapsed time.

        Parameters:
            time_elapsed (Quantity["time"]): The total elapsed time since the start of the mission.
        """
        # the max number of measurement sessions is defined to be PER DAY currently
        one_day = 86400 * u.s
        if (
            math.floor(self.nb_measurement_windows / self.nb_measurement_per_day)
            < (time_elapsed / one_day).value
        ):
            return True
        else:
            return False

    # This function is used by the switch tree algorithm to stop the measurement window if the max time was reached
    def campaign_finished(self) -> bool:
        """Check if the current measurement campaign is finished based on duration (measurement_max_duration)."""
        if not self.is_measuring:
            return True
        else:
            if self.measurement_duration >= self.measurement_max_duration:
                return True
            return False

    def update(
        self,
        old_mode: int,
        new_mode: int,
        rv: np.ndarray,
        com_window: bool,
        eclipse_status: bool,
        delta_t: TimeDelta,
    ) -> None:
        """
        Update the Payload subsystem state.

        Parameters:
            old_mode (int): Previous operating mode of the subsystem.
            new_mode (int): New operating mode of the subsystem.
            rv (np.ndarray): Position vector of the spacecraft.
            com_window (bool): Indicate if the spacecraft is in a communication window.
            eclipse_status (bool): Indicate if the spacecraft is in eclipse.
            delta_t (TimeDelta): Timestep for the update.
        """
        if new_mode == 5:  # MEASUREMENT
            if old_mode != new_mode:  # Just switched to measurement mode
                self.measurement_duration = 0.0 * u.s
                self.nb_measurement_windows += 1
            self.measurement_duration += delta_t
            self.is_measuring = True

        else:
            self.measurement_duration = 0.0 * u.s
            self.is_measuring = False

        # SAFE FLAG HANDLING
        # Check safe flag triggers (cannot generate a safe flag if already in safe mode)
        if new_mode != 1 and self.safe_flag == False:
            pass  # Not implemented yet for this subsystem
        # Check safe flag resolution
        if self.safe_flag == True:
            pass  # Not implemented yet for this subsystem

    def compute_power_consumed(self, mode: int) -> Quantity["power"]:
        """Compute the power consumed in the specified mode."""
        return self.consumption_mean_gnss[mode] + self.consumption_mean_tof[mode]

    def __str__(self) -> str:
        """Return a string representation of the Payload subsystem."""
        return f"Payload:"

    def compute_data_update(
        self, new_mode: int, delta_t: TimeDelta
    ) -> tuple[Quantity["data quantity"]]:
        """Calculate how much data was generated by the payloads during the current timestep.

        Returns:
            tuple which first term is generated scientific data and second term is generated housekeeping data (0 in this case as it is not generated by this subsystem)
        """
        data = 0.0 * u.Mbit
        if (
            new_mode == 5
            and self.measurement_duration >= self.start_measurement
            and self.measurement_duration < self.stop_measurement
        ):  # Only add data if in measurement mode and after/before pre/post conditioning (power budget)
            data += self.measurement_TOF_rate * delta_t

        if new_mode != 1:  # If not in safe mode, GNSS continuously add data
            data += self.measurement_GNSS_rate * delta_t

        return (
            data,
            0.0 * u.Mbit,
        )  # Second term for HouseKeeping data which is 0 in this case

    def raise_safe_flag(self) -> bool:
        """Return true if a safe flag is raised by this subsystem in order to trigger safe mode."""
        return self.safe_flag

    # Getters
    def get_name(self) -> str:
        return self.name

    def get_measurement_variables(self) -> tuple[Quantity["time"], bool]:
        return self.measurement_duration, self.is_measuring
