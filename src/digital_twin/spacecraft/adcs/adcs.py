"""
File for the ADCS (Attitude Determination and Control Subsystem) implementation.

This module defines the ADCS class used for managing the spacecraft's attitude.

Classes:
    Adcs: Represents the Attitude Determination and Control Subsystem (ADCS) for 
          managing spacecraft attitude.
"""

from astropy import units as u
from astropy.time import TimeDelta
from astropy.units import Quantity
import numpy as np

from digital_twin.spacecraft import SubSystem
from digital_twin.constants import attitude_mode_dict, attitude_dict


class Adcs(SubSystem):
    """
    Represent the Attitude Determination and Control Subsystem (ADCS) of the spacecraft.

    Attributes:
        name (str): The name of the subsystem.
        consumption_mean (dict[int, Quantity["power"]]): Mean power consumption for different operating modes.
        safe_flag (bool): Indicates whether the subsystem is in safe mode.
        attitude (int): The current attitude mode of the spacecraft (defined in the file constants.py)
    """

    def __init__(
        self, params: dict, init_operating_mode: int, verbose: bool = False
    ) -> None:
        print("Initializing ADCS subsystem... ") if verbose else None
        self.name = "ADCS"

        super(Adcs, self).__init__()

        self.consumption_mean = {
            int(k): v * u.W for k, v in params["consumption"].items()
        }

        self.safe_flag = False if params["init_safe_flag"] == "false" else True

        self.attitude = attitude_mode_dict[init_operating_mode]

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
        Update the ADCS state, including attitude mode and safe flag handling.

        Parameters:
            old_mode (int): Previous operating mode of the spacecraft.
            new_mode (int): New operating mode of the spacecraft.
            rv (np.ndarray): Position vector of the spacecraft.
            com_window (bool): Indicate if the spacecraft is in a communication window.
            eclipse_status (bool): Indicate if the spacecraft is in eclipse.
            delta_t (TimeDelta): Timestep for the update.
        """
        # Update the attitude
        self.attitude = attitude_mode_dict[new_mode]

        # SAFE FLAG HANDLING
        # Check safe flag triggers (cannot generate a safe flag if already in safe mode)
        if new_mode != 1 and self.safe_flag == False:
            pass  # Not implemented yet for this subsystem
        # Check safe flag resolution
        if self.safe_flag == True:
            pass  # Not implemented yet for this subsystem

    def get_cross_section(
        self, old_cross_section: Quantity["area"]
    ) -> Quantity["area"]:
        """Retrieve the cross-sectional area of the spacecraft in its current attitude. Currently using a constant cross section."""
        return old_cross_section

    def compute_power_consumed(self, mode: int) -> Quantity["power"]:
        """Compute the power consumed in the specified mode."""
        return self.consumption_mean[mode]

    def __str__(self) -> str:
        """Return a string representation of the ADCS subsystem."""
        string1 = "; ".join(list(attitude_dict.values()))
        return f"ADCS: \n- List of attitudes: {string1}"

    def raise_safe_flag(self) -> bool:
        """Return true if a safe flag is raised by this subsystem in order to trigger safe mode."""
        return self.safe_flag

    # Getters
    def get_name(self) -> str:
        return self.name

    def get_attitude(self) -> int:
        return self.attitude
