"""
File for the On-Board Computer (OBC) subsystem implementation.

This module defines the DataStorage and Obc classes used for managing the spacecraft's 
on-board computer operations, including data storage, generation of housekeeping (HK) data, 
and handling safe mode flags.

Classes:
    DataStorage: Represents the data storage system responsible for managing data.
    Obc: Represents the On-Board Computer subsystem.
"""

from astropy import units as u
from astropy.time import TimeDelta
from astropy.units import Quantity
import numpy as np

from digital_twin.spacecraft import SubSystem


class DataStorage:
    """
    Represent the data storage system for Housekeeping and science (GNSS And TOF) data.

    Attributes:
        max_storage (Quantity["data quantity"]): Maximum storage capacity of the system.
        data_storage (Quantity["data quantity"]): Current total data stored.
        data_TOF_GNSS (Quantity["data quantity"]): Current data from TOF and GNSS systems stored.
        data_HK (Quantity["data quantity"]): Cuerrent Housekeeping data stored.
        data_to_downlink (Quantity["data quantity"]): Data queued for downlink during X-band communication.
    """

    def __init__(self, params: dict, verbose: bool = False):
        print("Initializing the data storage... ") if verbose else None
        self.max_storage = float(params["max_storage"]) * u.Mbit

        # Data_storage gathers all 3 types of data (TOF, GNSS, HK = Housekeeping)
        self.data_storage = float(params["init_data"]) * u.Mbit
        self.data_TOF_GNSS = float(params["init_data_TOF_GNSS"]) * u.Mbit
        self.data_HK = float(params["init_data_HK"]) * u.Mbit

        # We need data_TOF_GNSS + data_HK = data_storage
        if abs(self.data_TOF_GNSS + self.data_HK - self.data_storage).to_value() > 1e-6:
            raise (RuntimeError)

        self.data_to_downlink = (
            float(params["init_data_to_downlink"]) * u.Mbit
        )  # Updated at the start of an X-band comm window

    def update_data_storage(
        self,
        data_update_TOF_GNSS: Quantity["data quantity"],
        data_update_HK: Quantity["data quantity"],
    ) -> None:
        """Update the data storage with new TOF/GNSS and HK data."""
        # TOF/GNSS data: cannot remove more than what it has
        if (self.data_TOF_GNSS + data_update_TOF_GNSS) < 0:
            self.data_storage -= self.data_TOF_GNSS
            self.data_to_downlink -= self.data_TOF_GNSS
            self.data_TOF_GNSS = 0.0 * u.Mbit
        else:
            self.data_storage += data_update_TOF_GNSS
            self.data_TOF_GNSS += data_update_TOF_GNSS
            self.data_to_downlink += data_update_TOF_GNSS

        # Housekeeping data
        if (self.data_HK + data_update_HK) < 0:
            self.data_storage -= self.data_HK
            self.data_HK = 0.0 * u.Mbit
        else:
            self.data_storage += data_update_HK
            self.data_HK += data_update_HK

    def x_band_data_empty(self) -> bool:
        """Check if all data to downlink (defined at the start of the XBAND-COM window) has been downlinked."""
        if self.data_to_downlink <= 0:
            return True
        else:
            return False

    def data_storage_empty(self, type: str = "all") -> bool:
        """
        Check if the data storage is empty for a specific type of data.

        Parameters:
            type (str): Type of data to check ('all', 'x_band', or 'uhf').
        """
        if type == "all":
            if self.data_storage <= 0:
                return True
            else:
                return False
        elif type == "x_band":
            if self.data_TOF_GNSS <= 0:
                return True
            else:
                False
        elif type == "uhf":  # housekeeping
            if self.data_HK <= 0:
                return True
            else:
                return False
        else:
            raise NotImplementedError("This type of data does not exist!")

    def data_storage_full(self) -> bool:
        """Check if the data storage is full."""
        if self.data_storage >= self.max_storage:
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
        """Update the class state during the simulation."""
        if new_mode == 4:  # X_BAND
            if new_mode != old_mode:  # Just switched to x_band_comm mode
                # Update how much data needs to be downlinked
                self.data_to_downlink = (
                    self.data_TOF_GNSS.value * u.Mbit
                )  # Not a direct = or it does not do a deepcopy, which introduces a bug

    # Getters
    def get_data_storage(self) -> tuple[Quantity["data quantity"]]:
        return self.data_storage, self.data_TOF_GNSS, self.data_HK

    def get_data_to_downlink(self) -> Quantity["data quantity"]:
        return self.data_to_downlink

    def __str__(self) -> str:
        """Return a string representation of the Data Storage instance."""
        return f"   - Maximum storage: {self.max_storage}"


class Obc(SubSystem):
    """
    Represent the On-Board Computer subsystem.

    Attributes:
        name (str): The name of the subsystem.
        consumption_mean (dict[int, Quantity["power"]]): Mean power consumption for different operating modes.
        HK_rate (Quantity["bandwidth"]): Rate of housekeeping data generation.
        data_storage (DataStorage): Instance of the DataStorage class.
        safe_flag_spacecraft_raised (bool): Indicates if a safe flag was raised by any subsystem.
        safe_flag_spacecraft_subsystem (Optional[str]): Subsystem that raised the safe flag, if any.
        safe_flag (bool): Indicates whether the OBC subsystem is in safe mode.
    """

    def __init__(
        self, params: dict, init_operating_mode: int, verbose: bool = False
    ) -> None:
        print("Initializing OBC subsystem... ") if verbose else None
        self.name = "OBC"

        super(Obc, self).__init__()

        self.consumption_mean = {
            int(k): v * u.W for k, v in params["consumption"].items()
        }

        # Obc generates housekeeping data
        self.HK_rate = float(params["housekeeping_data_rate"]) * (u.Mbit / u.s)

        self.data_storage = DataStorage(params["data_storage"], verbose)

        # Safe flag handling
        # These 2 will be automatically updated at initialization since the Spacecraft calls update_safe_flag()
        self.safe_flag_spacecraft_raised = (
            False  # If a safe flag was raised by ANY subsystem
        )
        self.safe_flag_spacecraft_subsystem = (
            None  # Which subsystem raised the same flag
        )

        # Safe flag regarding current subsystem (OBC)
        self.safe_flag = False if params["init_safe_flag"] == "false" else True

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
        Update the OBC state.

        Parameters:
            old_mode (int): Previous operating mode of the spacecraft.
            new_mode (int): New operating mode of the spacecraft.
            rv (np.ndarray): Position vector of the spacecraft.
            com_window (bool): Indicate if the spacecraft is in a communication window.
            eclipse_status (bool): Indicate if the spacecraft is in eclipse.
            delta_t (TimeDelta): Timestep for the update.
        """
        self.data_storage.update(
            old_mode, new_mode, rv, com_window, eclipse_status, delta_t
        )

        # SAFE FLAG HANDLING
        # Check safe flag triggers (cannot generate a safe flag if already in safe mode)
        if new_mode != 1 and self.safe_flag == False:
            pass  # Not implemented yet for this subsystem
        # Check safe flag resolution
        if self.safe_flag == True:
            pass  # Not implemented yet for this subsystem

    def compute_power_consumed(self, mode: int) -> Quantity["power"]:
        """Compute the power consumed in the specified mode."""
        return self.consumption_mean[mode]

    def update_data_storage(
        self,
        data_update_TOF_GNSS: Quantity["data quantity"],
        data_update_HK: Quantity["data quantity"],
        delta_t: TimeDelta,
    ) -> None:
        """
        Update the data storage with new data and housekeeping data.

        Parameters:
            data_update_TOF_GNSS (Quantity["data quantity"]): Update for TOF/GNSS data.
            data_update_HK (Quantity["data quantity"]): Update for housekeeping data.
            delta_t (TimeDelta): Timestep for the update.
        """
        # Housekeeping (HK) data also added (dooesn't depend on the operating mode, always generated)
        data_update_HK += self.HK_rate * delta_t
        self.data_storage.update_data_storage(data_update_TOF_GNSS, data_update_HK)

    def __str__(self) -> str:
        """Return a string representation of the OBC subsystem."""
        string_data_storage = str(self.data_storage)
        return f"Obc: \n- Housekeeping data generation rate: {self.HK_rate} \n- Data storage:\n{string_data_storage}"

    # This is a method checking safe flags for the OBC subsystem
    def raise_safe_flag(self) -> bool:
        """Return true if a safe flag is raised by this subsystem in order to trigger safe mode."""
        return self.safe_flag

    def update_safe_flag(
        self, safe_flag_raised: bool, safe_flag_subsystem: str
    ) -> None:
        """Check if a safe flag was raised by any other subsystem and update the value returned by raise_spacecraft_safe_flag()"""
        if safe_flag_raised:
            self.safe_flag_spacecraft_raised = True
            self.safe_flag_spacecraft_subsystem = safe_flag_subsystem
        else:
            self.safe_flag_spacecraft_raised = False
            self.safe_flag_spacecraft_subsystem = None

    def raise_spacecraft_safe_flag(self) -> bool:
        """Check if safe flag raised by a subsystem in the spacecraft"""
        return self.safe_flag_spacecraft_raised

    # Getters
    def get_name(self) -> str:
        return self.name

    def get_data(self) -> tuple[Quantity["data quantity"]]:
        return self.data_storage.get_data_storage()

    def get_data_storage(self) -> DataStorage:
        return self.data_storage
