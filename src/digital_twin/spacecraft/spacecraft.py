"""File for the spacecraft object."""

from typing import Optional

from astropy import units as u
from astropy.time import TimeDelta
from astropy.units import Quantity
import numpy as np

from digital_twin.spacecraft.adcs import Adcs
from digital_twin.spacecraft.eps import Eps
from digital_twin.spacecraft.obc import Obc, DataStorage
from digital_twin.spacecraft.payload import Payload
from digital_twin.spacecraft.telecom import Telecom


class Spacecraft:
    """Represent a spacecraft with various subsystems.

    The Spacecraft class models the behavior and attributes of a spacecraft,
    including its general properties and interactions among its subsystems.

    Attributes:
        params (dict): dictionary of spacecraft parameters given by the user.
        name (str): Name of the spacecraft.
        C_D (Quantity["no_unit"]): Drag coefficient of the spacecraft.
        cross_section (Quantity["area"]): Cross-sectional area of the spacecraft (km²).
        mass (Quantity["mass"]): Mass of the spacecraft (kg).
        A_over_m (Quantity): Area-to-mass ratio of the spacecraft (km²/kg).
        B (Quantity): Ballistic coefficient (C_D * A_over_m).
        eps_subsystem (Eps): Electrical Power System (EPS) subsystem.
        telecom_subsystem (Telecom): Telecommunications subsystem.
        adcs_subsystem (Adcs): Attitude Determination and Control System (ADCS).
        payload_subsystem (Payload): Payload subsystem.
        obc_subsystem (Obc): On-Board Computer (OBC) subsystem.
        subsystems (list): List of all initialized subsystems.
    """

    def __init__(
        self, params: dict, init_operating_mode: int, verbose: bool = False
    ) -> None:
        print("Initializing the spacecraft...") if verbose else None
        self.params = params
        self.name = params["general"]["name"]
        self.C_D = params["general"]["C_D"] * u.one
        init_cross_section = params["general"]["cross_section_area"] * u.km**2
        self.cross_section = (
            init_cross_section  # NOTE: Can be updated with dynamic cross section
        )
        self.mass = params["general"]["mass"] * u.kg
        self.A_over_m = ((self.cross_section) / (self.mass)).to_value(
            u.km**2 / u.kg
        ) * (u.km**2 / u.kg)
        self.B = self.C_D * self.A_over_m  # Approximation of the ballistic coefficient

        # SUBSYSTEMS INITIALIZATION
        self.eps_subsystem = Eps(params["eps"], init_operating_mode, verbose)
        self.telecom_subsystem = Telecom(
            params["telecom"], init_operating_mode, verbose
        )
        self.adcs_subsystem = Adcs(params["adcs"], init_operating_mode, verbose)
        self.payload_subsystem = Payload(
            params["payload"], init_operating_mode, verbose
        )
        self.obc_subsystem = Obc(params["obc"], init_operating_mode, verbose)

        self.subsystems = [
            self.eps_subsystem,
            self.telecom_subsystem,
            self.adcs_subsystem,
            self.payload_subsystem,
            self.obc_subsystem,
        ]

        # Check if any safe flag is raised by subsystem and update OBC
        safe_flag_raised = False
        safe_flag_subsystem = None
        for subsystem in self.subsystems:
            if subsystem.raise_safe_flag():
                safe_flag_raised = True
                safe_flag_subsystem = subsystem.get_name()
        self.obc_subsystem.update_safe_flag(safe_flag_raised, safe_flag_subsystem)

    def update_subsystems(
        self,
        old_mode: int,
        new_mode: int,
        rv: np.ndarray,
        com_window: bool,
        eclipse_status: bool,
        delta_t: TimeDelta,
        r_earth_sun: Quantity["length"],
        gs_coords: Optional[np.ndarray],
    ) -> None:
        """Update all subsystems and spacecraft properties based on current conditions.

        Args:
            old_mode (int): Previous operational mode of the spacecraft.
            new_mode (int): New operational mode of the spacecraft.
            rv (np.ndarray): Position and velocity vector of the spacecraft (km).
            com_window (bool): Indicates whether the spacecraft is in a communication window.
            eclipse_status (bool): Indicates whether the spacecraft is in Earth's shadow.
            delta_t (TimeDelta): Time step for the update.
            r_earth_sun (Quantity["length"]): Vector from Earth to the Sun.
            gs_coords (Optional[np.ndarray]): Coordinates of the ground station if satellite is in visibility window.

        Updates:
            - Individual updates for each subsystem.
            - EPS battery level.
            - Data storage.
            - Safe flags.
        """
        # Update each subsystem independently
        for subsystem in self.subsystems:
            subsystem.update(
                old_mode, new_mode, rv, com_window, eclipse_status, delta_t
            )

        # Get new cross section from ADCS and update related spacecraft attributes
        # Currently the cross section does not change => will change if dynamic cross section is implemented
        new_cross_section = self.adcs_subsystem.get_cross_section(self.cross_section)
        if new_cross_section != self.cross_section:
            self.cross_section = new_cross_section
            self.A_over_m = ((self.cross_section) / (self.mass)).to_value(
                u.km**2 / u.kg
            ) * (u.km**2 / u.kg)
            self.B = self.C_D * self.A_over_m

        # Compute the power consumed by each subsystem at current timestep
        power_consumed = 0
        for subsystem in self.subsystems:
            power_consumed += subsystem.compute_power_consumed(new_mode)
        # Update EPS based on data gathered for all other subsystems
        attitude = self.adcs_subsystem.get_attitude()
        self.eps_subsystem.update_batteries(
            power_consumed,
            delta_t,
            eclipse_status,
            attitude,
            r_earth_sun,
            rv[:3] * u.km,  # position,
            gs_coords,
            new_mode,
        )

        # Compute data change (generation by payload or removal by telecom) at current timestep
        data_update_TOF_GNSS = 0 * u.Mbit
        data_update_HK = 0.0 * u.Mbit
        TOF_GNSS_telecom, HK_telecom = self.telecom_subsystem.compute_data_update(
            new_mode, delta_t
        )
        data_update_TOF_GNSS += TOF_GNSS_telecom
        data_update_HK += HK_telecom

        TOF_GNSS_payload, HK_payload = self.payload_subsystem.compute_data_update(
            new_mode, delta_t
        )
        data_update_TOF_GNSS += TOF_GNSS_payload
        data_update_HK += HK_payload
        self.obc_subsystem.update_data_storage(
            data_update_TOF_GNSS, data_update_HK, delta_t
        )

        # Check if any safe flag is raised by subsystem and update OBC
        safe_flag_raised = False
        safe_flag_subsystem = None
        for subsystem in self.subsystems:
            if subsystem.raise_safe_flag():
                safe_flag_raised = True
                safe_flag_subsystem = subsystem.get_name()
        self.obc_subsystem.update_safe_flag(safe_flag_raised, safe_flag_subsystem)

    def __str__(self) -> str:
        """Return a string representation of the spacecraft and its subsystems."""
        string1 = "\n".join(
            [
                f"- mass: {self.mass}",
                f"- ballistic coeff: {self.B}",
                f"- cross section: {self.cross_section}",
            ]
        )
        tmp = []
        for subsystem in self.subsystems:
            tmp.append(str(subsystem))
        string2 = "\n".join(tmp)
        return f'Spacecraft "{self.name}":\n{string1} \nWith subsystems:\n{string2}'

    # Getters
    @property
    def C_D(self) -> Quantity:
        return self._C_D

    @C_D.setter
    def C_D(self, C_D):
        self._C_D = C_D

    @property
    def A_over_m(self) -> Quantity:
        return self._A_over_m

    @A_over_m.setter
    def A_over_m(self, A_over_m):
        self._A_over_m = A_over_m

    def get_eps(self) -> Eps:
        return self.eps_subsystem

    def get_telecom(self) -> Telecom:
        return self.telecom_subsystem

    def get_payload(self) -> Payload:
        return self.payload_subsystem

    def get_obc(self) -> Obc:
        return self.obc_subsystem

    def get_data_storage(self) -> DataStorage:
        return self.obc_subsystem.get_data_storage()

    def save_state(self) -> dict:
        """Save the current state of the spacecraft and its subsystems. File has the same structure as the input file for the spaceraft parameters.

        Returns:
            dict: Updated dictionary representing the spacecraft's state, including subsystem flags,
                  battery levels, and data storage values. Can be used to start another simulation from this final state.
        """

        spacecraft_state = self.params  # The initial one, which is now updated

        # Update safe flags
        spacecraft_state["eps"]["init_safe_flag"] = (
            "true" if self.eps_subsystem.raise_safe_flag() else "false"
        )
        spacecraft_state["adcs"]["init_safe_flag"] = (
            "true" if self.adcs_subsystem.raise_safe_flag() else "false"
        )
        spacecraft_state["obc"]["init_safe_flag"] = (
            "true" if self.obc_subsystem.raise_safe_flag() else "false"
        )
        spacecraft_state["payload"]["init_safe_flag"] = (
            "true" if self.payload_subsystem.raise_safe_flag() else "false"
        )
        spacecraft_state["telecom"]["init_safe_flag"] = (
            "true" if self.telecom_subsystem.raise_safe_flag() else "false"
        )

        spacecraft_state["eps"]["init_battery_level"] = (
            self.eps_subsystem.get_battery_energy().to(u.W * u.hour).to_value()
        )

        # Update data storage variables
        full, TOF_GNSS, HK = self.obc_subsystem.get_data()
        spacecraft_state["obc"]["data_storage"]["init_data"] = full.to_value()
        spacecraft_state["obc"]["data_storage"][
            "init_data_TOF_GNSS"
        ] = TOF_GNSS.to_value()
        spacecraft_state["obc"]["data_storage"]["init_data_HK"] = HK.to_value()
        spacecraft_state["obc"]["data_storage"]["init_data_to_downlink"] = (
            self.obc_subsystem.get_data_storage().get_data_to_downlink().to_value()
        )

        # Update measurement variables (payload class)
        dur, is_measuring = self.payload_subsystem.get_measurement_variables()
        spacecraft_state["payload"]["init_measurement_duration"] = dur.to_value()
        spacecraft_state["payload"]["init_is_measuring"] = (
            "true" if is_measuring else "false"
        )

        # Update telecom-related variables
        (
            t_no_com,
            com_duration,
            is_communicating,
            alternating,
            safe_flag_reason,
            safe_flag_duration_left,
            is_visible,
        ) = self.telecom_subsystem.get_telecom_variables()

        spacecraft_state["telecom"]["init_time_no_com"] = t_no_com.to_value()
        spacecraft_state["telecom"]["init_com_duration"] = com_duration.to_value()
        spacecraft_state["telecom"]["init_is_communicating"] = (
            "true" if is_communicating else "false"
        )
        spacecraft_state["telecom"]["init_is_alternating"] = (
            "true" if alternating else "false"
        )
        spacecraft_state["telecom"]["init_safe_flag_reason"] = safe_flag_reason
        spacecraft_state["telecom"]["init_safe_flag_duration_left"] = (
            "none"
            if safe_flag_duration_left is None
            else safe_flag_duration_left.to_value()
        )
        spacecraft_state["telecom"]["init_is_visible"] = (
            "true" if is_visible else "false"
        )
        return spacecraft_state
