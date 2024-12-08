"""File for the spacecraft object.
"""

from typing import Dict, Tuple

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
    def __init__(self, params: Dict, init_operating_mode: int) -> None:
        print("Initializing the spacecraft...")
        self.params = params
        self.name = params["general"]["name"]
        self.C_D = params["general"]["C_D"] * u.one
        init_cross_section = params["general"]["cross_section_area"] * u.km**2
        self.cross_section = init_cross_section  # NOTE: will have to update this when I implement ADCS module!
        self.mass = params["general"]["mass"] * u.kg
        self.A_over_m = ((self.cross_section) / (self.mass)).to_value(
            u.km**2 / u.kg
        ) * (u.km**2 / u.kg)
        self.B = self.C_D * self.A_over_m  # Approximation of the ballistic coefficient

        # SUBSYSTEMS INITIALIZATION (main ones)
        self.eps_subsystem = Eps(params["eps"], init_operating_mode)
        self.telecom_subsystem = Telecom(params["telecom"], init_operating_mode)
        self.adcs_subsystem = Adcs(params["adcs"], init_operating_mode)
        self.payload_subsystem = Payload(params["payload"], init_operating_mode)
        self.obc_subsystem = Obc(params["obc"], init_operating_mode)

        self.subsystems = [
            self.eps_subsystem,
            self.telecom_subsystem,
            self.adcs_subsystem,
            self.payload_subsystem,
            self.obc_subsystem,
        ]

        # check if any safe flag is raised by subsystem and update OBC
        safe_flag_raised = False
        safe_flag_subsystem = None
        for subsystem in self.subsystems:
            if subsystem.raise_safe_flag():
                safe_flag_raised = True
                safe_flag_subsystem = subsystem.get_name()
        self.obc_subsystem.update_safe_flag(safe_flag_raised, safe_flag_subsystem)

    def update_subsystems(
        self,
        old_mode: str,
        new_mode: str,
        rv: np.ndarray[Quantity],
        com_window: bool,
        eclipse_status: bool,
        delta_t: TimeDelta,
        r_earth_sun: Quantity,
        gs_coords: Quantity,
    ) -> None:
        # Update each subsystem based on old mode, new mode, location, communication_window and eclispe_status boolean
        for subsystem in self.subsystems:
            subsystem.update(
                old_mode, new_mode, rv, com_window, eclipse_status, delta_t
            )

        # get new cross section from ADCS change of orientation, and update related spacecraft attributes
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
            rv[:3] * u.km,
            gs_coords,
        )

        # Compute data change (generation or removal) at current timestep
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

        # check if any safe flag is raised by subsystem and update OBC
        safe_flag_raised = False
        safe_flag_subsystem = None
        for subsystem in self.subsystems:
            if subsystem.raise_safe_flag():
                safe_flag_raised = True
                safe_flag_subsystem = subsystem.get_name()
        self.obc_subsystem.update_safe_flag(safe_flag_raised, safe_flag_subsystem)

    def __str__(self):
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

    def save_state(self) -> Dict:
        spacecraft_state = self.params  # the initial one, which is now updated

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

        full, TOF_GNSS, HK = self.obc_subsystem.get_data()
        spacecraft_state["obc"]["data_storage"]["init_data"] = full.to_value()
        spacecraft_state["obc"]["data_storage"][
            "init_data_TOF_GNSS"
        ] = TOF_GNSS.to_value()
        spacecraft_state["obc"]["data_storage"]["init_data_HK"] = HK.to_value()
        spacecraft_state["obc"]["data_storage"]["init_data_to_downlink"] = (
            self.obc_subsystem.get_data_storage().get_data_to_downlink().to_value()
        )

        dur, is_measuring = self.payload_subsystem.get_measurement_variables()
        spacecraft_state["payload"]["init_measurement_duration"] = dur.to_value()
        spacecraft_state["payload"]["init_is_measuring"] = (
            "true" if is_measuring else "false"
        )

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
