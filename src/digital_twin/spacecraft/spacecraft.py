"""File for the spacecraft object.
"""

from typing import Dict, Tuple

from astropy import units as u
from astropy.time import TimeDelta
from astropy.units import Quantity
import numpy as np

from digital_twin.spacecraft.adcs import Adcs
from digital_twin.spacecraft.eps import Eps
from digital_twin.spacecraft.obc import Obc
from digital_twin.spacecraft.payload import Payload
from digital_twin.spacecraft.telecom import Telecom


class Spacecraft:
    def __init__(self, params: Dict, init_operating_mode: int) -> None:
        print("Initializing the spacecraft...")
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

    def update_subsystems(
        self,
        old_mode: str,
        new_mode: str,
        rv: np.ndarray,
        com_window: bool,
        eclipse_status: bool,
        delta_t: TimeDelta,
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
        self.eps_subsystem.update_batteries(
            power_consumed, delta_t, eclipse_status, new_mode
        )

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
