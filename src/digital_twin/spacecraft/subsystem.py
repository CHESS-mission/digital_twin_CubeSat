"""Interface for subsystems to make sure all susbsystems have necessary functions implemented
"""

from typing import Dict

import numpy as np
from astropy.units import Quantity
from astropy.time import TimeDelta


# TODO: check if it is better to implement with an @abstractclass decorator instead
class SubSystem:
    # this init function is call by all subsystems
    def __init__(self):
        pass

    def update(
        self,
        old_mode: str,
        new_mode: str,
        rv: np.array,
        com_window: bool,
        eclipse_status: bool,
        delta_t: TimeDelta,
    ) -> None:
        raise (NotImplementedError)

    def compute_power_consumed(self, mode: int) -> Quantity:
        raise (NotImplementedError)

    def __str__(self) -> str:
        raise (NotImplementedError)
