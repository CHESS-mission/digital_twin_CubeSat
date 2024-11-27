from typing import Dict

from astropy.time import TimeDelta
from astropy.units import Quantity
import numpy as np


# TODO: check if it is better to implement with an @abstractclass decorator instead
class SubSystem:
    """Interface for subsystems to make sure all susbsystems have necessary functions implemented."""

    # this init function is call by all subsystems
    def __init__(self):
        pass

    def update(
        self,
        old_mode: str,
        new_mode: str,
        rv: np.ndarray,
        com_window: bool,
        eclipse_status: bool,
        delta_t: TimeDelta,
    ) -> None:
        raise (NotImplementedError)

    def compute_power_consumed(self, mode: int) -> Quantity:
        raise (NotImplementedError)

    def raise_safe_flag(self) -> bool:
        raise (NotImplementedError)

    def get_name(self) -> str:
        raise (NotImplementedError)

    def __str__(self) -> str:
        raise (NotImplementedError)
