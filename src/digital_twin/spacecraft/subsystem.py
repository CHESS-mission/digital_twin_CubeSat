"""
Subsystem Interface Module

This module defines the `SubSystem` class, which serves as an abstract interface 
for all spacecraft subsystems.
"""

from astropy.time import TimeDelta
from astropy.units import Quantity
import numpy as np


class SubSystem:
    """Abstract interface for spacecraft subsystems.

    This class ensures that all subsystems implement the necessary methods to
    interact with the spacecraft. Each subsystem is expected to define its
    behavior for updates, power consumption, safety flags, and identification.

    Methods defined here must be implemented by all derived classes.
    """

    # This init function is called by all subsystems
    def __init__(self) -> None:
        pass

    def update(
        self,
        old_mode: int,
        new_mode: int,
        rv: np.ndarray,
        com_window: bool,
        eclipse_status: bool,
        delta_t: TimeDelta,
    ) -> None:
        """Update the subsystem's state based on current spacecraft conditions."""
        raise (NotImplementedError)

    def compute_power_consumed(self, mode: int) -> Quantity["power"]:
        """Calculate the power consumed by the subsystem in the given mode."""
        raise (NotImplementedError)

    def raise_safe_flag(self) -> bool:
        """Determine whether the subsystem raises a safety flag."""
        raise (NotImplementedError)

    def get_name(self) -> str:
        """Get the name of the subsystem."""
        raise (NotImplementedError)

    def __str__(self) -> str:
        """Provide a string representation of the subsystem."""
        raise (NotImplementedError)
