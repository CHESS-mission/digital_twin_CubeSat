"""Defininition of the logic for switching between different operating modes."""

from digital_twin.constants import mode_dict
from digital_twin.spacecraft.eps import Eps
from digital_twin.spacecraft.obc import DataStorage
from digital_twin.spacecraft.payload import Payload
from digital_twin.spacecraft.telecom import Telecom


class ModeSwitch:
    """Implement the mode switch algorithm defined in the MOP for the CHESS mission."""

    def __init__(self, init_mode: int, verbose: bool = False) -> None:
        print("Initializing the mode switch algorithm...") if verbose else None
        self.mode_dict = mode_dict
        self.operating_mode = init_mode

    def print_operating_mode(self) -> str:
        """Return the human-readable description of the current operating mode."""
        return self.mode_dict[self.operating_mode]

    def switch_mode(
        self,
        eps: Eps,
        telecom: Telecom,
        payload: Payload,
        data_storage: DataStorage,
        com_window: bool,
        eclipse_status: bool,
        measurement_session: bool,
        safe_flag=False,
    ) -> None:
        """
        Determine and set the next operating mode based on the current spacecraft state.

        Args:
            eps (Eps): The Electrical Power Subsystem (EPS) instance.
            telecom (Telecom): The Telecommunication subsystem instance.
            payload (Payload): The Payload subsystem instance.
            data_storage (DataStorage): The onboard data storage instance.
            com_window (bool): Indicates whether a communication window is available.
            eclipse_status (bool): Indicates whether the spacecraft is currently in an eclipse.
            measurement_session (bool): Indicates whether a measurement session can be started.
            safe_flag (bool): Indicates whether a safe flag has been raised.
        """
        if safe_flag:
            self.operating_mode = 1
        else:
            match self.operating_mode:
                case 0:
                    self.mode_switch_from_IDLE(
                        eps,
                        data_storage,
                        com_window,
                        eclipse_status,
                        measurement_session,
                    )
                case 1:
                    self.mode_switch_from_SAFE()
                case 2:
                    self.mode_switch_from_CHARGING(
                        eps,
                        data_storage,
                        com_window,
                        eclipse_status,
                        measurement_session,
                    )
                case 3:
                    self.mode_switch_from_UHF_COM(
                        eps, telecom, data_storage, com_window
                    )
                case 4:
                    self.mode_switch_from_X_BAND_COM(data_storage, com_window)
                case 5:
                    self.mode_switch_from_MEASUREMENT(payload, data_storage)
                case _:
                    raise AssertionError(
                        f"The current operating mode ({self.operating_mode}) does not exist"
                    )

    @property
    def operating_mode(self) -> int:
        return self._operating_mode

    @operating_mode.setter
    def operating_mode(self, mode: int) -> None:
        self._operating_mode = mode

    def mode_switch_from_IDLE(
        self,
        eps: Eps,
        data_storage: DataStorage,
        com_window: bool,
        eclipse_status: bool,
        measurement_session: bool,
    ) -> None:
        """Handle transitions from IDLE mode."""
        # Charging in priority if low battery (if possible)
        if not (eclipse_status) and (eps.battery_level <= eps.min_battery):
            self.operating_mode = 2
        else:
            # Else try to measure
            if (
                (eps.battery_level >= eps.measure_threshold)
                and (not data_storage.data_storage_full())
                and (measurement_session)
            ):
                self.operating_mode = 5
            else:
                # Else try to communicate
                if (eps.battery_level >= eps.com_threshold) and (com_window):
                    self.operating_mode = 3
                else:
                    # Else try to charge (if possible)
                    if (not eclipse_status) and (eps.battery_level < eps.max_battery):
                        self.operating_mode = 2
                    else:
                        self.operating_mode = 0

    def mode_switch_from_SAFE(self) -> None:
        """Handle transitions from SAFE mode."""
        self.operating_mode = 0

    def mode_switch_from_CHARGING(
        self,
        eps: Eps,
        data_storage: DataStorage,
        com_window: bool,
        eclipse_status: bool,
        measurement_session: bool,
    ) -> None:
        """Handle transitions from CHARGING mode."""
        # Try to measure
        if (
            (eps.battery_level >= eps.measure_threshold)
            and (not data_storage.data_storage_full())
            and (measurement_session)
        ):
            self.operating_mode = 5

        else:
            # Else try to communicate
            if (eps.battery_level >= eps.com_threshold) and (com_window):
                self.operating_mode = 3
            else:
                # Else stop charging if in eclipse or battery is full
                if (eclipse_status) or (eps.battery_level >= eps.max_battery):
                    self.operating_mode = 0
                else:
                    self.operating_mode = 2

    def mode_switch_from_UHF_COM(
        self,
        eps: Eps,
        telecom: Telecom,
        data_storage: DataStorage,
        com_window: bool,
    ) -> None:
        """Handle transitions from UHF-COM mode."""
        # X-BAND DOWNLINK in priority if possible
        if (
            (eps.battery_level >= eps.xb_threshold)
            and (telecom.handshake())
            and (com_window)
            and not data_storage.data_storage_empty(type="x_band")
            and (telecom.can_downlink())
        ):
            self.operating_mode = 4
        else:
            # Else try to go to idle
            if (not com_window) or (telecom.com_finished()):
                self.operating_mode = 0
            else:
                self.operating_mode = 3  # Stay in current mode

    def mode_switch_from_X_BAND_COM(
        self,
        data_storage: DataStorage,
        com_window: bool,
    ) -> None:
        """Handle transitions from XBAND-COM mode."""
        # Go back to UHF-COM if downlink is finished (if COM is still possible)
        if com_window and data_storage.x_band_data_empty():
            self.operating_mode = 3
        else:
            # Else try to go to idle
            if not com_window:
                self.operating_mode = 0
            else:
                self.operating_mode = 4  # Stay in x-band

    def mode_switch_from_MEASUREMENT(
        self, payload: Payload, data_storage: DataStorage
    ) -> None:
        """Handle transitions from MEASUREMENT mode."""
        # Go to idle after measurement ends
        if data_storage.data_storage_full() or payload.campaign_finished():
            self.operating_mode = 0
        else:
            self.operating_mode = 5  # Stay in measurement mode!
