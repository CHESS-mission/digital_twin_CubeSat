"""File in which all functions for switching between different operation mode are implemented
"""

from digital_twin.spacecraft.eps import Eps
from digital_twin.spacecraft.telecom import Telecom


class ModeSwitch:
    def __init__(self, init_mode: str) -> None:
        self.operating_mode = init_mode

    def switch_mode(
        self,
        eps: Eps,
        telecom: Telecom,
        com_window: bool,
        eclipse_status: bool,
        measurement_session: bool,
        safe_flag=False,
    ) -> None:
        if safe_flag:
            self.operating_mode = "SAFE"
        else:
            match self.operating_mode:
                case "IDLE":
                    self.mode_switch_from_IDLE(
                        eps, telecom, com_window, eclipse_status, measurement_session
                    )
                case "SAFE":
                    self.mode_switch_from_SAFE()
                case "CHARGING":
                    self.mode_switch_from_CHARGING(
                        eps, telecom, com_window, eclipse_status, measurement_session
                    )
                case "UHF_COM":
                    self.mode_switch_from_UHF_COM(eps, telecom, com_window)
                case "X_BAND_COM":
                    self.mode_switch_from_X_BAND_COM(telecom, com_window)
                case "MEASUREMENT":
                    self.mode_switch_from_MEASUREMENT(eps, com_window, eclipse_status)
                case _:
                    raise AssertionError(
                        f"The current operating mode ({self.operating_mode}) does not exist"
                    )

    @property
    def operating_mode(self) -> str:
        return self._operating_mode

    @operating_mode.setter
    def operating_mode(self, mode):
        self._operating_mode = mode

    def mode_switch_from_IDLE(
        self,
        eps: Eps,
        telecom: Telecom,
        com_window: bool,
        eclipse_status: bool,
        measurement_session: bool,
    ):
        # Charging in priority if low battery (if possible)
        if not (eclipse_status) and (eps.battery_level <= eps.min_battery):
            self.operating_mode = "CHARGING"
        else:
            # Else try to measure
            if (
                (eps.battery_level >= eps.measure_threshold)
                and (not telecom.data_storage_full)
                and (measurement_session)
            ):
                self.operating_mode = "MEASUREMENT"
                # eps.MeasurementDuration = 0
                # eps.DownlinkComplete = false
                # eps.CampaignFinished = false
            else:
                # Else try to communicate
                if (eps.battery_level >= eps.com_threshold) and (com_window):
                    self.operating_mode = "UHF_COM"
                    # eps.OperatingMode = eps.UHF_COM
                    # eps.ComFinished = false
                    # eps.ComDuration = 0
                else:
                    # Else try to charge (if possible)
                    if (not eclipse_status) and (eps.battery_level < eps.max_battery):
                        self.operating_mode = "CHARGING"

    def mode_switch_from_SAFE(self):
        # the case of a safe flag was treated in the switch_mode() function already
        # so, if the safe flag disappeared, switch to IDLE mode
        self.operating_mode = "IDLE"

    def mode_switch_from_CHARGING(
        self,
        eps: Eps,
        telecom: Telecom,
        com_window: bool,
        eclipse_status: bool,
        measurement_session: bool,
    ):
        # Try to measure
        if (
            (eps.battery_level >= eps.measure_threshold)
            and (not telecom.data_storage_full)
            and (measurement_session)
        ):
            self.operating_mode = "MEASUREMENT"
            # eps.MeasurementDuration = 0
            # eps.DownlinkComplete = false
            # eps.CampaignFinished = false

        else:
            # Else try to communicate
            if (eps.battery_level >= eps.com_threshold) and (com_window):
                self.operating_mode = "UHF_COM"
                # eps.OperatingMode = eps.UHF_COM
                # eps.ComFinished = false
                # eps.ComDuration = 0
            else:
                # Else stop charging if in eclipse or battery is full
                if (eclipse_status) or (eps.battery_level >= eps.max_battery):
                    self.operating_mode = "IDLE"

    def mode_switch_from_UHF_COM(self, eps: Eps, telecom: Telecom, com_window: bool):
        # X-BAND DOWNLINK in priority if possible
        if (
            (eps.battery_level >= eps.xb_threshold)
            and (telecom.handshake)
            and (com_window)
            and (not telecom.downlink_complete)
        ):
            self.operating_mode = "X_BAND_COM"
            # eps.data_storage_full = False
        else:
            # Else try to go to idle
            if (not com_window) or (telecom.com_finished):
                self.operating_mode = "IDLE"
                # eps.com_duration = 0
                # eps.com_finished = False
            # else:
            # Increment communication duration
            # eps.com_duration += dt
            # if eps.com_duration >= eps.com_max_duration:
            # eps.com_finished = True

    def mode_switch_from_X_BAND_COM(self, telecom: Telecom, com_window: bool):
        # Go back to UHF-COM if downlink is finished (if COM is still possible)
        if com_window and telecom.downlink_complete:
            self.operating_mode = "UHF_COM"
        else:
            # Else try to go to idle
            if not com_window:
                self.operating_mode = "IDLE"
            # else:
            # Else decrement data storage
            # eps.data_storage -= eps.x_band_rate * dt  # adapt to data rate
            # if eps.data_storage <= 0:
            # eps.downlink_complete = True
            # eps.campaign_finished = False
            # eps.data_storage = 0

    def mode_switch_from_MEASUREMENT(
        self, eps: Eps, com_window: bool, eclipse_status: bool
    ):
        pass

    def mode_switch_from_MEASUREMENT(
        self,
        eps: Eps,
        telecom: Telecom,
        com_window: bool,
        eclipse_status: bool,
        dt: float,
    ):
        # Go to idle after measurement normally ends
        if telecom.data_storage_full or telecom.campaign_finished:
            self.operating_mode = "IDLE"
            # eps.measurement_duration = 0
        # else:
        # Increment data stored
        # eps.data_storage += eps.measure_rate * dt  # adapt to data rate of ScienceInstrument
        # if eps.data_storage >= eps.max_storage:
        # eps.data_storage_full = True

        # Increment campaign duration
        # eps.measurement_duration += dt
        # if eps.measurement_duration >= eps.science_max_duration:
        # eps.campaign_finished = True
