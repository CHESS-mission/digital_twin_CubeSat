"""File in which all functions for switching between different operation mode are implemented
"""

from digital_twin.spacecraft.eps import Eps
from digital_twin.spacecraft.telecom import Telecom
from digital_twin.spacecraft.payload import Payload
from digital_twin.constants import mode_dict


class ModeSwitch:
    def __init__(self, init_mode: int) -> None:
        self.mode_dict = mode_dict
        self.operating_mode = int(init_mode)

    def print_operating_mode(self) -> str:
        return self.mode_dict[self.operating_mode]

    def switch_mode(
        self,
        eps: Eps,
        telecom: Telecom,
        payload: Payload,
        com_window: bool,
        eclipse_status: bool,
        measurement_session: bool,
        safe_flag=False,
    ) -> None:
        if safe_flag:
            self.operating_mode = 1
        else:
            match self.operating_mode:
                case 0:
                    self.mode_switch_from_IDLE(
                        eps,
                        telecom,
                        payload,
                        com_window,
                        eclipse_status,
                        measurement_session,
                    )
                case 1:
                    self.mode_switch_from_SAFE()
                case 2:
                    self.mode_switch_from_CHARGING(
                        eps,
                        telecom,
                        payload,
                        com_window,
                        eclipse_status,
                        measurement_session,
                    )
                case 3:
                    self.mode_switch_from_UHF_COM(eps, telecom, payload, com_window)
                case 4:
                    self.mode_switch_from_X_BAND_COM(telecom, payload, com_window)
                case 5:
                    self.mode_switch_from_MEASUREMENT(telecom, payload)
                case _:
                    raise AssertionError(
                        f"The current operating mode ({self.operating_mode}) does not exist"
                    )

    @property
    def operating_mode(self) -> int:
        return self._operating_mode

    @operating_mode.setter
    def operating_mode(self, mode):
        self._operating_mode = mode

    def mode_switch_from_IDLE(
        self,
        eps: Eps,
        telecom: Telecom,
        payload: Payload,
        com_window: bool,
        eclipse_status: bool,
        measurement_session: bool,
    ):
        # Charging in priority if low battery (if possible)
        if not (eclipse_status) and (eps.battery_level <= eps.min_battery):
            self.operating_mode = 2
        else:
            # Else try to measure
            if (
                (eps.battery_level >= eps.measure_threshold)
                and (not payload.data_storage_full())
                and (measurement_session)
            ):
                self.operating_mode = 5
                # eps.MeasurementDuration = 0
                # eps.DownlinkComplete = false
                # eps.CampaignFinished = false
            else:
                # Else try to communicate
                if (eps.battery_level >= eps.com_threshold) and (com_window):
                    self.operating_mode = 3
                    # eps.OperatingMode = eps.UHF_COM
                    # eps.ComFinished = false
                    # eps.ComDuration = 0
                else:
                    # Else try to charge (if possible)
                    if (not eclipse_status) and (eps.battery_level < eps.max_battery):
                        self.operating_mode = 2
                    else:
                        self.operating_mode = 0

    def mode_switch_from_SAFE(self):
        # the case of a safe flag was treated in the switch_mode() function already
        # so, if the safe flag disappeared, switch to IDLE mode
        self.operating_mode = 0

    def mode_switch_from_CHARGING(
        self,
        eps: Eps,
        telecom: Telecom,
        payload: Payload,
        com_window: bool,
        eclipse_status: bool,
        measurement_session: bool,
    ):
        # Try to measure
        if (
            (eps.battery_level >= eps.measure_threshold)
            and (not payload.data_storage_full())
            and (measurement_session)
        ):
            self.operating_mode = 5
            # eps.MeasurementDuration = 0
            # eps.DownlinkComplete = false
            # eps.CampaignFinished = false

        else:
            # Else try to communicate
            if (eps.battery_level >= eps.com_threshold) and (com_window):
                self.operating_mode = 3
                # eps.OperatingMode = eps.UHF_COM
                # eps.ComFinished = false
                # eps.ComDuration = 0
            else:
                # Else stop charging if in eclipse or battery is full
                if (eclipse_status) or (eps.battery_level >= eps.max_battery):
                    self.operating_mode = 0
                else:
                    self.operating_mode = 2

    def mode_switch_from_UHF_COM(
        self, eps: Eps, telecom: Telecom, payload: Payload, com_window: bool
    ):
        # X-BAND DOWNLINK in priority if possible
        if (
            (eps.battery_level >= eps.xb_threshold)
            and (telecom.handshake())
            and (com_window)
            and (not telecom.downlink_complete(payload))
        ):
            self.operating_mode = 4
            # eps.data_storage_full() = False
        else:
            # Else try to go to idle
            if (not com_window) or (telecom.com_finished()):
                self.operating_mode = 0
                # eps.com_duration = 0
                # eps.com_finished() = False
            else:
                self.operating_mode = 3  # stays un current mode
                # Increment communication duration
                # eps.com_duration += dt
                # if eps.com_duration >= eps.com_max_duration:
                # eps.com_finished() = True

    def mode_switch_from_X_BAND_COM(
        self, telecom: Telecom, payload: Payload, com_window: bool
    ):
        # Go back to UHF-COM if downlink is finished (if COM is still possible)
        if com_window and telecom.downlink_complete(payload):
            self.operating_mode = 3
        else:
            # Else try to go to idle
            if not com_window:
                self.operating_mode = 0
            else:
                self.operating_mode = 4
                # decrement data storage
                # eps.data_storage -= eps.x_band_rate * dt  # adapt to data rate
                # if eps.data_storage <= 0:
                # eps.downlink_complete() = True
                # eps.campaign_finished() = False
                # eps.data_storage = 0

    def mode_switch_from_MEASUREMENT(self, telecom: Telecom, payload: Payload):
        # Go to idle after measurement normally ends
        if payload.data_storage_full() or payload.campaign_finished():
            self.operating_mode = 0
            # eps.measurement_duration = 0
        else:
            self.operating_mode = 5  # stay in measurement mode!
            # Increment data stored
            # eps.data_storage += eps.measure_rate * dt  # adapt to data rate of ScienceInstrument
            # if eps.data_storage >= eps.max_storage:
            # eps.data_storage_full() = True

            # Increment campaign duration
            # eps.measurement_duration += dt
            # if eps.measurement_duration >= eps.science_max_duration:
            # eps.campaign_finished() = True
