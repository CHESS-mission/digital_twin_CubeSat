"""file that defines the wrapper around poliastro orbit propagator."""

import copy

from astropy import units as u
from astropy.coordinates import (
    GCRS,
    ITRS,
    get_body_barycentric_posvel,
)
from astropy.time import Time
from astropy.units import Quantity
import numpy as np
from poliastro.core.events import eclipse_function
from poliastro.core.perturbations import (
    J2_perturbation,
    atmospheric_drag,
)
from poliastro.core.propagation import func_twobody
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import CowellPropagator

from digital_twin.constants import earth_R, J2, sun_R, earth_k
from digital_twin.orbit_propagator import (
    AtmosphereModelCOESA76,
    AtmosphereModelNRLMSISE00,
    AtmosphereModelSolarActivity,
    AtmosphereModelJB2008,
)
from digital_twin.orbit_propagator.constants import (
    attractor,
    star_string,
    attractor_string,
)
from digital_twin.utils import angle_between_vectors


class OrbitPropagator:
    """
    Class for propagating the orbit of a satellite in the Earth's atmosphere, taking into account
    the perturbative forces of atmospheric drag and the J2 effect.

    The class provides functionality to propagate the orbit of the satellite using Cowell's method
    with perturbations, including atmospheric drag and J2 effect. It updated the air density based
    on an atmospheric model chosen by the user every constant interval of time (also chosen by the
    user). It includes capabilities for calculating satellite visibility from ground stations and eclipse
    status based on its position relative to the Earth and Sun.

    Attributes:
        init_orbit (Orbit): The initial orbit of the satellite.
        current_orbit (Orbit): The current propagated orbit.
        orbit_type (str): Type of the orbit, e.g., "SSO" for Sun-synchronous orbits.
        atmosphere_model (AtmosphereModel): The model used to calculate the air density at a given altitude.
        update_rho_time_interval (Quantity["time"]): The time interval to update air density.
        countdown_rho (Quantity["time"]): Countdown timer for air density update.
        rho (Quantity["mass density"]): The air density at the satellite's current position.
        drag_params (dict): Parameters for calculating atmospheric drag (C_D and A/m).
        method (CowellPropagator): The numerical method used for orbit propagation.
    """

    def __init__(
        self,
        orbit_params: dict,
        epoch: Time,
        atmosphere_model: str,
        update_air_density_timestep: Quantity["time"],
    ) -> None:
        print("Initializing the propagator...")

        # INITAL CHECKS
        if star_string != "Sun":
            raise ValueError(
                "Simulation only implemented for attractor orbiting about the SUN!"
            )
        if attractor_string != "Earth":
            raise ValueError(
                "Simulation only implemented for satellite orbiting about the EARTH!"
            )

        # INITIALIZE ORBITAL ELEMENTS
        init_SMA = orbit_params["orbital_elements"]["SMA"] * u.km
        init_INC = orbit_params["orbital_elements"]["INC"] * u.deg
        init_ECC = orbit_params["orbital_elements"]["ECC"] * u.one
        init_RAAN = orbit_params["orbital_elements"]["RAAN"] * u.deg
        init_AOP = orbit_params["orbital_elements"]["AOP"] * u.deg
        init_TA = orbit_params["orbital_elements"]["TA"] * u.deg

        self.orbit_type = orbit_params["orbit_type"]

        if self.orbit_type == "SSO":
            self.init_orbit = Orbit.heliosynchronous(
                attractor=attractor,
                a=init_SMA,
                ecc=init_ECC,
                # inc = init_INC,
                raan=init_RAAN,
                argp=init_AOP,
                nu=init_TA,
                epoch=epoch,
            )
        else:
            self.init_orbit = Orbit.from_classical(
                attractor=attractor,
                a=init_SMA,
                ecc=init_ECC,
                inc=init_INC,
                raan=init_RAAN,
                argp=init_AOP,
                nu=init_TA,
                epoch=epoch,
            )

        self.current_orbit = copy.deepcopy(self.init_orbit)

        # Choose atmospheric model
        if atmosphere_model == "nrlmsise00":
            self.atmosphere_model = AtmosphereModelNRLMSISE00()
        elif atmosphere_model == "solar_activity":
            self.atmosphere_model = AtmosphereModelSolarActivity()
        elif atmosphere_model == "jb2008":
            self.atmosphere_model = AtmosphereModelJB2008()
        else:
            self.atmosphere_model = AtmosphereModelCOESA76()

        # Rho is air density
        # Approximation updated every __update_rho_time_interval__
        self.update_rho_time_interval = update_air_density_timestep
        self.countdown_rho = copy.deepcopy(
            self.update_rho_time_interval
        )  # Rho will be updated when countdown reaches 0
        self.rho = self.atmosphere_model.get_density(
            iso_date_str=str(self.current_orbit.epoch),
            position=self.current_orbit.r.value,
        ).to_value(u.kg / u.km**3) * (u.kg / u.km**3)

        self.method = CowellPropagator(f=self.f)  # Numerical propagator

        self.drag_params = {
            "C_D": 0.0 * u.one,
            "A_over_m": 0.0 * (u.km**2 / u.kg),
        }  # For drag calculations, updated every time propagate() is called for dynamic drag calculations

    def __str__(self) -> str:
        """Return a string representation of the initial orbit."""
        return f"Initial orbit: {self.init_orbit} with initial period {self.init_orbit.period}"

    @property
    def r(self) -> np.ndarray:
        return self.current_orbit.r

    @property
    def v(self) -> np.ndarray:
        return self.current_orbit.v

    def get_initial_orbit(self) -> Orbit:
        return self.init_orbit

    def propagate(
        self, delta_t: Quantity["time"], C_D: Quantity, A_over_m: Quantity
    ) -> np.ndarray:
        """Propagate the satellite by 1 timestep.

        Args:
            delta_t (Quantity["time"]): Timestep used to propagate.
            C_D (Quantity): Spacecraft new grad coefficient.
            A_over_m (Quantity): Spacecraft new ratio of mass and area.

        Returns:
            np.ndarray: New position of the satellite.
        """
        # Update drag parameters
        self.drag_params["C_D"] = C_D
        self.drag_params["A_over_m"] = A_over_m
        # Update air density
        if self.countdown_rho > 0 * u.s:
            self.countdown_rho -= delta_t
        else:
            self.update_rho(
                self.current_orbit.r.to(u.km).to_value()
            )  # Update air density based on current altitude
            self.countdown_rho = copy.deepcopy(
                self.update_rho_time_interval
            )  # Restart countdown

        new_orbit = self.current_orbit.propagate(delta_t, method=self.method)

        self.current_orbit = new_orbit
        rv = np.zeros(6)
        rv[:3] = self.current_orbit.r
        rv[3:] = self.current_orbit.v
        return rv

    def calculate_vis_window(
        self, ground_stations: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """For each ground station, check if the satellite is visible.

        Args:
            ground_stations (np.ndarray): Array of ground stations.

        Returns:
            np.ndarray: array of Boolean values for each ground station (visible: True).
        """
        visibility = []
        gs_coords_array = []
        satellite_coords = np.reshape(np.array(self.current_orbit.r), (1, 3))
        obstime = self.current_orbit.epoch

        for gs in ground_stations:
            cartesian_coords = gs.cartesian_coords
            itrs_object = ITRS(cartesian_coords, obstime=obstime)
            gcrs_object = itrs_object.transform_to(GCRS(obstime=obstime))
            gs_coords = np.reshape(np.array(gcrs_object.data.xyz.value), (1, 3))
            gs_coords_array.append(gs_coords)
            distance_vec = satellite_coords - gs_coords
            angle = angle_between_vectors(distance_vec, gs_coords)
            if (
                angle >= (90 * u.deg - gs.elev_angle).value
            ):  # Take into account ground station elevation angle
                visibility.append(False)
            else:
                visibility.append(True)

        return np.array(visibility), np.array(gs_coords_array)

    def calculate_eclipse_status(self) -> tuple[bool, Quantity["length"]]:
        """Calculates if the satellite receives light from the sun.
        If the return value is negative, the satellite is in eclipse and does not receive sunlight.
        """
        # Taking into account umbra and penumbra for now
        epoch = self.current_orbit.epoch
        r = self.current_orbit.r.value
        v = self.current_orbit.v.value

        # Position vector of Sun wrt Solar System Barycenter
        r_sec_ssb = get_body_barycentric_posvel("Sun", epoch)[0]
        r_pri_ssb = get_body_barycentric_posvel("Earth", epoch)[0]
        r_sec = ((r_sec_ssb - r_pri_ssb).xyz << u.km).value  # vector from earth to sun

        eclipse = eclipse_function(
            earth_k.value,
            np.hstack((r, v)),
            r_sec,
            sun_R.value,
            earth_R.value,
            umbra=True,
        )

        if eclipse <= 0:
            return True, r_sec * u.km
        else:
            return False, r_sec * u.km
        # If <= 0, the satellite is in eclipse and doesn't get sunlight

    def update_rho(self, position: np.ndarray) -> None:
        """Update air density using the atmosphere model.

        Args:
            position (np.ndarray): Position of the satellite.
        """
        self.rho = self.atmosphere_model.get_density(
            iso_date_str=str(self.current_orbit.epoch), position=position
        ).to_value(u.kg / u.km**3) * (u.kg / u.km**3)

    # Function used by self.f() to compute acceleration
    def a_d(self, t0, state, k, J2, R, C_D, A_over_m, rho):
        return J2_perturbation(t0, state, k, J2, R) + atmospheric_drag(
            t0, state, k, C_D, A_over_m, rho
        )

    # Function used by propagate()
    def f(self, t0, state, k):
        du_kep = func_twobody(t0, state, k)

        ax, ay, az = self.a_d(
            t0,
            state,
            k,
            J2=J2.value,
            R=earth_R.value,
            C_D=self.drag_params["C_D"].value,
            A_over_m=self.drag_params["A_over_m"].value,
            rho=self.rho.value,
        )
        du_ad = np.array([0, 0, 0, ax, ay, az])

        return du_kep + du_ad

    def save_state(self) -> dict:
        """Save the current state of the orbit. File has the same structure as the input file for the orbit parameters.

        Returns:
            dict: Updated dictionary representing the orbit's state. Can be used to start another simulation from this final state.
        """
        save_orbit_params = {}
        save_orbit_params["orbital_elements"] = {}
        save_orbit_params["orbital_elements"]["SMA"] = self.current_orbit.a.to(
            u.km
        ).to_value()
        save_orbit_params["orbital_elements"]["INC"] = self.current_orbit.inc.to(
            u.deg
        ).to_value()
        save_orbit_params["orbital_elements"]["ECC"] = self.current_orbit.ecc.to_value()
        save_orbit_params["orbital_elements"]["RAAN"] = self.current_orbit.raan.to(
            u.deg
        ).to_value()
        save_orbit_params["orbital_elements"]["AOP"] = self.current_orbit.argp.to(
            u.deg
        ).to_value()
        save_orbit_params["orbital_elements"]["TA"] = self.current_orbit.nu.to(
            u.deg
        ).to_value()
        save_orbit_params["orbit_type"] = self.orbit_type
        save_orbit_params["epoch"] = str(self.current_orbit.epoch)

        return save_orbit_params
