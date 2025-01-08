"""
Atmosphere models for calculating atmospheric density based on different models and solar activity data.

This module defines various atmosphere models to compute the mass density at a given altitude, considering
solar activity and space weather. The models include COESA76, NRLMSISE00, and JB2008, each providing different
ways to calculate atmospheric density based on the provided date and position.

Classes:
    AtmosphereModel: Base class defining the interface for different atmosphere models.
    AtmosphereModelCOESA76: Implements the COESA76 model for atmospheric density calculation.
    AtmosphereModelSolarActivity: Implements a model using solar activity values (F10.7, Ap) to calculate density.
    AtmosphereModelNRLMSISE00: Implements the NRLMSISE00 model for atmospheric density calculation.
    AtmosphereModelJB2008: Implements the JB2008 model for atmospheric density calculation.
"""

from astropy.coordinates import (
    GCRS,
    ITRS,
    SphericalRepresentation,
    CartesianRepresentation,
)
from astropy.time import Time
import astropy.units as u
from astropy.units import Quantity
from datetime import datetime
from dateutil.relativedelta import relativedelta
import numpy as np
from poliastro.earth.atmosphere import COESA76
from pyatmos import (
    download_sw_nrlmsise00,
    read_sw_nrlmsise00,
    download_sw_jb2008,
    read_sw_jb2008,
)
from pyatmos import nrlmsise00, jb2008

from digital_twin.constants import earth_R


class AtmosphereModel:
    """
    Interface for atmosphere models.

    This is an abstract base class, and each subclass must implement the `get_density` method to compute
    the atmospheric density for a given date and position in space.
    """

    # this init function is called by all almosphere models
    def __init__(self, verbose) -> None:
        print("Initializing the atmosphere model... ") if verbose else None
        pass

    def get_density(
        self, iso_date_str: str, position: np.ndarray
    ) -> Quantity["mass density"]:
        raise (NotImplementedError)

    def __str__(self) -> str:
        """Return a string representation of the atmosphere model."""
        raise (NotImplementedError)


class Coesa76(AtmosphereModel):
    """
    COESA76 model for atmospheric density calculation.

    This model is used for calculating atmospheric density based on the COESA76 model, which is
    commonly used as a "simple" model for calculating atmospheric properties at various altitudes.
    """

    def __init__(self, verbose: bool = False) -> None:
        super(Coesa76, self).__init__(verbose)
        self.model = COESA76()

    def get_density(
        self, iso_date_str: str, position: np.ndarray
    ) -> Quantity["mass density"]:
        """Calculate the atmospheric density based on the given position (this model does not use the date)"""
        alt = (np.linalg.norm(position) - earth_R.value) * u.km
        return self.model.density(alt)

    def __str__(self) -> str:
        """Return a string representation of the atmosphere model."""
        return "Coesa76 atmosphere model (Poliastro wrapper)"


class SolarActivity(AtmosphereModel):
    """
    Solar activity-based model for atmospheric density calculation.

    This model calculates atmospheric density based on solar activity values (F10.7 and Ap indices)
    and altitude. Solar activity values are stored in preloaded arrays, and density is calculated
    using empirical formulas.
    """

    def __init__(
        self,
        path_f10: str = "data/atmosphere_data/solar_activity/F10.npy",
        path_ap: str = "data/atmosphere_data/solar_activity/Ap.npy",
        path_times: str = "data/atmosphere_data/solar_activity/dates_sa.npy",
        verbose: bool = False,
    ) -> None:
        super(SolarActivity, self).__init__(verbose)
        # SOLAR ACTIVITY VALUES
        self.F10_7 = np.load(path_f10)
        self.Ap = np.load(path_ap)
        self.times = np.load(path_times, allow_pickle=True)

    def get_solar_activity_values(self, iso_date_str: str) -> tuple[float]:
        """Retrieve the solar activity values (F10.7 and Ap indices) for a given date."""
        target_date = Time(
            iso_date_str, format="iso"
        ).ymdhms  # Extracts year, month, day tuple
        for i, time in enumerate(self.times):
            time_date = time.ymdhms
            if (time_date.year, time_date.month, time_date.day) == (
                target_date.year,
                target_date.month,
                target_date.day,
            ):
                return self.F10_7[i], self.Ap[i]
        raise ValueError(f"Date {iso_date_str} not found in times array.")

    def calculate_density(
        self, F10_7: float, Ap: float, alt: float
    ) -> Quantity["mass density"]:
        """Calculate the atmospheric density based on solar activity values and altitude."""
        T = 900 + 2.5 * (F10_7 - 70) + 1.5 * Ap  # [Kelvin]
        m = 27 - 0.012 * (alt - 200)
        H = T / m  # [km]
        rho = 6 * 10 ** (-10) * np.exp(-(alt - 175) / H)  # [kg m‐3]
        return rho * (u.kg / u.m**3)

    def get_density(
        self, iso_date_str: str, position: np.ndarray
    ) -> Quantity["mass density"]:
        """Calculate the atmospheric density at a given position and date using solar activity data."""
        alt = np.linalg.norm(position) - earth_R.value
        F10_7, Ap = self.get_solar_activity_values(iso_date_str)
        return self.calculate_density(F10_7, Ap, alt)

    def __str__(self) -> str:
        """Return a string representation of the atmosphere model."""
        return "Analytical atmosphere model based on solar activity data"


class NRLMSISE00(AtmosphereModel):
    """
    NRLMSISE-00 model for atmospheric density calculation.

    This model calculates atmospheric density based on the NRLMSISE-00 empirical model, which incorporates
    space weather data such as solar activity to estimate the density at a given altitude.
    """

    def __init__(
        self, path: str = "data/atmosphere_data/NRLMSISE00/", verbose: bool = False
    ) -> None:
        super(NRLMSISE00, self).__init__(verbose)
        swfile = download_sw_nrlmsise00(path)
        # Read the space weather data
        self.swdata = read_sw_nrlmsise00(swfile)
        self.t_max = datetime.strptime(
            "2024-12-01 00:00:00.000", "%Y-%m-%d %H:%M:%S.%f"
        )  # Maximum time for space weather data (checked experimentally)
        self.cycle = relativedelta(years=11)  # Solar activity cycle duration

    def get_density(
        self, iso_date_str: str, position: np.ndarray
    ) -> Quantity["mass density"]:
        """Calculate the atmospheric density at a given position and date using the NRLMSISE-00 model."""
        t_wanted = datetime.strptime(iso_date_str, "%Y-%m-%d %H:%M:%S.%f")

        # Do back in time with solar cycles in order to ask for a valid date
        t = t_wanted
        while t > self.t_max:
            t = t - self.cycle
        epoch = t.strftime("%Y-%m-%d %H:%M:%S.%f")

        raw_xyz = CartesianRepresentation(position)
        raw_obstime = epoch
        gcrs_xyz = GCRS(
            raw_xyz, obstime=raw_obstime, representation_type=CartesianRepresentation
        )
        itrs_xyz = gcrs_xyz.transform_to(ITRS(obstime=raw_obstime))
        itrs_latlon = itrs_xyz.represent_as(SphericalRepresentation)
        lat = itrs_latlon.lat.to_value(u.deg)
        lon = itrs_latlon.lon.to_value(u.deg)
        alt = np.linalg.norm(position) - earth_R.value

        nrl00 = nrlmsise00(t, (lat, lon, alt), self.swdata)
        return nrl00.rho * (u.kg / u.m**3)

    def __str__(self) -> str:
        """Return a string representation of the atmosphere model."""
        return "NRLMSISE00 atmosphere model"


class JB2008(AtmosphereModel):
    """
    JB2008 model for atmospheric density calculation.

    This model calculates atmospheric density using the JB2008 empirical model, which is based on
    solar activity and space weather data to estimate the density at a given altitude.
    """

    def __init__(
        self, path: str = "data/atmosphere_data/JB2008/", verbose: bool = False
    ) -> None:
        super(JB2008, self).__init__(verbose)
        swfile = download_sw_jb2008(path)
        # Read the space weather data
        self.swdata = read_sw_jb2008(swfile)
        self.t_max = datetime.strptime(
            "2024-12-01 00:00:00.000", "%Y-%m-%d %H:%M:%S.%f"
        )  # Maximum time for space weather data
        self.cycle = relativedelta(years=11)  # Solar activity cycle duration

    def get_density(
        self, iso_date_str: str, position: np.ndarray
    ) -> Quantity["mass density"]:
        """Calculate the atmospheric density at a given position and date using the JB2008 model."""
        t_wanted = datetime.strptime(iso_date_str, "%Y-%m-%d %H:%M:%S.%f")

        # Do back in time with solar cycles in order to ask for a valid date
        t = t_wanted
        while t > self.t_max:
            t = t - self.cycle
        epoch = t.strftime("%Y-%m-%d %H:%M:%S.%f")

        raw_xyz = CartesianRepresentation(position)
        raw_obstime = epoch
        gcrs_xyz = GCRS(
            raw_xyz, obstime=raw_obstime, representation_type=CartesianRepresentation
        )
        itrs_xyz = gcrs_xyz.transform_to(ITRS(obstime=raw_obstime))
        itrs_latlon = itrs_xyz.represent_as(SphericalRepresentation)
        lat = itrs_latlon.lat.to_value(u.deg)
        lon = itrs_latlon.lon.to_value(u.deg)
        alt = np.linalg.norm(position) - earth_R.value

        jc08 = jb2008(t, (lat, lon, alt), self.swdata)
        return jc08.rho * (u.kg / u.m**3)

    def __str__(self) -> str:
        """Return a string representation of the atmosphere model."""
        return "JB2008 atmosphere model"
