from typing import Dict

from astropy.coordinates import (
    GCRS,
    ITRS,
    get_body_barycentric_posvel,
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
    """Interface for atmosphere models."""

    # this init function is call by all almosphere models
    def __init__(self):
        pass

    def get_density(self, iso_date_str: str, position: np.ndarray) -> Quantity:
        raise (NotImplementedError)

    def __str__(self) -> str:
        raise (NotImplementedError)


class AtmosphereModelCOESA76(AtmosphereModel):
    def __init__(self):
        super(AtmosphereModelCOESA76, self).__init__()
        self.model = COESA76()

    def get_density(self, iso_date_str: str, position: np.ndarray) -> Quantity:
        alt = (np.linalg.norm(position) - earth_R.value) * u.km
        return self.model.density(alt)

    def __str__(self):
        return "COESA76 atmosphere model (Poliastro wrapper)"


class AtmosphereModelSolarActivity(AtmosphereModel):
    def __init__(self):
        super(AtmosphereModelSolarActivity, self).__init__()
        # SOLAR ACTIVITY VALUES
        self.F10_7 = np.load("../data/atmosphere_data/solar_activity/F10.npy")
        self.Ap = np.load("../data/atmosphere_data/solar_activity/Ap.npy")
        self.times = np.load(
            "../data/atmosphere_data/solar_activity/dates_sa.npy", allow_pickle=True
        )

    def get_solar_activity_values(self, iso_date_str: str):
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

    def calculate_density(self, F10_7: float, Ap: float, alt: float):
        T = 900 + 2.5 * (F10_7 - 70) + 1.5 * Ap  # [Kelvin]
        m = 27 - 0.012 * (alt - 200)
        H = T / m  # [km]
        rho = 6 * 10 ** (-10) * np.exp(-(alt - 175) / H)  # [kg m‐3]
        return rho * (u.kg / u.m**3)

    def get_density(self, iso_date_str: str, position: np.ndarray) -> Quantity:
        alt = np.linalg.norm(position) - earth_R.value
        F10_7, Ap = self.get_solar_activity_values(iso_date_str)
        return self.calculate_density(F10_7, Ap, alt)


class AtmosphereModelNRLMSISE00(AtmosphereModel):
    def __init__(self):
        super(AtmosphereModelNRLMSISE00, self).__init__()
        swfile = download_sw_nrlmsise00("../data/atmosphere_data/NRLMSISE00/")
        # Read the space weather data
        self.swdata = read_sw_nrlmsise00(swfile)
        self.t_max = datetime.strptime(
            "2024-12-01 00:00:00.000", "%Y-%m-%d %H:%M:%S.%f"
        )  # Maximum time for space weather data (checked experimentally)
        self.cycle = relativedelta(years=11)  # solar activity cycle duration

    def get_density(self, iso_date_str: str, position: np.ndarray) -> Quantity:
        t_wanted = datetime.strptime(iso_date_str, "%Y-%m-%d %H:%M:%S.%f")

        # go back in time with solar cycles
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


class AtmosphereModelJB2008(AtmosphereModel):
    def __init__(self):
        super(AtmosphereModelJB2008, self).__init__()
        swfile = download_sw_jb2008("../data/atmosphere_data/JB2008/")
        # Read the space weather data
        self.swdata = read_sw_jb2008(swfile)
        self.t_max = datetime.strptime(
            "2024-12-01 00:00:00.000", "%Y-%m-%d %H:%M:%S.%f"
        )  # Maximum time for space weather data
        self.cycle = relativedelta(years=11)  # solar activity cycle duration

    def get_density(self, iso_date_str: str, position: np.ndarray) -> Quantity:
        t_wanted = datetime.strptime(iso_date_str, "%Y-%m-%d %H:%M:%S.%f")

        # go back in time with solar cycles
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
