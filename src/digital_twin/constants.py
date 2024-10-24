"""Constant for the simulation."""

from astropy import units as u
from poliastro.bodies import Earth, Sun
from poliastro.earth.atmosphere import COESA76


# EARTH CONSTANTS
earth_R = Earth.R.to(u.km)  # Earth's radius
earth_k = Earth.k.to(u.km**3 / u.s**2)  # GM constant for earth
J2 = Earth.J2

# SUN CONSTANTS
sun_R = Sun.R.to(u.km)
solar_power = 1350 * (u.W / u.m**2)

# ATMOSPHERIC MODEL
atmosphere_model = COESA76()

# SPACECRAFT MODES
mode_dict = {
    0: "IDLE",
    1: "SAFE",
    2: "CHARGING",
    3: "UHF_COM",
    4: "X_BAND_COM",
    5: "MEASUREMENT",
}
