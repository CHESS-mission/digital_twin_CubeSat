"""Constants for the simulation."""

from astropy import units as u
from poliastro.bodies import Earth, Sun

# SIMULATION TIME CONSTANTS
simulation_unit = u.s
simulation_unit_string = "second"

# EARTH CONSTANTS
earth_R = Earth.R.to(u.km)  # Earth's radius
earth_k = Earth.k.to(u.km**3 / u.s**2)  # GM constant for earth
J2 = Earth.J2

# SUN CONSTANTS
sun_R = Sun.R.to(u.km)
solar_flux = 1367 * (u.W / u.m**2)

# SPACECRAFT MODES
mode_dict = {
    0: "IDLE",
    1: "SAFE",
    2: "CHARGING",
    3: "UHF_COM",
    4: "X_BAND_COM",
    5: "MEASUREMENT",
}

# SPACECRAFT ATTITUDES
attitude_dict = {
    0: "thomson spin",
    1: "sun spin",
    2: "nadir spin",
    3: "ground tracking",
    4: "nadir pointing",
}

attitude_mode_dict = {
    0: 0,
    1: 1,
    2: 1,
    3: 2,
    4: 3,
    5: 4,
}  # key is mode, value is corresponding attitude
