from poliastro.bodies import Earth, Sun
from poliastro.earth.atmosphere import COESA76

from astropy import units as u

# EARTH CONSTANTS
earth_R = Earth.R.to(u.km)  # Earth's radius (km)
earth_k = Earth.k.to(u.km**3 / u.s**2)  # GM constant for earth
J2 = Earth.J2

# SUN CONSTANTS
sun_R = Sun.R.to(u.km)

# ATMOSPHERIC MODEL
atmosphere_model = COESA76()
