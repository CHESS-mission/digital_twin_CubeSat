from astropy import units as u
from astropy.units import Unit


def get_astropy_unit(unit_string: str) -> Unit:
    units = {"second": u.s, "hour": u.h, "day": u.day, "year": u.year}
    return units[unit_string]
