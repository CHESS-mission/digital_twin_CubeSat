import os
import shutil

import numpy as np

from astropy import units as u
from astropy.units import Unit
from poliastro.core.elements import rv2coe

from digital_twin.constants import earth_R, earth_k


def get_astropy_unit_time(unit_string: str) -> Unit:
    units = {"second": u.s, "hour": u.h, "day": u.day, "year": u.year}
    return units[unit_string]


def get_astropy_units_angle(unit_string: str) -> Unit:
    units = {"degree": u.deg, "radian": u.rad}
    return units[unit_string]


def check_and_empty_folder(folder_path: str) -> None:
    # Check if folder exists, if not, create it
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    else:
        # Check if the folder is empty
        if len(os.listdir(folder_path)) > 0:
            # Remove all files and subdirectories in the folder
            for filename in os.listdir(folder_path):
                file_path = os.path.join(folder_path, filename)
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)  # Remove file or symlink
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)  # Remove directory and its contents


def extract_propagation_data_from_ephemeris(eph: np.array):
    rr = eph[:, :3]
    vv = eph[:, 3:]
    orbital_params = np.array([rv2coe(earth_k, r, v) for r, v in zip(rr, vv)])
    ps = orbital_params[:, 0]
    ECCs = orbital_params[:, 1]
    INCs = orbital_params[:, 2]
    RAANs = orbital_params[:, 3]
    AOPs = orbital_params[:, 4]
    TAs = orbital_params[:, 5]
    SMAs = np.divide(
        ps, 1 - np.multiply(ECCs, ECCs)
    )  # Formula linking semi-latus rectum to semi-major axis
    altitudes = np.linalg.norm(eph, axis=1) - earth_R.value

    return rr, vv, SMAs, ECCs, INCs, RAANs, AOPs, TAs, altitudes
