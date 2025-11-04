"""Utilitiy functions."""

import json
import os
import shutil
from typing import Any

from astropy import units as u
from astropy.units import Unit
from poliastro.core.elements import rv2coe
import numpy as np

from astropy.coordinates import (
    GCRS,
    ITRS,
    CartesianRepresentation,
    SphericalRepresentation,
)
from astropy.time import Time

from digital_twin.constants import earth_R, earth_k


def get_astropy_unit_time(unit_string: str) -> Unit:
    """Convert a time unit string to an instance of the astropy Unit class."""
    units = {"second": u.s, "hour": u.h, "day": u.day, "year": u.year}
    return units[unit_string]


def get_astropy_units_angle(unit_string: str) -> Unit:
    """Convert an angle unit string to an instance of the astropy Unit class."""
    units = {"degree": u.deg, "radian": u.rad}
    return units[unit_string]


def check_and_empty_folder(folder_path: str) -> None:
    """For a specific folder, check if it exists (if not, create it) and empty it."""
    # If the folder does not exist, create it
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    else:
        # Check if the folder is empty
        if len(os.listdir(folder_path)) > 0:
            # Remove all files and subdirectories in the folder
            for filename in os.listdir(folder_path):
                file_path = os.path.join(folder_path, filename)
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)


def extract_propagation_data_from_ephemeris(eph: np.ndarray) -> tuple[np.ndarray]:
    """Uses Poliastro rv2coe() function to extract orbital elements at each timestep using position data.

    Args:
        eph (np.ndarray): Position data

    Returns:
        np.array (9): Orbital elements
            rr: position
            vv: velocity
            SMAs: semi-major aixs
            ECCs: eccentricity
            INCs: inclination
            RAANs: right ascension of the ascending node
            AOPs: Argument of Periapsis
            TAs: True Anomaly
            altitudes: altitude from attractor's surface
    """
    rr = eph[:, :3]
    vv = eph[:, 3:]
    
    # Create vectorized version of rv2coe
    def rv2coe_wrapper(r, v):
        return np.array(rv2coe(earth_k, r, v))
    
    vectorized_rv2coe = np.vectorize(
        rv2coe_wrapper, 
        signature='(3),(3)->(6)'
    )
    orbital_params = vectorized_rv2coe(rr, vv)
    
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

    return (rr, vv, SMAs, ECCs, INCs, RAANs, AOPs, TAs, altitudes)


def angle_between_vectors(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """
    Calculates the angle (in degrees) between corresponding 3D vectors
    in matrices A and B. Each row represents a vector in 3D space.

    Args:
        A : (np.ndarray): Matrix of shape (n, 3) representing n 3D vectors.
        B : (np.ndarray): Matrix of shape (n, 3) representing n 3D vectors.

    Returns:
        np.ndarray: A vector of shape (n,) containing the angles (in degrees) between corresponding vectors.
    """
    # Dot product of corresponding rows in A and B
    dot_product = np.einsum("ij,ij->i", A, B)

    # Norm (magnitude) of vectors in A and B
    norm_A = np.linalg.norm(A, axis=1)
    norm_B = np.linalg.norm(B, axis=1)

    # Cosine of the angle
    cos_theta = dot_product / (norm_A * norm_B)

    # Clip to avoid numerical errors leading to values outside the range [-1, 1]
    cos_theta = np.clip(cos_theta, -1.0, 1.0)

    # Calculate the angle in radians and then convert to degrees
    angle_rad = np.arccos(cos_theta)  # KEEP values between 0 and pi/2
    angle_deg = np.degrees(angle_rad)

    return angle_deg


def parse_data_file(file_path: str) -> dict:
    """Parse JSON file provided by user

    Args:
        file_path (str): path to JSON file

    Returns:
        dict: Data from file in a dictionary
    """
    f = open(file_path)
    data = json.load(f)
    f.close()
    return data
def convert_cartesian_to_spherical(rr: np.ndarray, epochs: np.ndarray, target_seconds:np.ndarray = None) -> tuple[np.ndarray]:
    """Convert Cartesian coordinates to latitude and longitude in degrees and optionally interpolate to target times."""

    # Convert trajectory to ITRS at original times, then interpolate lat/lon directly
    raw_xyz = CartesianRepresentation(rr, xyz_axis=-1)
    raw_obstime = Time(epochs)

    # transform once to ITRS and get spherical representation (lat/lon) at raw times
    gcrs_xyz = GCRS(raw_xyz, obstime=raw_obstime, representation_type=CartesianRepresentation)
    itrs_xyz = gcrs_xyz.transform_to(ITRS(obstime=raw_obstime))
    itrs_sph = itrs_xyz.represent_as(SphericalRepresentation)

    lat = itrs_sph.lat.to(u.deg).value
    lon = itrs_sph.lon.to(u.deg).value

    # interpolate to target times if provided
    if target_seconds is not None:
        raw_seconds = (raw_obstime - raw_obstime[0]).sec
        # interpolate latitude directly
        lat = np.interp(target_seconds, raw_seconds, lat)

        # handle longitude wrap-around by unwrapping in radians, interpolating, then re-wrapping
        lon = np.rad2deg(np.interp(target_seconds, raw_seconds, np.unwrap(np.deg2rad(lon))))
        lon = ((lon + 180.0) % 360.0) - 180.0  # normalize to [-180, 180)

    return lat, lon