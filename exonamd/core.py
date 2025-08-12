import numpy as np
import pandas as pd


def compute_amdk(m, e, di, a):
    """
    Compute the angular momentum deficit (AMD) for a planet in a system or for multiple planets.

    Parameters
    ----------
    m : float, np.array
        Mass of the planet.
    e : float, np.array
        Eccentricity of the planet's orbit.
    di : float, np.array
        Relative angle (inclination or obliquity) of the planet's orbit in degrees w.r.t. a reference (e.g. the most massive planet in the system).
    a : float, np.array
        Semi-major axis of the planet's orbit.

    Returns
    -------
    amdk : float
        Angular momentum deficit.
    """
    return m * np.sqrt(a) * (1 - np.sqrt(1 - e**2) * np.cos(np.deg2rad(di)))


def compute_namd(amdk, m, sqrt_a):
    """
    Compute the normalized angular momentum deficit (NAMD) for a given system.

    Parameters
    ----------
    amdk : np.array
        Angular momentum deficit for the planets in a system.
    m : np.array
        Mass of the planets in the system
    sqrt_a : np.array
        Square root of the semi-major axis of the planets in the system

    Returns
    -------
    namd : float
        Normalized angular momentum deficit.
    """

    return np.sum(amdk) / np.sum(m * sqrt_a)
