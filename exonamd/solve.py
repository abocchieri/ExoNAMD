import numpy as np
import astropy.constants as cc
import astropy.units as u



def solve_a_rs(sma, rstar, ars):
    """
    Solve semi-major axis -- stellar radius system of equations.

    Parameters
    ----------
    sma: Float
        Orbital semi-major axis (AU).
    rstar: Float
        Stellar radius (r_sun).
    ars: Float
        sma--rstar ratio.
    """
    missing = np.isnan(sma) + np.isnan(ars) + np.isnan(rstar)
    # Know everything or not enough:
    if missing != 1:
        return sma, rstar, ars

    if np.isnan(sma):
        sma = ars * rstar * u.R_sun
        sma = sma.to(u.au).value
    elif np.isnan(rstar):
        rstar = sma * u.au / ars
        rstar = rstar.to(u.R_sun).value
    elif np.isnan(ars):
        ars = sma * u.au / (rstar * u.R_sun).to(u.au)
        ars = ars.value

    return sma, rstar, ars


def solve_rprs(rplanet, rstar, rprs):
    """
    Solve planet radius -- stellar radius system of equations.

    Parameters
    ----------
    rplanet: Float
        Planet radius (r_earth).
    rstar: Float
        Stellar radius (r_sun).
    rprs: Float
        Planet--star radius ratio.
    """
    missing = np.isnan(rplanet) + np.isnan(rstar) + np.isnan(rprs)
    # Know everything or not enough:
    if missing != 1:
        return rplanet, rstar, rprs

    if np.isnan(rplanet):
        rplanet = rprs * (rstar * u.R_sun)
        rplanet = rplanet.to(u.R_earth).value
    elif np.isnan(rstar):
        rstar = rplanet * u.R_earth / rprs
        rstar = rstar.to(u.R_sun).value
    elif np.isnan(rprs):
        rprs = rplanet * u.R_earth / (rstar * u.R_sun).to(u.R_earth)
        rprs = rprs.value

    return rplanet, rstar, rprs


def solve_a_period(period, sma, mstar):
    """
    Solve period-sma-mstar system of equations.

    Parameters
    ----------
    period: Float
        Orbital period (days).
    sma: Float
        Orbital semi-major axis (AU).
    mstar: Float
        Stellar mass (m_sun).
    """
    missing = np.isnan(period) + np.isnan(sma) + np.isnan(mstar)
    # Know everything or not enough:
    if missing != 1:
        return period, sma, mstar

    two_pi_G = 2.0 * np.pi / np.sqrt(cc.G)
    if np.isnan(mstar):
        mstar = (sma * u.au) ** 3.0 / (period * u.day / two_pi_G) ** 2.0
        mstar = mstar.to(u.M_sun).value
    elif np.isnan(period):
        period = np.sqrt((sma * u.au) ** 3.0 / (mstar * u.M_sun)) * two_pi_G
        period = period.to(u.day).value
    elif np.isnan(sma):
        sma = ((period * u.day / two_pi_G) ** 2.0 * (mstar * u.M_sun)) ** (1 / 3)
        sma = sma.to(u.au).value

    return period, sma, mstar


def solve_values(row):

    sma = row["pl_orbsmax"]
    ars = row["pl_ratdor"]
    rstar = row["st_rad"]
    rplanet = row["pl_rade"]
    rprs = row["pl_ratror"]
    period = row["pl_orbper"]
    mstar = row["st_mass"]

    # Rank groups
    a_rs_ = np.isnan(sma) + np.isnan(ars) + np.isnan(rstar)
    rprs_ = np.isnan(rplanet) + np.isnan(rprs) + np.isnan(rstar)
    a_period_ = np.isnan(period) + np.isnan(sma) + np.isnan(mstar)
    solve_order = np.argsort([a_rs_, rprs_, a_period_])
    for i in solve_order:
        if i == 0:
            # Solve semi-major axis -- stellar radius system of equations.
            solution = solve_a_rs(sma, rstar, ars)
            sma, rstar, ars = solution
        elif i == 1:
            # Solve planet radius -- stellar radius system of equations.
            solution = solve_rprs(rplanet, rstar, rprs)
            rplanet, rstar, rprs = solution
        elif i == 2:
            # Solve period-sma-mstar system of equations.
            solution = solve_a_period(period, sma, mstar)
            period, sma, mstar = solution

    return sma, ars, rstar, rplanet, rprs, period, mstar
