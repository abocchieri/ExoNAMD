import numpy as np
import astropy.constants as cc
import astropy.units as u
from spright import RMRelation


rmr = RMRelation()


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


def solve_eccentricity(row):
    ecc = row["pl_orbeccen"]
    eccerr1 = row["pl_orbeccenerr1"]
    eccerr2 = row["pl_orbeccenerr2"]
    flag = row["flag"]

    if np.isnan(ecc):
        ecc = 0.63 * row["sy_pnum"] ** (-1.02)
        eccerr1 = 0.0
        eccerr2 = 0.0
        flag += "1+-"

    if np.isnan(eccerr1):
        eccerr1 = 0.0
        flag += "1+"

    if np.isnan(eccerr2):
        eccerr2 = 0.0
        flag += "1-"

    return ecc, eccerr1, eccerr2, flag


def solve_mass(row, min_radius=0.5, max_radius=6.0):

    mass = row["pl_bmasse"]
    masserr1 = row["pl_bmasseerr1"]
    masserr2 = row["pl_bmasseerr2"]
    radius = row["pl_rade"]
    radiuserr1 = row["pl_radeerr1"]
    radiuserr2 = row["pl_radeerr2"]
    flag = row["flag"]

    if (
        np.isnan(mass)
        and radius > min_radius
        and radius < max_radius
        and not np.isnan(radiuserr1)
        and not np.isnan(radiuserr2)
    ):
        mds = rmr.predict_mass(radius=(radius, 0.5 * (radiuserr1 - radiuserr2)))
        q16, q50, q84 = np.percentile(mds.samples, [16, 50, 84])
        mass = q50
        masserr1 = q84 - q50
        masserr2 = q16 - q50
        flag += "2+-"

    if np.isnan(masserr1):
        masserr1 = 0.0
        flag += "2+"

    if np.isnan(masserr2):
        masserr2 = 0.0
        flag += "2-"

    return mass, masserr1, masserr2, flag


def solve_inclination(row, df):

    hostname = row["hostname"]
    incl = row["pl_orbincl"]
    inclerr1 = row["pl_orbinclerr1"]
    inclerr2 = row["pl_orbinclerr2"]
    flag = row["flag"]

    host = df[df["hostname"] == hostname]
    mass_max = host["pl_bmasse"].idxmax()
    inclnan = host[host["pl_orbincl"].isnull()]

    # if inclination is not nan, check the errors
    # if they are nan, set them to 0
    if not np.isnan(incl):
        if np.isnan(inclerr1):
            inclerr1 = 0.0
            flag += "3+"
        if np.isnan(inclerr2):
            inclerr2 = 0.0
            flag += "3-"

    elif len(inclnan) == len(host):  
        # i.e., all inclinations in the system are nan
        incl = 90.0
        inclerr1 = 0.0
        inclerr2 = 0.0
        flag += "3+-"

    elif mass_max in inclnan.index:
        # i.e., the most massive planet has inclination nan
        # Note: we could change this to the next most massive planet later on
        incl = 90.0
        inclerr1 = 0.0
        inclerr2 = 0.0
        flag += "3+-"

    else:
        # i.e., the most massive planet reports inclination
        incl = host.loc[mass_max, "pl_orbincl"]
        inclerr1 = host.loc[mass_max, "pl_orbinclerr1"]
        inclerr2 = host.loc[mass_max, "pl_orbinclerr2"]
        flag += "3+-"

    return incl, inclerr1, inclerr2, flag


def solve_sma(row):

    smaerr1 = row["pl_orbsmaxerr1"]
    smaerr2 = row["pl_orbsmaxerr2"]
    flag = row["flag"]

    if np.isnan(smaerr1):
        smaerr1 = 0.0
        flag += "4+"

    if np.isnan(smaerr2):
        smaerr2 = 0.0
        flag += "4-"

    return smaerr1, smaerr2, flag