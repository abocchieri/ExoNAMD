import numpy as np
import astropy.constants as cc
import astropy.units as u
import pandas as pd

# import modin.pandas as pd

from exonamd.core import compute_amd
from exonamd.core import compute_namd


# Functions solve_a_rs, solve_rprs, solve_a_period, solve_values
# below are originally from gen_tso
# by Patricio Cubillos: https://github.com/pcubillos/gen_tso
# modified to work in this project


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

    return pd.Series([sma, ars, rstar, rplanet, rprs, period, mstar])


def solve_relincl(row, df):
    # Computes the inclination w.r.t. the most massive planet in the system

    hostname = row["hostname"]
    incl = row["pl_orbincl"]
    inclerr1 = row["pl_orbinclerr1"]
    inclerr2 = row["pl_orbinclerr2"]

    host = df[df["hostname"] == hostname]
    mass_max = host["pl_bmasse"].idxmax()
    max_mass_data = host.loc[
        mass_max, ["pl_orbincl", "pl_orbinclerr1", "pl_orbinclerr2"]
    ]

    relincl = max_mass_data["pl_orbincl"] - incl
    relinclerr1 = np.sqrt(inclerr1**2 + max_mass_data["pl_orbinclerr1"] ** 2)
    relinclerr2 = -np.sqrt(inclerr2**2 + max_mass_data["pl_orbinclerr2"] ** 2)

    return pd.Series([relincl, relinclerr1, relinclerr2])


def solve_amd(row, kind: str):

    mass = row["pl_bmasse"]
    eccen = row["pl_orbeccen"]
    di = {"rel": row["pl_relincl"], "abs": row["pl_trueobliq"]}
    di = di[kind]
    sma = row["pl_orbsmax"]

    amd = compute_amd(mass, eccen, di, sma)

    return pd.Series([amd])


def solve_namd(row, df, kind: str):

    if f"namd_{kind}" in row.index:
        return row[f"namd_{kind}"]

    hostname = row["hostname"]
    host = df[df["hostname"] == hostname]

    amd = host[f"amd_{kind}"]
    mass = host["pl_bmasse"]
    sqrt_sma = np.sqrt(host["pl_orbsmax"])

    namd = compute_namd(amd, mass, sqrt_sma)

    return pd.Series([namd])


def solve_amd_mc(row, kind, Npt, threshold, namd=False):

    if not namd and f"amd_{kind}_q50" in row.index:
        return row[[f"amd_{kind}_q16", f"amd_{kind}_q50", f"amd_{kind}_q84"]]

    mass = row["pl_bmasse"]
    masserr1 = row["pl_bmasseerr1"]
    masserr2 = row["pl_bmasseerr2"]
    eccen = row["pl_orbeccen"]
    eccenerr1 = row["pl_orbeccenerr1"]
    eccenerr2 = row["pl_orbeccenerr2"]

    di = {
        "rel": row["pl_relincl"],
        "relerr1": row["pl_relinclerr1"],
        "relerr2": row["pl_relinclerr2"],
        "abs": row["pl_trueobliq"],
        "abserr1": row["pl_trueobliqerr1"],
        "abserr2": row["pl_trueobliqerr2"],
    }

    sma = row["pl_orbsmax"]
    smaerr1 = row["pl_orbsmaxerr1"]
    smaerr2 = row["pl_orbsmaxerr2"]

    # Sample the parameters
    mass_mc = np.random.normal(mass, 0.5 * (masserr1 - masserr2), Npt)
    eccen_mc = np.random.normal(eccen, 0.5 * (eccenerr1 - eccenerr2), Npt)
    di_mc = np.random.normal(
        di[kind], 0.5 * (di[f"{kind}err1"] - di[f"{kind}err2"]), Npt
    )
    sma_mc = np.random.normal(sma, 0.5 * (smaerr1 - smaerr2), Npt)

    # Mask unphysical values
    mask = (
        (mass_mc > 0)
        & (eccen_mc >= 0)
        & (eccen_mc < 1)
        & np.abs(di_mc >= 0)
        & (di_mc < 180)
        & (sma_mc > 0)
    )

    mass_mc = np.ma.MaskedArray(mass_mc, mask=~mask)
    eccen_mc = np.ma.MaskedArray(eccen_mc, mask=~mask)
    di_mc = np.ma.MaskedArray(di_mc, mask=~mask)
    sma_mc = np.ma.MaskedArray(sma_mc, mask=~mask)

    # Check number of valid samples
    if len(mass_mc.compressed()) < threshold:

        out = {
            f"amd_{kind}_q16": np.nan,
            f"amd_{kind}_q50": np.nan,
            f"amd_{kind}_q84": np.nan,
        }

        return pd.Series(out)

    # Compute the amd
    amd = compute_amd(mass_mc, eccen_mc, di_mc, sma_mc)

    if namd:
        out = {
            f"amd_{kind}_mc": amd,
            f"mass_{kind}_mc": mass_mc,
            f"sqrt_sma_{kind}_mc": np.sqrt(sma_mc),
        }

        return pd.Series(out)

    amd_p = np.percentile(amd.compressed(), [0.16, 0.5, 0.84])

    out = {
        f"amd_{kind}_q16": amd_p[0],
        f"amd_{kind}_q50": amd_p[1],
        f"amd_{kind}_q84": amd_p[2],
    }

    return pd.Series(out)
