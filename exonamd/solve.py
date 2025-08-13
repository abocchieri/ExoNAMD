import numpy as np
import astropy.constants as cc
import astropy.units as u
import pandas as pd
from loguru import logger

from exonamd.utils import get_value
from exonamd.utils import sample_trunc_normal
from exonamd.core import compute_amdk
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
    missing = np.isnan(sma) + np.isnan(rstar) + np.isnan(ars)
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


@logger.catch
def solve_values(row):
    sma = get_value(row["pl_orbsmax"])
    ars = get_value(row["pl_ratdor"])
    rstar = get_value(row["st_rad"])
    rplanet = get_value(row["pl_rade"])
    rprs = get_value(row["pl_ratror"])
    period = get_value(row["pl_orbper"])
    mstar = get_value(row["st_mass"])

    # Rank groups
    a_rs_ = np.isnan(sma) + np.isnan(ars) + np.isnan(rstar)
    rprs_ = np.isnan(rplanet) + np.isnan(rprs) + np.isnan(rstar)
    a_period_ = np.isnan(period) + np.isnan(sma) + np.isnan(mstar)
    solve_order = np.argsort([a_rs_, rprs_, a_period_])
    for i in solve_order:
        if i == 0:
            logger.trace(
                "Solving semi-major axis -- stellar radius system of equations."
            )
            solution = solve_a_rs(sma, rstar, ars)
            sma, rstar, ars = solution
        elif i == 1:
            logger.trace("Solving planet radius -- stellar radius system of equations.")
            solution = solve_rprs(rplanet, rstar, rprs)
            rplanet, rstar, rprs = solution
        elif i == 2:
            logger.trace("Solving period-sma-mstar system of equations.")
            solution = solve_a_period(period, sma, mstar)
            period, sma, mstar = solution

    out = {
        "pl_orbsmax": sma,
        "pl_ratdor": ars,
        "st_rad": rstar,
        "pl_rade": rplanet,
        "pl_ratror": rprs,
        "pl_orbper": period,
        "st_mass": mstar,
    }

    return pd.Series(out)


@logger.catch
def solve_relincl(row, df):
    """
    Computes the relative inclination of a planet with respect to the most massive planet in the same system.

    Parameters
    ----------
    row: pandas.Series
        A row from the planet table.
    df: pandas.DataFrame
        The planet table.

    Returns
    -------
    pandas.Series
        A pandas Series containing the relative inclination and the associated uncertainties.
    """
    hostname = get_value(row["hostname"])
    incl = get_value(row["pl_orbincl"])
    inclerr1 = get_value(row["pl_orbinclerr1"])
    inclerr2 = get_value(row["pl_orbinclerr2"])

    host = df[df["hostname"] == hostname]
    max_mass = host["pl_bmasse"].idxmax()
    max_mass = host.loc[max_mass, ["pl_orbincl", "pl_orbinclerr1", "pl_orbinclerr2"]]

    relincl = max_mass["pl_orbincl"] - incl
    relinclerr1 = np.sqrt(inclerr1**2 + max_mass["pl_orbinclerr1"] ** 2)
    relinclerr2 = -np.sqrt(inclerr2**2 + max_mass["pl_orbinclerr2"] ** 2)

    out = {
        "pl_relincl": relincl,
        "pl_relinclerr1": relinclerr1,
        "pl_relinclerr2": relinclerr2,
    }

    return pd.Series(out)


def solve_amdk(row, kind: str):
    """
    Wrapper of the function **compute_amdk** to compute the angular momentum deficit (AMD) for a planet (or planets) in a system.

    Parameters
    ----------
    row: pandas.Series
        A row from the planet table.
    kind: str
        Which type of AMD to compute. One of 'rel' (relative, using the relative inclination) or 'abs' (absolute, using the obliquity).

    Returns
    -------
    pandas.Series
        A pandas Series containing the angular momentum deficit and the associated mass and sqrt of the semi-major axis.
    """
    mass = get_value(row["pl_bmasse"])
    eccen = get_value(row["pl_orbeccen"])
    di_ = {
        "rel": get_value(row["pl_relincl"]),
        "abs": get_value(row["pl_trueobliq"])
        / 2.0,  # divided by 2 to ensure the normalization of the absolute NAMD is between 0 and 1
    }
    di = di_[kind]
    sma = get_value(row["pl_orbsmax"])

    amdk = compute_amdk(mass, eccen, di, sma)

    out = {
        f"amdk_{kind}": amdk,
        "mass": mass,
        "sqrt_sma": np.sqrt(sma),
    }

    return pd.Series(out)


@logger.catch
def solve_namd(host, kind: str):
    """
    Wrapper of the functions **solve_amdk** and **compute_namd** to compute the normalized angular momentum deficit (NAMD) for a given system.

    Parameters
    ----------
    host: pandas.DataFrame
        A DataFrame containing the planet table.
    kind: str
        Which type of NAMD to compute. One of 'rel' (relative, using the relative inclination) or 'abs' (absolute, using the obliquity).

    Returns
    -------
    pandas.Series
        A pandas Series containing the normalized angular momentum deficit.
    """
    retval = host.apply(solve_amdk, args=(kind,), axis=1)

    amdk = retval[f"amdk_{kind}"]
    mass = retval["mass"]
    sqrt_sma = retval["sqrt_sma"]

    namd = compute_namd(amdk, mass, sqrt_sma)

    out = {f"namd_{kind}": namd}
    return pd.Series(out)


def solve_amdk_mc(row, kind, Npt, threshold, use_trunc_normal=True):
    """
    Compute the absolute or relative angular momentum deficit (AMD) for a given planet using a Monte Carlo approach.

    Parameters
    ----------
    row: pandas.Series
        A row from the planet table.
    kind: str
        Which type of AMD to compute. One of 'rel' (relative, using the relative inclination) or 'abs' (absolute, using the obliquity).
    Npt: int
        Number of Monte Carlo samples.
    threshold: int
        Minimum number of valid samples required.
    use_trunc_normal: bool
        If True, use a truncated normal distribution to sample the parameters. If False, use a normal distribution with rejection sampling.

    Returns
    -------
    pandas.Series
        A pandas Series containing the absolute or relative AMD, the associated mass and sqrt of the semi-major axis.
    """
    mass = get_value(row["pl_bmasse"])
    masserr1 = get_value(row["pl_bmasseerr1"])
    masserr2 = get_value(row["pl_bmasseerr2"])
    eccen = get_value(row["pl_orbeccen"])
    eccenerr1 = get_value(row["pl_orbeccenerr1"])
    eccenerr2 = get_value(row["pl_orbeccenerr2"])

    di_ = {
        "rel": np.abs(get_value(row["pl_relincl"])),  # sign is not important as this is the argument of the cosine
        "relerr1": get_value(row["pl_relinclerr1"]),
        "relerr2": get_value(row["pl_relinclerr2"]),
        "abs": get_value(row["pl_trueobliq"]),
        "abserr1": get_value(row["pl_trueobliqerr1"]),
        "abserr2": get_value(row["pl_trueobliqerr2"]),
    }

    di = di_[kind]
    dierr1 = di_[f"{kind}err1"]
    dierr2 = di_[f"{kind}err2"]
    di_upper = 180.0 if kind == "abs" else 90.0

    sma = get_value(row["pl_orbsmax"])
    smaerr1 = get_value(row["pl_orbsmaxerr1"])
    smaerr2 = get_value(row["pl_orbsmaxerr2"])

    # Sample the parameters
    if use_trunc_normal:
        mass_mc = sample_trunc_normal(
            mu=mass,
            sigma=0.5 * (masserr1 - masserr2),
            lower=0.0,
            upper=np.inf,
            n=Npt,
        )
        eccen_mc = sample_trunc_normal(
            mu=eccen,
            sigma=0.5 * (eccenerr1 - eccenerr2),
            lower=0.0,
            upper=1.0,
            n=Npt,
        )
        di_mc = sample_trunc_normal(
            mu=di,
            sigma=0.5 * (dierr1 - dierr2),
            lower=0.0,
            upper=di_upper,
            n=Npt,
        )
        sma_mc = sample_trunc_normal(
            mu=sma,
            sigma=0.5 * (smaerr1 - smaerr2),
            lower=0.0,
            upper=np.inf,
            n=Npt,
        )

    else:
        mass_mc = np.random.normal(mass, 0.5 * (masserr1 - masserr2), Npt)
        eccen_mc = np.random.normal(eccen, 0.5 * (eccenerr1 - eccenerr2), Npt)
        di_mc = np.random.normal(di, 0.5 * (dierr1 - dierr2), Npt)
        sma_mc = np.random.normal(sma, 0.5 * (smaerr1 - smaerr2), Npt)

        # Mask unphysical values
        mask = (
            (mass_mc > 0)
            & (eccen_mc >= 0)
            & (eccen_mc < 1)
            & (di_mc >= 0.0)
            & (di_mc < di_upper)
            & (sma_mc > 0)
        )

        # Check number of valid samples
        if np.sum(mask) < threshold:
            out = {
                f"amdk_{kind}_mc": np.nan,
                "mass_mc": np.nan,
                "sqrt_sma_mc": np.nan,
            }

            return pd.Series(out)

        mass_mc = np.ma.MaskedArray(mass_mc, mask=~mask)
        eccen_mc = np.ma.MaskedArray(eccen_mc, mask=~mask)
        di_mc = np.ma.MaskedArray(di_mc, mask=~mask)
        sma_mc = np.ma.MaskedArray(sma_mc, mask=~mask)

    if kind == "abs":
        di_mc /= 2.0  # divided by 2 to ensure the normalization of the absolute NAMD is between 0 and 1

    # Compute the amdk
    amdk = compute_amdk(mass_mc, eccen_mc, di_mc, sma_mc)

    out = {
        f"amdk_{kind}_mc": amdk,
        "mass_mc": mass_mc,
        "sqrt_sma_mc": np.sqrt(sma_mc),
    }

    return pd.Series(out)


@logger.catch
def solve_namd_mc(host, kind, Npt, threshold, use_trunc_normal, full=False):
    """
    Wrapper of the functions **solve_amdk_mc** and **compute_namd** to compute the normalized angular momentum deficit (NAMD) for a given system using a Monte Carlo approach.

    Parameters
    ----------
    host: pandas.DataFrame
        A DataFrame containing the system table.
    kind: str
        Which type of NAMD to compute. One of 'rel' (relative, using the relative inclination) or 'abs' (absolute, using the obliquity).
    Npt: int
        Number of Monte Carlo samples.
    threshold: int
        Minimum number of valid samples required.
    use_trunc_normal: bool
        If True, use a truncated normal distribution to sample the parameters. If False, use a normal distribution with rejection sampling.
    full: bool
        If True, return the full array of NAMD values. Otherwise, return only the 16th, 50th and 84th percentiles.

    Returns
    -------
    pandas.Series
        A pandas Series containing the normalized angular momentum deficit results.
    """
    retval = host.apply(
        solve_amdk_mc, args=(kind, Npt, threshold, use_trunc_normal), axis=1
    )

    amdk = retval[f"amdk_{kind}_mc"]
    mass = retval["mass_mc"]
    sqrt_sma = retval["sqrt_sma_mc"]

    namd = compute_namd(amdk, mass, sqrt_sma)
    if isinstance(namd, np.ma.MaskedArray):
        namd = namd.compressed()

    if len(namd) < threshold:
        out = {
            f"namd_{kind}_mc": np.nan,
            f"namd_{kind}_q16": np.nan,
            f"namd_{kind}_q50": np.nan,
            f"namd_{kind}_q84": np.nan,
        }

        return pd.Series(out)

    namd_p = np.percentile(namd, [16, 50, 84])

    out = {
        f"namd_{kind}_mc": namd if full else np.nan,
        f"namd_{kind}_q16": namd_p[0],
        f"namd_{kind}_q50": namd_p[1],
        f"namd_{kind}_q84": namd_p[2],
    }

    return pd.Series(out)
