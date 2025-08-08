import numpy as np
import pandas as pd
from spright import RMRelation
from loguru import logger

from exonamd.utils import get_value


rmr = RMRelation()


@logger.catch
def interp_eccentricity(row):
    """
    Interpolates missing eccentricity values using the relation found in
    ...

    Also, sets flags for the interpolation.

    Parameters
    ----------
    row: pandas.Series
        The planet row in the sytem table.

    Returns
    -------
    pandas.Series
        A pandas Series containing the interpolated eccentricity, the associated uncertainties and flag.
    """
    ecc = get_value(row["pl_orbeccen"])
    eccerr1 = get_value(row["pl_orbeccenerr1"])
    eccerr2 = get_value(row["pl_orbeccenerr2"])
    flag = get_value(row["flag"])
    sy_pnum = get_value(row["sy_pnum"])

    if np.isnan(ecc):
        ecc = 0.63 * sy_pnum ** (-1.02)
        eccerr1 = 0.0
        eccerr2 = 0.0
        flag += "1+-"

    if np.isnan(eccerr1):
        eccerr1 = 0.0
        flag += "1+"

    if np.isnan(eccerr2):
        eccerr2 = 0.0
        flag += "1-"

    out = {
        "pl_orbeccen": ecc,
        "pl_orbeccenerr1": eccerr1,
        "pl_orbeccenerr2": eccerr2,
        "flag": flag,
    }

    return pd.Series(out)


@logger.catch
def interp_mass(row, min_radius=0.5, max_radius=6.0):
    """
    Interpolate missing mass values by using the mass-radius relation implemented in the spright package, which is described in https://spright.readthedocs.io/en/latest/.

    Parameters
    ----------
    row : pandas.Series
        The planet row in the sytem table.
    min_radius : float, optional
        Minimum radius in units of R_earth below which the interpolation is not performed, by default 0.5.
    max_radius : float, optional
        Maximum radius in units of R_earth  above which the interpolation is not performed, by default 6.0.

    Returns
    -------
    pandas.Series
        A pandas Series containing the interpolated mass, the associated uncertainties and flag.
    """
    mass = get_value(row["pl_bmasse"])
    masserr1 = get_value(row["pl_bmasseerr1"])
    masserr2 = get_value(row["pl_bmasseerr2"])
    radius = get_value(row["pl_rade"])
    radiuserr1 = get_value(row["pl_radeerr1"])
    radiuserr2 = get_value(row["pl_radeerr2"])
    flag = get_value(row["flag"])

    if (
        np.isnan(mass)
        and not np.isnan(radius)
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

    out = {
        "pl_bmasse": mass,
        "pl_bmasseerr1": masserr1,
        "pl_bmasseerr2": masserr2,
        "flag": flag,
    }

    return pd.Series(out)


@logger.catch
def interp_sma(row):
    """
    Interpolate missing semi-major axis uncertainties by setting them to zero and adding the flag.

    Parameters
    ----------
    row : pandas.Series
        The planet row in the sytem table.

    Returns
    -------
    pandas.Series
        A pandas Series containing the interpolated semi-major axis uncertainties and the associated flag.
    """
    smaerr1 = get_value(row["pl_orbsmaxerr1"])
    smaerr2 = get_value(row["pl_orbsmaxerr2"])
    flag = get_value(row["flag"])

    if np.isnan(smaerr1):
        smaerr1 = 0.0
        flag += "4+"

    if np.isnan(smaerr2):
        smaerr2 = 0.0
        flag += "4-"

    out = {
        "pl_orbsmaxerr1": smaerr1,
        "pl_orbsmaxerr2": smaerr2,
        "flag": flag,
    }

    return pd.Series(out)


def interpolate_angle(row, df, value_type):
    """
    Interpolate missing values for inclination or obliquity by using the values from the most massive planet in the same system.

    Parameters
    ----------
    row : pandas.Series
        The planet row in the sytem table.
    df : pandas.DataFrame
        The system table.
    value_type : str
        Either "inclination" or "obliquity".

    Returns
    -------
    pandas.Series
        A pandas Series containing the interpolated value, the associated uncertainties and flag.
    """
    hostname = get_value(row["hostname"])
    flag = get_value(row["flag"])

    if value_type == "inclination":
        value_col = "pl_orbincl"
        err1_col = "pl_orbinclerr1"
        err2_col = "pl_orbinclerr2"
        backup_value = 90.0
        backup_err1 = 0.0
        backup_err2 = 0.0
        flag_suffix = "3"
    elif value_type == "obliquity":
        value_col = "pl_trueobliq"
        err1_col = "pl_trueobliqerr1"
        err2_col = "pl_trueobliqerr2"
        backup_value = get_value(row["pl_relincl"])
        backup_err1 = get_value(row["pl_relinclerr1"])
        backup_err2 = get_value(row["pl_relinclerr2"])
        flag_suffix = "5"
    else:
        raise ValueError("Invalid value_type provided")

    value = get_value(row[value_col])
    err1 = get_value(row[err1_col])
    err2 = get_value(row[err2_col])

    host = df[df["hostname"] == hostname]
    max_mass_idx = host["pl_bmasse"].idxmax()
    valuenan = host[host[value_col].isnull()]

    # if value is not nan, check the errors
    # if they are nan, set them to 0
    if not np.isnan(value):
        if np.isnan(err1):
            err1 = 0.0
            flag += f"{flag_suffix}+"
        if np.isnan(err2):
            err2 = 0.0
            flag += f"{flag_suffix}-"

    elif len(valuenan) == len(host):
        # if all planets have value nan
        value = backup_value
        err1 = backup_err1
        err2 = backup_err2
        flag += f"{flag_suffix}d+-"

    else:
        if max_mass_idx in valuenan.index:
            # if the most massive planet has value nan
            # we go to the next most massive planet
            nonan_values = host[host[value_col].notnull()]
            max_mass_idx = nonan_values["pl_bmasse"].idxmax()

        max_mass = host.loc[max_mass_idx, [value_col, err1_col, err2_col]]
        value = max_mass[value_col]
        err1 = max_mass[err1_col] if not np.isnan(max_mass[err1_col]) else 0.0
        err2 = max_mass[err2_col] if not np.isnan(max_mass[err2_col]) else 0.0
        flag += f"{flag_suffix}+-"

    out = {
        value_col: value,
        err1_col: err1,
        err2_col: err2,
        "flag": flag,
    }

    return pd.Series(out)


@logger.catch
def interp_inclination(row, df):
    """
    Wrapper of the function **interpolate_angle** to interpolate inclination by using the values from the most massive planet in the same system.

    Parameters
    ----------
    row : pandas.Series
        The planet row in the sytem table.
    df : pandas.DataFrame
        The system table.

    Returns
    -------
    pandas.Series
        A pandas Series containing the interpolated value, the associated uncertainties and flag.
    """
    return interpolate_angle(row, df, "inclination")


@logger.catch
def interp_trueobliq(row, df):
    """
    Wrapper of the function **interpolate_angle** to interpolate obliquity by using the values from the most massive planet in the same system.

    Parameters
    ----------
    row : pandas.Series
        The planet row in the sytem table.
    df : pandas.DataFrame
        The system table.

    Returns
    -------
    pandas.Series
        A pandas Series containing the interpolated value, the associated uncertainties and flag.
    """
    return interpolate_angle(row, df, "obliquity")
