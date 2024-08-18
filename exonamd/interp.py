import numpy as np
import pandas as pd

# import modin.pandas as pd
from spright import RMRelation


rmr = RMRelation()


def interp_eccentricity(row):
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

    out = {
        "pl_orbeccen": ecc,
        "pl_orbeccenerr1": eccerr1,
        "pl_orbeccenerr2": eccerr2,
        "flag": flag,
    }

    return pd.Series(out)


def interp_mass(row, min_radius=0.5, max_radius=6.0):
    mass = row["pl_bmasse"]
    masserr1 = row["pl_bmasseerr1"]
    masserr2 = row["pl_bmasseerr2"]
    radius = row["pl_rade"]
    radiuserr1 = row["pl_radeerr1"]
    radiuserr2 = row["pl_radeerr2"]
    flag = row["flag"]

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


def interp_inclination(row, df):
    hostname = row["hostname"]
    incl = row["pl_orbincl"]
    inclerr1 = row["pl_orbinclerr1"]
    inclerr2 = row["pl_orbinclerr2"]
    flag = row["flag"]

    host = df[df["hostname"] == hostname]
    max_mass = host["pl_bmasse"].idxmax()
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
        # if all planets have inclination nan
        incl = 90.0
        inclerr1 = 0.0
        inclerr2 = 0.0
        flag += "3+-"

    else:
        if max_mass in inclnan.index:
            # if the most massive planet has inclination nan
            # we go to the next most massive planet
            inclnotnan = host[host["pl_orbincl"].notnull()]
            max_mass = inclnotnan["pl_bmasse"].idxmax()

        max_mass = host.loc[
            max_mass, ["pl_orbincl", "pl_orbinclerr1", "pl_orbinclerr2"]
        ]

        incl = max_mass["pl_orbincl"]
        inclerr1 = max_mass["pl_orbinclerr1"]
        inclerr1 = inclerr1 if not np.isnan(inclerr1) else 0.0
        inclerr2 = max_mass["pl_orbinclerr2"]
        inclerr2 = inclerr2 if not np.isnan(inclerr2) else 0.0
        flag += "3+-"

    out = {
        "pl_orbincl": incl,
        "pl_orbinclerr1": inclerr1,
        "pl_orbinclerr2": inclerr2,
        "flag": flag,
    }

    return pd.Series(out)


def interp_sma(row):
    smaerr1 = row["pl_orbsmaxerr1"]
    smaerr2 = row["pl_orbsmaxerr2"]
    flag = row["flag"]

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


def interp_trueobliq(row, df):
    hostname = row["hostname"]
    obliq = row["pl_trueobliq"]
    obliqerr1 = row["pl_trueobliqerr1"]
    obliqerr2 = row["pl_trueobliqerr2"]
    relincl = row["pl_relincl"]
    relinclerr1 = row["pl_relinclerr1"]
    relinclerr2 = row["pl_relinclerr2"]
    flag = row["flag"]

    host = df[df["hostname"] == hostname]
    max_mass = host["pl_bmasse"].idxmax()
    obliqnan = host[host["pl_trueobliq"].isnull()]

    # if obliquity is not nan, check the errors
    # if they are nan, set them to 0
    if not np.isnan(obliq):
        if np.isnan(obliqerr1):
            obliqerr1 = 0.0
            flag += "5+"
        if np.isnan(obliqerr2):
            obliqerr2 = 0.0
            flag += "5-"

    elif len(obliqnan) == len(host):
        # if all planets have obliquity nan
        obliq = relincl
        obliqerr1 = relinclerr1
        obliqerr2 = relinclerr2
        flag += "5+-"

    else:
        if max_mass in obliqnan.index:
            # if the most massive planet has obliquity nan
            # we go to the next most massive planet
            obliqnotnan = host[host["pl_trueobliq"].notnull()]
            max_mass = obliqnotnan["pl_bmasse"].idxmax()

        max_mass = host.loc[
            max_mass, ["pl_trueobliq", "pl_trueobliqerr1", "pl_trueobliqerr2"]
        ]
        obliq = max_mass["pl_trueobliq"]
        obliqerr1 = max_mass["pl_trueobliqerr1"]
        obliqerr1 = obliqerr1 if not np.isnan(obliqerr1) else 0.0
        obliqerr2 = max_mass["pl_trueobliqerr2"]
        obliqerr2 = obliqerr2 if not np.isnan(obliqerr2) else 0.0
        flag += "5+-"

    out = {
        "pl_trueobliq": obliq,
        "pl_trueobliqerr1": obliqerr1,
        "pl_trueobliqerr2": obliqerr2,
        "flag": flag,
    }

    return pd.Series(out)
