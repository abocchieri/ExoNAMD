import numpy as np
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

    return ecc, eccerr1, eccerr2, flag


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

    return mass, masserr1, masserr2, flag


def interp_inclination(row, df):

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
        mass_max_data = host.loc[mass_max, ["pl_orbincl", "pl_orbinclerr1", "pl_orbinclerr2"]]

        incl = mass_max_data["pl_orbincl"]
        inclerr1 = 0.0 if np.isnan(mass_max_data["pl_orbinclerr1"]) else mass_max_data["pl_orbinclerr1"]
        inclerr2 = 0.0 if np.isnan(mass_max_data["pl_orbinclerr2"]) else mass_max_data["pl_orbinclerr2"]
        flag += "3+-"

    return incl, inclerr1, inclerr2, flag


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

    return smaerr1, smaerr2, flag