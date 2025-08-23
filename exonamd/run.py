import os
import numpy as np
import pandas as pd
import swifter
import warnings
from loguru import logger

from exonamd.catalog import download_nasa_confirmed_planets
from exonamd.utils import ROOT
from exonamd.utils import fetch_aliases
from exonamd.utils import update_host
from exonamd.utils import update_planet
from exonamd.utils import check_name
from exonamd.solve import solve_values
from exonamd.interp import interp_eccentricity
from exonamd.interp import interp_mass
from exonamd.interp import interp_inclination
from exonamd.interp import interp_sma
from exonamd.solve import solve_relincl
from exonamd.interp import interp_trueobliq
from exonamd.utils import groupby_apply_merge
from exonamd.solve import solve_namd
from exonamd.solve import solve_namd_mc
from exonamd.plot import plot_sample_namd


warnings.filterwarnings("ignore", category=RuntimeWarning, module="pandas")
pd.options.display.max_columns = 20
pd.options.display.max_rows = 30
pd.options.mode.copy_on_write = True
swifter.set_defaults(
    npartitions=None,
    dask_threshold=1,
    scheduler="processes",
    progress_bar=True,
    progress_bar_desc=None,
    allow_dask_on_strings=False,
    force_parallel=False,
)


__all__ = [
    "create_db",
    "interp_db",
    "calc_namd",
]


def create_db(from_scratch=True, out_path=""):
    """
    Downloads the NASA Exoplanet Archive confirmed planets table, deals with aliases, computes missing values from simple relations, and stores the curated database.

    Parameters
    ----------
    from_scratch : bool, optional
        If True, downloads the entire table. If False, downloads only the rows newer than the latest update date in the current table. Defaults to True.
    out_path : str, optional
        The path where the curated database will be stored. If empty, the default path will be used.

    Returns
    -------
    df : pandas.DataFrame
        The curated database.
    """
    # Task 1: get the data
    df, df_old = download_nasa_confirmed_planets(
        min_sy_pnum=2,
        from_scratch=from_scratch,
    )

    # Task 2: deal with the aliases
    out_path = os.path.join(ROOT, "data", "aliases.pkl")
    aliases = fetch_aliases(
        hosts=df["hostname"].unique(),
        output_file=out_path,
    )

    logger.info("Updating host and planet names")
    df["hostname"] = df.swifter.apply(update_host, args=(aliases,), axis=1)
    df["pl_name"] = df.swifter.apply(update_planet, args=(aliases,), axis=1)
    logger.info("Names updated")

    logger.info("Checking consistency of planet names")
    name_ok = df.groupby("hostname")["pl_name"].apply(check_name)
    for hostname in name_ok[~name_ok].index:
        logger.error(f"Inconsistent planet names for {hostname}")
    logger.info("Consistency check done")

    # Task 3: compute missing values (if any) from simple relations
    logger.info("Computing missing values from simple relations")
    df[
        [
            "pl_orbsmax",
            "pl_ratdor",
            "st_rad",
            "pl_rade",
            "pl_ratror",
            "pl_orbper",
            "st_mass",
        ]
    ] = df.swifter.apply(solve_values, axis=1)
    logger.info("Missing values computed")

    # Task 4: store the curated database
    logger.debug("Dropping columns that are no longer needed")
    df.drop(
        columns=[
            "pl_ratdor",
            "st_rad",
            "pl_ratror",
            "pl_orbper",
            "st_mass",
        ],
        inplace=True,
    )
    logger.debug("Columns dropped")

    logger.info("Storing the curated database")
    if df_old is not None:
        df = pd.concat([df.copy(), df_old], ignore_index=True)
        df = df.drop_duplicates(keep="last")
        df.reset_index(drop=True)

    if out_path == "":
        out_path = os.path.join(ROOT, "data", "exo.csv")
    df.to_csv(out_path, index=False)
    logger.info(f"Database stored at {out_path}")

    return df


def interp_db(df: pd.DataFrame, out_path=""):
    """
    Curates and interpolates missing values in the NASA exoplanet database.

    This function reloads the database, thins it down by using the median of the values for each parameter of each planet, and then interpolates missing values of:
    - eccentricities
    - planetary masses
    - inclinations
    - semi-major axis uncertainties
    - true obliquity

    The function then stores the curated+interpolated database in a new file.

    Parameters
    ----------
    df : pd.DataFrame
        The database to be curated and interpolated. If None, the function will reload the database.
    out_path : str, optional
        The path where the curated+interpolated database will be stored. If empty, the default path will be used.

    Returns
    -------
    pd.DataFrame
        The curated+interpolated database.

    Raises
    ------
    ValueError
        If there are duplicated rows for a given hostname+planet name pair.

    Notes
    -----
    The output database will have the same columns as the input database, with the addition of the "flag" column which indicates whether the value was interpolated or not (flag=0).
    """
    # Task 1: reload database
    if df is None:
        logger.info("Reloading the database")
        df = pd.read_csv(os.path.join(ROOT, "data", "exo.csv"))
        logger.info("Database reloaded")

    # Task 2: input missing values (if any) by interpolation
    logger.info("Thinning down the data with nanmedian")
    NaT_idx = df[df["rowupdate"].isna()].index
    if len(NaT_idx) > 0:
        logger.warning(f"NaT values in rowupdate: {len(NaT_idx)}")
        df.loc[NaT_idx, "rowupdate"] = pd.to_datetime("1900-01-01")
    df["rowupdate"] = pd.to_datetime(df["rowupdate"], format="%Y-%m-%d")
    latest_update_idx = df.groupby("pl_name")["rowupdate"].idxmax()

    cols = df.columns.difference(["hostname", "pl_name", "default_flag", "rowupdate"])
    medians = df.groupby("pl_name")[cols].transform(np.nanmedian)
    df.loc[latest_update_idx, cols] = medians.loc[latest_update_idx]
    df = df.loc[latest_update_idx].drop(columns="default_flag")
    logger.info("Data thinned down")

    logger.info("Checking for duplicates")
    dp = df[df.duplicated(subset=["hostname", "pl_name"], keep=False)].sort_values(
        by=["hostname", "pl_name"]
    )

    if len(dp) > 0:
        logger.error(f"Duplicated rows for {dp['hostname'].unique()}")
        raise ValueError(f"Duplicated rows for {dp['hostname'].unique()}")
    logger.info("No duplicates found")

    logger.info("Instantiating the flags")
    df["flag"] = "0"
    logger.info("Flags instantiated")

    logger.info("Interpolating missing eccentricity values")
    df[
        [
            "pl_orbeccen",
            "pl_orbeccenerr1",
            "pl_orbeccenerr2",
            "flag",
        ]
    ] = df.swifter.apply(interp_eccentricity, axis=1)
    logger.info("Values interpolated")

    logger.info("Interpolating missing planetary mass values")
    df[
        [
            "pl_bmasse",
            "pl_bmasseerr1",
            "pl_bmasseerr2",
            "flag",
        ]
    ] = df.swifter.apply(interp_mass, axis=1)
    logger.info("Values interpolated")

    logger.debug("Dropping columns that are no longer needed")
    df.drop(columns=["pl_rade", "pl_radeerr1", "pl_radeerr2"], inplace=True)
    logger.debug("Columns dropped")

    logger.info(
        "Removing systems where at least one planet has no mass or semi-major axis"
    )
    mask = (
        df.groupby("hostname")[["pl_bmasse", "pl_orbsmax"]]
        .transform(lambda x: x.isnull().any())
        .any(axis=1)
    )
    rm_systems = df[mask]["hostname"].unique()
    logger.info(f"Removing {len(rm_systems)} systems: {rm_systems}")
    df = df[~mask]
    logger.info("Systems removed")

    logger.info("Interpolating missing values in inclinations")
    df[
        [
            "pl_orbincl",
            "pl_orbinclerr1",
            "pl_orbinclerr2",
            "flag",
        ]
    ] = df.swifter.apply(interp_inclination, args=(df,), axis=1)
    logger.info("Values interpolated")

    logger.info("Interpolating missing values in semi-major axis uncertainties")
    df[
        [
            "pl_orbsmaxerr1",
            "pl_orbsmaxerr2",
            "flag",
        ]
    ] = df.swifter.apply(interp_sma, axis=1)
    logger.info("Values interpolated")

    # Task 3: compute the parameters for the NAMD calculation
    logger.info("Computing the relative inclinations")
    df[
        [
            "pl_relincl",
            "pl_relinclerr1",
            "pl_relinclerr2",
        ]
    ] = df.swifter.apply(solve_relincl, args=(df,), axis=1)
    logger.info("Values computed")

    logger.info("Interpolating missing values in true obliquity")
    df[
        [
            "pl_trueobliq",
            "pl_trueobliqerr1",
            "pl_trueobliqerr2",
            "flag",
        ]
    ] = df.swifter.apply(interp_trueobliq, args=(df,), axis=1)
    logger.info("Values interpolated")

    # Task 4: store the curated+interpolated database
    logger.info("Storing the curated+interpolated database")

    if out_path == "":
        out_path = os.path.join(ROOT, "data", "exo_interp.csv")
    df.to_csv(out_path, index=False)
    logger.info(f"Database stored at {out_path}")

    return df


def calc_namd(
    df: pd.DataFrame,
    save=True,
    plot=False,
    core=True,
    filt=None,
    which=["rel", "abs"],
    threshold=100,
    use_trunc_normal=True,
    Npt=10000,
    out_path="",
):
    """
    Compute the NAMD for a given sample of planetary systems.

    Parameters
    ----------
    df : pd.DataFrame
        The input database. If None, it will be reloaded from the default location.
    save : bool, optional
        Whether to save the output database. Default is True.
    plot : bool, optional
        Whether to plot some diagnostic plots. Default is False.
    core : bool, optional
        Whether to select only the "core" sample. Default is True.
    filt : callable, optional
        The filter to obtain the "core" sample. Default is None.

    Notes
    -----
    The output database will have the same columns as the input database, with the addition of the NAMD parameters.

    The "core" sample is defined by default as the systems with all planets having a flag of either 0, 05+, 05-, or 05+-, i.e. nothing or only the obliquity has been interpolated. Nonetheless, a custom filter can be provided via the `filt` parameter.
    """
    # Task 1: reload database
    if df is None:
        logger.info("Reloading the database")
        df = pd.read_csv(os.path.join(ROOT, "data", "exo_interp.csv"))
        logger.info("Database reloaded")

    logger.debug("Dropping columns that are no longer needed")
    df = df.drop(columns=["pl_orbincl", "pl_orbinclerr1", "pl_orbinclerr2"])
    logger.debug("Columns dropped")

    # Task 2: compute the NAMD
    if not set(which).issubset({"rel", "abs"}):
        raise ValueError(
            "Invalid 'which' parameter. Must be a subset of {'rel', 'abs'}."
        )

    if "rel" in which:
        logger.info("Computing the relative NAMD")
        df = groupby_apply_merge(
            df,
            "hostname",
            solve_namd,
            kind="rel",
            allow_overwrite=True,
        )
        logger.info("Relative NAMD computed")

    if "abs" in which:
        logger.info("Computing the absolute NAMD")
        df = groupby_apply_merge(
            df,
            "hostname",
            solve_namd,
            kind="abs",
            allow_overwrite=True,
        )
        logger.info("Absolute NAMD computed")

    if plot and (which == ["rel", "abs"]):
        (
            df.groupby("hostname")[["namd_rel", "namd_abs"]]
            .transform("mean")
            .plot(
                kind="scatter",
                x="namd_rel",
                y="namd_abs",
                loglog=True,
                title="Full sample",
            )
        )

    if plot and ("rel" in which):
        (
            df.groupby("hostname")[["sy_pnum", "namd_rel"]]
            .mean()
            .reset_index()
            .plot(
                kind="scatter",
                x="sy_pnum",
                y="namd_rel",
                logy=True,
                title="Full sample",
            )
        )

    if plot and ("abs" in which):
        (
            df.groupby("hostname")[["sy_pnum", "namd_abs"]]
            .mean()
            .reset_index()
            .plot(
                kind="scatter",
                x="sy_pnum",
                y="namd_abs",
                logy=True,
                title="Full sample",
            )
        )

    if core and (filt is None):
        logger.info("Defining the core sample")
        core_flags = ["0", "05+", "05-", "05+-", "05d+-"]
        df = df.groupby("hostname").filter(lambda x: all(x["flag"].isin(core_flags)))
        logger.info("Core sample defined")
    elif core and (filt is not None):
        logger.info("Defining the core sample using the custom filter")
        df = df.groupby("hostname").filter(filt)

    if plot and ("rel" in which) and core:
        (
            df.groupby("hostname")[["sy_pnum", "namd_rel"]]
            .transform("mean")
            .plot(
                kind="scatter",
                x="sy_pnum",
                y="namd_rel",
                logy=True,
                title="Core sample",
            )
        )

    if plot and ("abs" in which) and core:
        (
            df.groupby("hostname")[["sy_pnum", "namd_abs"]]
            .transform("mean")
            .plot(
                kind="scatter",
                x="sy_pnum",
                y="namd_abs",
                logy=True,
                title="Core sample",
            )
        )

    # Task 3: compute the NAMD and associated confidence intervals
    if "rel" in which:
        logger.info("Computing the Monte Carlo relative NAMD")
        df = groupby_apply_merge(
            df,
            "hostname",
            solve_namd_mc,
            kind="rel",
            Npt=Npt,
            threshold=threshold,
            use_trunc_normal=use_trunc_normal,
            allow_overwrite=True,
        )
        logger.info("Relative NAMD computed")

    if "abs" in which:
        logger.info("Computing the Monte Carlo absolute NAMD")
        df = groupby_apply_merge(
            df,
            "hostname",
            solve_namd_mc,
            kind="abs",
            Npt=Npt,
            threshold=threshold,
            use_trunc_normal=use_trunc_normal,
            allow_overwrite=True,
        )
        logger.info("Absolute NAMD computed")

    if plot and (which == ["rel", "abs"]):
        (
            df.groupby("hostname")[["namd_rel_q50", "namd_abs_q50"]]
            .transform("mean")
            .plot(kind="scatter", x="namd_rel_q50", y="namd_abs_q50", loglog=True)
        )

    # Task 4: store the namd database
    if save:
        logger.info("Storing the NAMD database")
        if out_path == "":
            out_path = os.path.join(ROOT, "data", "exo_namd.csv")
        df.to_csv(out_path, index=False)
        logger.info(f"Database stored at {out_path}")

    return df


def run(from_scratch=True):
    df = create_db(from_scratch)

    df = interp_db(df)

    df = calc_namd(df, core=True)

    plot_sample_namd(df, title="Sample systems")
