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
from exonamd.plot import simple_plot
from exonamd.plot import pop_plot


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
    "plot_sample_namd",
    "plot_host_namd",
]


def create_db(from_scratch=True):
    # Task 1: get the data
    df, df_old = download_nasa_confirmed_planets(
        sy_pnum=1,
        from_scratch=from_scratch,
    )

    # Task 2: deal with the aliases
    aliases = fetch_aliases(df["hostname"].unique())

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

    out_path = os.path.join(ROOT, "data", "exo.csv")
    df.to_csv(out_path, index=False)
    logger.info(f"Database stored at {out_path}")

    return df


def interp_db(df: pd.DataFrame):
    # Task 1: reload database
    if not df:
        logger.info("Reloading the database")
        df = pd.read_csv(os.path.join(ROOT, "data", "exo.csv"))
        logger.info("Database reloaded")

    # Task 2: input missing values (if any) by interpolation
    logger.info("Thinning down the data with nanmedian")
    cols = df.columns.difference(["hostname", "pl_name", "default_flag", "rowupdate"])
    medians = df.groupby("pl_name")[cols].transform(np.nanmedian)
    df.loc[df["default_flag"] == 1, cols] = medians.loc[df["default_flag"] == 1]
    df = df[df["default_flag"] == 1].drop(columns="default_flag")
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
    out_path = os.path.join(ROOT, "data", "exo_interp.csv")
    df.to_csv(out_path, index=False)
    logger.info(f"Database stored at {out_path}")

    return df


def calc_namd(df: pd.DataFrame, plot=False, core=True):
    # Task 1: reload database
    if not df:
        logger.info("Reloading the database")
        df = pd.read_csv(os.path.join(ROOT, "data", "exo_interp.csv"))
        logger.info("Database reloaded")

    logger.debug("Dropping columns that are no longer needed")
    df.drop(columns=["pl_orbincl", "pl_orbinclerr1", "pl_orbinclerr2"], inplace=True)
    logger.debug("Columns dropped")

    # Task 2: compute the NAMD
    logger.info("Computing the relative NAMD")
    df = groupby_apply_merge(
        df,
        "hostname",
        solve_namd,
        kind="rel",
        allow_overwrite=True,
    )
    logger.info("Relative NAMD computed")

    logger.info("Computing the absolute NAMD")
    df = groupby_apply_merge(
        df,
        "hostname",
        solve_namd,
        kind="abs",
        allow_overwrite=True,
    )
    logger.info("Absolute NAMD computed")

    if plot:
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

    if core:
        logger.info("Defining the core sample")
        core_flags = ["0", "05+", "05-", "05+-"]
        df = df.groupby("hostname").filter(lambda x: all(x["flag"].isin(core_flags)))
        logger.info("Core sample defined")

    if plot:
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

    # Task 3: compute the NAMD and associated confidence intervals
    Npt = 200000
    threshold = 100

    logger.info("Computing the Monte Carlo relative NAMD")
    df = groupby_apply_merge(
        df,
        "hostname",
        solve_namd_mc,
        kind="rel",
        Npt=Npt,
        threshold=threshold,
        allow_overwrite=True,
    )
    logger.info("Relative NAMD computed")

    logger.info("Computing the Monte Carlo absolute NAMD")
    df = groupby_apply_merge(
        df,
        "hostname",
        solve_namd_mc,
        kind="abs",
        Npt=Npt,
        threshold=threshold,
        allow_overwrite=True,
    )
    logger.info("Absolute NAMD computed")

    if plot:
        (
            df.groupby("hostname")[["namd_rel_q50", "namd_abs_q50"]]
            .transform("mean")
            .plot(kind="scatter", x="namd_rel_q50", y="namd_abs_q50", loglog=True)
        )

    # Task 4: store the namd database
    logger.info("Storing the NAMD database")
    out_path = os.path.join(ROOT, "data", "exo_namd.csv")
    df.to_csv(out_path, index=False)
    logger.info(f"Database stored at {out_path}")

    return df


def plot_sample_namd(df: pd.DataFrame):
    # Task 1: reload database
    if not df:
        logger.info("Reloading the database")
        df = pd.read_csv(os.path.join(ROOT, "data", "exo_namd.csv"))
        logger.info("Database reloaded")

    # Task 2: plot the sample NAMD
    logger.info("Plotting the NAMD vs. multiplicity")
    pop_plot(
        df=df.groupby("hostname").apply(
            lambda g: g.select_dtypes(exclude=["object"]).mean(),
            include_groups=False,
        ),
        kind="rel",
        title="Core planets",
        which="namd",
        yscale="log",
        xoffs=True,
    )
    logger.info("Plot done")


def plot_host_namd(
    df: pd.DataFrame,
    hostname: str,
    which: str,
    Npt: int = 100000,
    threshold: int = 1000,
):
    # Task 1: reload database
    if not df:
        logger.info("Reloading the database")
        df = pd.read_csv(os.path.join(ROOT, "data", "exo_namd.csv"))
        logger.info("Database reloaded")

    # Task 2: sample the NAMD for a given host
    logger.info(f"Selecting the host: {hostname}")
    host = df[df["hostname"] == hostname]
    logger.info("Host selected")

    logger.info("Computing the Monte Carlo NAMD")
    retval = solve_namd_mc(
        host=host,
        kind="f{which}",
        Npt=Npt,
        threshold=threshold,
        full=True,
    )

    # Task 2: plot the NAMD for a given host
    simple_plot(
        df=retval,
        kind="f{which}",
        title=hostname,
        which="namd",
        scale="log",
    )
    logger.info("Plot done")


def run(from_scratch=True):
    df = create_db(from_scratch)

    df = interp_db(df)

    df = calc_namd(df)

    plot_sample_namd(df)
