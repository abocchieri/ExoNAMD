import os
import numpy as np
import pandas as pd
import swifter
import warnings

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


def create_db(from_scratch=True):

    # Task 1: get the data
    df, df_old = download_nasa_confirmed_planets(
        sy_pnum=1,
        from_scratch=from_scratch,
    )

    # Task 2: deal with the aliases
    aliases = fetch_aliases(df["hostname"].unique())

    df["hostname"] = df.swifter.apply(update_host, args=(aliases,), axis=1)
    df["pl_name"] = df.swifter.apply(update_planet, args=(aliases,), axis=1)

    name_ok = df.groupby("hostname")["pl_name"].apply(check_name)
    for hostname in name_ok[~name_ok].index:
        print(f"Inconsistent planet names for {hostname}")

    # Task 3: compute missing values (if any) from simple equations
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

    # Task 4: store the curated database
    if df_old is not None:
        df = pd.concat([df.copy(), df_old], ignore_index=True)
        df = df.drop_duplicates(keep="last")
        df.reset_index(drop=True)

    df.to_csv(os.path.join(ROOT, "data", "exo.csv"), index=False)

    return df


def interp_db(df: pd.DataFrame):

    # Task 1: reload database
    if not df:
        df = pd.read_csv(os.path.join(ROOT, "data", "exo.csv"))

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

    # Task 2: input missing values (if any) by interpolation
    cols = [
        c
        for c in df.columns
        if c not in ["hostname", "pl_name", "default_flag", "rowupdate"]
    ]

    medians = df.groupby("pl_name")[cols].transform(np.nanmedian)
    df.loc[df["default_flag"] == 1, cols] = medians.loc[df["default_flag"] == 1]
    df = df[df["default_flag"] == 1]
    df.drop(columns="default_flag", inplace=True)

    duplicates = df[
        df.duplicated(subset=["hostname", "pl_name"], keep=False)
    ].sort_values(by=["hostname", "pl_name"])

    if len(duplicates) > 0:
        raise ValueError("Duplicated entries found")

    df["flag"] = "0"

    df[
        [
            "pl_orbeccen",
            "pl_orbeccenerr1",
            "pl_orbeccenerr2",
            "flag",
        ]
    ] = df.swifter.apply(interp_eccentricity, axis=1)

    df[
        [
            "pl_bmasse",
            "pl_bmasseerr1",
            "pl_bmasseerr2",
            "flag",
        ]
    ] = df.swifter.apply(interp_mass, axis=1)

    df.drop(columns=["pl_rade", "pl_radeerr1", "pl_radeerr2"], inplace=True)

    mask = (
        df.groupby("hostname")[["pl_bmasse", "pl_orbsmax"]]
        .transform(lambda x: x.isnull().any())
        .any(axis=1)
    )
    df = df[~mask]

    df[
        [
            "pl_orbincl",
            "pl_orbinclerr1",
            "pl_orbinclerr2",
            "flag",
        ]
    ] = df.swifter.apply(interp_inclination, args=(df,), axis=1)

    df[
        [
            "pl_orbsmaxerr1",
            "pl_orbsmaxerr2",
            "flag",
        ]
    ] = df.swifter.apply(interp_sma, axis=1)

    # Task 3: compute the parameters for the NAMD calculation
    df[
        [
            "pl_relincl",
            "pl_relinclerr1",
            "pl_relinclerr2",
        ]
    ] = df.swifter.apply(solve_relincl, args=(df,), axis=1)

    df[
        [
            "pl_trueobliq",
            "pl_trueobliqerr1",
            "pl_trueobliqerr2",
            "flag",
        ]
    ] = df.swifter.apply(interp_trueobliq, args=(df,), axis=1)

    # Task 4: store the curated+interpolated database
    df.to_csv(os.path.join(ROOT, "data", "exo_interp.csv"), index=False)

    return df


def calc_namd(df: pd.DataFrame, plot=False, core=True):

    # Task 1: reload database
    if not df:
        df = pd.read_csv(os.path.join(ROOT, "data", "exo_interp.csv"))
    df.drop(columns=["pl_orbincl", "pl_orbinclerr1", "pl_orbinclerr2"], inplace=True)

    # Task 2: compute the NAMD
    df = groupby_apply_merge(
        df,
        "hostname",
        solve_namd,
        kind="rel",
        allow_overwrite=True,
    )

    df = groupby_apply_merge(
        df,
        "hostname",
        solve_namd,
        kind="abs",
        allow_overwrite=True,
    )

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
        core_flags = ["0", "05+", "05-", "05+-"]
        df = df.groupby("hostname").filter(lambda x: all(x["flag"].isin(core_flags)))

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

    df = groupby_apply_merge(
        df,
        "hostname",
        solve_namd_mc,
        kind="rel",
        Npt=Npt,
        threshold=threshold,
        allow_overwrite=True,
    )

    df = groupby_apply_merge(
        df,
        "hostname",
        solve_namd_mc,
        kind="abs",
        Npt=Npt,
        threshold=threshold,
        allow_overwrite=True,
    )

    if plot:
        (
            df.groupby("hostname")[["namd_rel_q50", "namd_abs_q50"]]
            .transform("mean")
            .plot(kind="scatter", x="namd_rel_q50", y="namd_abs_q50", loglog=True)
        )

    # Task 4: store the namd database
    df.to_csv(os.path.join(ROOT, "data", "exo_namd.csv"), index=False)


def plot_sample_namd(df: pd.DataFrame):

    # Task 1: reload database
    if not df:
        df = pd.read_csv(os.path.join(ROOT, "data", "exo_namd.csv"))

    # Task 2: plot the sample NAMD
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


def plot_host_namd(
    df: pd.DataFrame, hostname: str, Npt: int = 100000, threshold: int = 1000
):

    # Task 1: reload database
    if not df:
        df = pd.read_csv(os.path.join(ROOT, "data", "exo_namd.csv"))

    # Task 2: sample the NAMD for a given host
    host = df[df["hostname"] == hostname]
    retval = solve_namd_mc(
        host=host,
        kind="rel",
        Npt=Npt,
        threshold=threshold,
        full=True,
    )

    # Task 2: plot the NAMD for a given host
    simple_plot(
        df=retval,
        kind="rel",
        title=hostname,
        which="namd",
        scale="log",
    )


def run(from_scratch=True):

    df = create_db(from_scratch)

    df = interp_db(df)

    df = calc_namd(df)

    plot_sample_namd(df)
