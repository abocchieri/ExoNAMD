import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from loguru import logger

from exonamd.utils import ROOT
from exonamd.solve import solve_namd_mc


@logger.catch
def simple_plot(
    df,
    kind,
    title="",
    which="namd",
    ylabel="Frequency",
    scale="linear",
    bins=50,
    outpath=None,
):
    samples = df[f"{which}_{kind}_mc"]
    q50 = df[f"{which}_{kind}_q50"]
    q16 = df[f"{which}_{kind}_q16"]
    q84 = df[f"{which}_{kind}_q84"]

    xlabel = rf"{which.upper()}$_{kind[0].upper()}$"

    if scale == "log":
        samples = np.log10(samples)
        q16, q50, q84 = np.percentile(samples, [16, 50, 84])
        xlabel = rf"$\log_{{10}}$ {xlabel}"

    errup = q84 - q50
    errdown = q50 - q16
    title = f"{title}: " + rf"${q50:.2f}^{{+{errup:.2f}}}_{{-{errdown:.2f}}}$"

    plt.figure()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.hist(
        samples,
        bins=bins,
        histtype="step",
        weights=np.ones_like(samples) / len(samples),
    )
    plt.grid(which="both", linestyle="--", alpha=0.5)
    plt.vlines(
        [q16, q50, q84],
        0,
        plt.ylim()[1],
        color=["red", "black", "red"],
        linestyles="dashed",
    )
    if outpath:
        Path(outpath).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(outpath, bbox_inches="tight", dpi=300, format="pdf")

    plt.show()


@logger.catch
def pop_plot(df, kind, title="", which="namd", yscale="log", xoffs=False, outpath=None):
    # Plot the values vs multiplicity and color by their relative uncertainty

    df = df.sort_values(by="sy_pnum")
    sy_pnum = df["sy_pnum"]
    q50 = df[f"{which}_{kind}_q50"]
    q16 = df[f"{which}_{kind}_q16"]
    q84 = df[f"{which}_{kind}_q84"]

    nanidx = q50.isnull()
    sy_pnum = sy_pnum[~nanidx]
    q50 = q50[~nanidx]
    q16 = q16[~nanidx]
    q84 = q84[~nanidx]

    errup = q84 - q50
    errdown = q50 - q16
    iq = (q84 - q16) / 2
    sigma_rel = iq / q50

    ylabel = rf"{which.upper()}$_{kind[0].upper()}$"

    coeffs = np.polyfit(sy_pnum, q50, 1)
    line = np.polyval(coeffs, np.array(list(set(sy_pnum))))
    if yscale == "log":
        coeffs = np.polyfit(sy_pnum, np.log10(q50), 1)
        line = 10 ** np.polyval(coeffs, np.array(list(set(sy_pnum))))

    plt.figure()
    plt.plot(
        np.array(list(set(sy_pnum))),
        line,
        "k--",
        alpha=0.5,
        lw=1.5,
        zorder=10,
    )

    if xoffs:
        M = set(sy_pnum)
        n_list = []
        for m in M:
            idx = sy_pnum == m
            n_list.append(idx.sum())
        n_list = np.array(n_list)
        xoffs = 0.3 * n_list / n_list.max()
        for i, m in enumerate(M):
            idx = sy_pnum == m
            sy_pnum[idx] += np.linspace(-xoffs[i], xoffs[i], idx.sum())

    plt.errorbar(
        sy_pnum,
        q50,
        yerr=[errdown, errup],
        fmt="none",
        c="k",
        alpha=0.5,
        lw=0.5,
        capsize=2,
    )
    s = plt.scatter(
        sy_pnum,
        q50,
        c=sigma_rel,
        cmap="coolwarm",
        zorder=9,
    )
    plt.colorbar(s, label="Relative uncertainty")
    plt.clim(0, 1)
    plt.xlabel("Multiplicity")
    plt.ylabel(ylabel)
    plt.yscale(yscale)
    plt.title(title)

    plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.grid(which="both", linestyle="--", alpha=0.5)

    if outpath:
        Path(outpath).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(outpath, bbox_inches="tight", dpi=300, format="pdf")

    plt.show()


def plot_host_namd(
    df: pd.DataFrame,
    hostname: str,
    kind: str,
    Npt: int = 100000,
    threshold: int = 1000,
    outpath: str = None,
):
    """
    Plot the NAMD for a given host.

    Parameters
    ----------
    df : pd.DataFrame
        The NAMD database. If None, the function reloads the database from disk.
    hostname : str
        The hostname to plot.
    kind : str
        Which type of NAMD to plot. One of 'rel' (relative NAMD) or 'abs' (absolute NAMD).
    Npt : int
        Number of Monte Carlo samples.
    threshold : int
        Minimum number of valid samples required.

    Returns
    -------
    None
    """
    # Task 1: reload database
    if df is None:
        logger.info("Reloading the database")
        df = pd.read_csv(os.path.join(ROOT, "data", "exo_namd.csv"))
        logger.info("Database reloaded")

    # Task 1: sample the NAMD for a given host
    logger.info(f"Selecting the host: {hostname}")
    host = df[df["hostname"] == hostname]
    logger.info("Host selected")

    logger.info("Computing the Monte Carlo relative NAMD")
    retval = solve_namd_mc(
        host=host,
        kind=f"{kind}",
        Npt=Npt,
        threshold=threshold,
        full=True,
    )
    logger.info("Values computed")

    # Task 2: plot the NAMD for a given host
    logger.info("Plotting the relative NAMD distribution")
    simple_plot(
        df=retval,
        kind=f"{kind}",
        title=hostname,
        which="namd",
        scale="log",
        outpath=outpath,
    )
    logger.info("Plot done")


def plot_sample_namd(df: pd.DataFrame, title: str):
    """
    Plot the sample NAMD against the multiplicity.

    If df is None, the function reloads the database from disk.

    Parameters
    ----------
    df : pd.DataFrame
        The NAMD database.
    title : str
        The title of the plot.

    Returns
    -------
    None
    """
    # Task 1: reload database
    if df is None:
        logger.info("Reloading the database")
        df = pd.read_csv(os.path.join(ROOT, "data", "exo_namd.csv"))
        logger.info("Database reloaded")

    # Task 2: plot the sample NAMD
    logger.info("Plotting the NAMD vs. multiplicity")
    pop_plot(
        df=df.groupby("hostname").apply(
            lambda g: g.select_dtypes(exclude=["object"]).mean(),
        ),
        kind="rel",
        title=title,
        which="namd",
        yscale="log",
        xoffs=True,
    )
    logger.info("Plot done")
