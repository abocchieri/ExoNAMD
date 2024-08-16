import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def compute_amd(m: np.array, e: np.array, di: np.array, a: np.array) -> np.array:
    return m * np.sqrt(a) * (1 - np.sqrt(1 - e**2) * np.cos(np.deg2rad(di)))


def compute_namd(amd: np.array, m: np.array, sqrt_a: np.array) -> float:
    return np.sum(amd) / np.sum(m * sqrt_a)


def compute_namd_mcmc(
    userid: str,
    df: pd.DataFrame,
    Npt: int,
    do_plot: bool = True,
    seed: int = 42,
    return_samples: bool = False,
):
    np.random.seed(seed)

    user_df = df[df["hostname"] == userid]

    pl_df_mcmc_samples = pd.DataFrame(
        columns=np.array(
            [
                [
                    f"pl_bmasse_{pl_name}",
                    f"pl_orbsmax_{pl_name}",
                    f"pl_orbeccen_{pl_name}",
                    f"pl_relincl_{pl_name}",
                ]
                for pl_name in user_df["pl_name"]
            ]
        ).flatten(),
        index=range(Npt),
    )

    for pl_name in user_df["pl_name"]:
        pl_df = user_df[user_df["pl_name"] == pl_name]

        pl_df_mcmc_samples[f"pl_bmasse_{pl_name}"] = np.random.normal(
            loc=pl_df["pl_bmasse"],
            scale=0.5 * (pl_df["pl_bmasseerr1"] - pl_df["pl_bmasseerr2"]),
            size=Npt,
        )

        pl_df_mcmc_samples[f"pl_orbsmax_{pl_name}"] = np.random.normal(
            loc=pl_df["pl_orbsmax"],
            scale=0.5 * (pl_df["pl_orbsmaxerr1"] - pl_df["pl_orbsmaxerr2"]),
            size=Npt,
        )

        pl_df_mcmc_samples[f"pl_orbeccen_{pl_name}"] = np.random.normal(
            loc=pl_df["pl_orbeccen"],
            scale=0.5 * (pl_df["pl_orbeccenerr1"] - pl_df["pl_orbeccenerr2"]),
            size=Npt,
        )

        pl_df_mcmc_samples[f"pl_relincl_{pl_name}"] = np.random.normal(
            loc=pl_df["pl_relincl"],
            scale=0.5 * (pl_df["pl_relinclerr1"] - pl_df["pl_relinclerr2"]),
            size=Npt,
        )

        pl_df_mcmc_samples[f"good_idx_{pl_name}"] = (
            (pl_df_mcmc_samples[f"pl_orbeccen_{pl_name}"] >= 0)
            & (pl_df_mcmc_samples[f"pl_orbeccen_{pl_name}"] <= 1)
            & (np.abs(pl_df_mcmc_samples[f"pl_relincl_{pl_name}"]) >= 0)
            & (np.abs(pl_df_mcmc_samples[f"pl_relincl_{pl_name}"]) <= 180)
            & (pl_df_mcmc_samples[f"pl_orbsmax_{pl_name}"] >= 0)
            & (pl_df_mcmc_samples[f"pl_bmasse_{pl_name}"] >= 0)
        )

        pl_df_mcmc_samples[f"amd_{pl_name}"] = compute_amd(
            pl_df_mcmc_samples[f"pl_bmasse_{pl_name}"],
            pl_df_mcmc_samples[f"pl_orbeccen_{pl_name}"],
            pl_df_mcmc_samples[f"pl_relincl_{pl_name}"],
            pl_df_mcmc_samples[f"pl_orbsmax_{pl_name}"],
        )

    # keep rows of pl_df_mcmc_samples where good_idx is True for all planets
    pl_df_mcmc_samples = pl_df_mcmc_samples[
        pl_df_mcmc_samples.filter(like="good_idx").all(axis=1)
    ]

    # for each row, sum the AMDs of all planets and store the result in a new column called "namd_num"
    pl_df_mcmc_samples["namd_num"] = pl_df_mcmc_samples.filter(like="amd").sum(axis=1)

    # for each row, compute the sum of: (mass of each planet * sqrt(semi-major axis of each planet))
    pl_df_mcmc_samples["namd_den"] = np.sum(
        [
            pl_df_mcmc_samples[f"pl_bmasse_{pl_name}"]
            * np.sqrt(pl_df_mcmc_samples[f"pl_orbsmax_{pl_name}"])
            for pl_name in user_df["pl_name"]
        ],
        axis=0,
    )

    # compute the NAMD
    pl_df_mcmc_samples["namd"] = (
        pl_df_mcmc_samples["namd_num"] / pl_df_mcmc_samples["namd_den"]
    )

    if len(pl_df_mcmc_samples) < Npt // 20:
        return {
            "q16": np.nan,
            "q50": np.nan,
            "q84": np.nan,
            "N": user_df["sy_pnum"].values[0],
        }

    quantiles = pl_df_mcmc_samples["namd"].quantile([0.16, 0.5, 0.84])
    if do_plot:
        print("Plotting...")
        plt.figure()
        # plt.hist(pl_df_mcmc_samples["namd"], bins=50)
        # use log-spaced bins instead of linear-spaced bins
        # plt.hist(pl_df_mcmc_samples["namd"], bins=np.logspace(np.log10(pl_df_mcmc_samples["namd"].min()), np.log10(pl_df_mcmc_samples["namd"].max()), 50))
        log_namd = np.log10(pl_df_mcmc_samples["namd"])
        log_quantiles = log_namd.quantile([0.16, 0.5, 0.84])
        weights = np.ones_like(log_namd) / len(log_namd)
        plt.hist(log_namd, bins=50, alpha=0.5, weights=weights)
        plt.vlines(
            log_quantiles,
            0,
            plt.ylim()[1],
            color=["red", "black", "red"],
            linestyles="dashed",
        )
        plt.xlabel(rf"log(NAMD)")
        plt.ylabel("Relative frequency")
        title = (
            f"{userid}: "
            + rf"log(NAMD) = ${log_quantiles[0.5]:.2f}^{{+{log_quantiles[0.84]-log_quantiles[0.5]:.2f}}}_{{-{log_quantiles[0.5]-log_quantiles[0.16]:.2f}}}$"
        )
        plt.title(title)
        plt.show()

    return_dict = {
        "q16": quantiles[0.16],
        "q50": quantiles[0.5],
        "q84": quantiles[0.84],
        "N": user_df["sy_pnum"].values[0],
    }

    if return_samples:
        return_dict.update({"samples": pl_df_mcmc_samples})

    return return_dict


def namd_loop(df: pd.DataFrame, Npt: int = 100000):
    df_namd = pd.DataFrame(
        columns=[
            "hostname",
            "namd",
            "namd_err1",
            "namd_err2",
            "N",
            "flag",
            "relative_uncertainty",
        ]
    )
    for k, name in enumerate(df.hostname.unique()):
        df_namd.loc[k, "hostname"] = name
        df_namd.loc[k, "flag"] = df[df["hostname"] == name]["flag"].values
        retval = compute_namd_mcmc(name, df, Npt=Npt, do_plot=False)
        df_namd.loc[k, "namd"] = retval["q50"]
        df_namd.loc[k, "namd_err1"] = retval["q84"] - retval["q50"]
        df_namd.loc[k, "namd_err2"] = retval["q50"] - retval["q16"]
        df_namd.loc[k, "N"] = retval["N"]
        df_namd.loc[k, "relative_uncertainty"] = (
            (retval["q84"] - retval["q16"]) / 2 / retval["q50"]
        )

    # # drop nans
    df_namd = df_namd.dropna()
    # drop where relative uncertainty is greater than 1
    remove_large_uncertainty = df_namd["relative_uncertainty"] > 1
    print(len(df_namd[remove_large_uncertainty]))
    df_namd = df_namd[~remove_large_uncertainty]

    return df_namd
