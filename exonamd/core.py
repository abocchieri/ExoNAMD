import numpy as np
import pandas as pd


def compute_amdk(m, e, di, a):
    return m * np.sqrt(a) * (1 - np.sqrt(1 - e**2) * np.cos(np.deg2rad(di)))


def compute_namd(amdk, m, sqrt_a):
    return np.sum(amdk) / np.sum(m * sqrt_a)


# def namd_loop(df: pd.DataFrame, Npt: int = 100000):
#     df_namd = pd.DataFrame(
#         columns=[
#             "hostname",
#             "namd",
#             "namd_err1",
#             "namd_err2",
#             "N",
#             "flag",
#             "relative_uncertainty",
#         ]
#     )
#     for k, name in enumerate(df.hostname.unique()):
#         df_namd.loc[k, "hostname"] = name
#         df_namd.loc[k, "flag"] = df[df["hostname"] == name]["flag"].values
#         retval = compute_namd_mcmc(name, df, Npt=Npt, do_plot=False)
#         df_namd.loc[k, "namd"] = retval["q50"]
#         df_namd.loc[k, "namd_err1"] = retval["q84"] - retval["q50"]
#         df_namd.loc[k, "namd_err2"] = retval["q50"] - retval["q16"]
#         df_namd.loc[k, "N"] = retval["N"]
#         df_namd.loc[k, "relative_uncertainty"] = (
#             (retval["q84"] - retval["q16"]) / 2 / retval["q50"]
#         )

#     # # drop nans
#     df_namd = df_namd.dropna()
#     # drop where relative uncertainty is greater than 1
#     remove_large_uncertainty = df_namd["relative_uncertainty"] > 1
#     print(len(df_namd[remove_large_uncertainty]))
#     df_namd = df_namd[~remove_large_uncertainty]

#     return df_namd
