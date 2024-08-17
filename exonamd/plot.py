import numpy as np
import matplotlib.pyplot as plt


def simple_plot(
    data,
    kind,
    title="",
    which="namd",
    ylabel="Frequency",
    scale="linear",
    bins=50,
):
    samples = data[f"{which}_{kind}_mc"]
    q50 = data[f"{which}_{kind}_q50"]
    q16 = data[f"{which}_{kind}_q16"]
    q84 = data[f"{which}_{kind}_q84"]

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
    plt.show()


def pop_plot(df, kind, title="", which="namd", yscale="log", xoffs=False):
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

    plt.show()
