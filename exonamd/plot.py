import numpy as np
import matplotlib.pyplot as plt


def simple_plot(
    data,
    kind,
    title,
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
