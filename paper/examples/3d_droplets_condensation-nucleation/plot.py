#!/usr/bin/env python3
"""Plot the paper moment comparison for the condensing-droplet case."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


QUANTITIES = [
    ("M_000", r"$M_{000}$"),
    ("M_300/M_200", r"$M_{300}/M_{200}$"),
    ("M_010/M_000", r"$M_{010}/M_{000}$"),
    ("M_001/M_000", r"$M_{001}/M_{000}$"),
]


def main() -> None:
    data = np.genfromtxt(Path("output") / "moment_comparison.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    figure, axes = plt.subplots(2, 2, figsize=(8, 6), constrained_layout=True)
    for axis, (quantity, label) in zip(axes.flat, QUANTITIES):
        selected = data[data["quantity"] == quantity]
        axis.plot(selected["time"], selected["qmomkit"], label="QMOMKit")
        axis.plot(selected["time"], selected["monte_carlo"], "o", fillstyle="none", label="Monte Carlo")
        axis.set_title(label)
        axis.set_ylabel("quantity")
    figure.supxlabel("time")
    axes.flat[0].legend()
    figure.savefig("condensing_droplets_moments.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
