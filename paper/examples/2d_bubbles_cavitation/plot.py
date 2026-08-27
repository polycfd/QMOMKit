#!/usr/bin/env python3
"""Plot the paper moment comparison for the polydisperse-cavitation case."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


QUANTITIES = [
    ("M_00", r"$M_{00}$"),
    ("M_10/M_00", r"$M_{10}/M_{00}$"),
    ("M_01/M_00", r"$M_{01}/M_{00}$"),
    ("M_11/M_00", r"$M_{11}/M_{00}$"),
    ("M_20/M_00", r"$M_{20}/M_{00}$"),
    ("M_02/M_00", r"$M_{02}/M_{00}$"),
]


def main() -> None:
    data = np.genfromtxt(Path("output") / "moment_comparison.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    figure, axes = plt.subplots(2, 3, figsize=(10, 5.5), constrained_layout=True)
    for axis, (quantity, label) in zip(axes.flat, QUANTITIES):
        selected = data[data["quantity"] == quantity]
        axis.plot(selected["time"], selected["qmomkit"], label="QMOMKit")
        axis.plot(selected["time"], selected["monte_carlo"], "o", fillstyle="none", label="Monte Carlo")
        axis.set_title(label)
        axis.set_ylabel("quantity")
    figure.supxlabel("time")
    axes.flat[0].legend()
    figure.savefig("polydisperse_cavitation_moments.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
