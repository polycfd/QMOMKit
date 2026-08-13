#!/usr/bin/env python3
"""Plot the paper moment comparison for the non-spherical-particle case."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


QUANTITIES = [
    ("M_00000", r"$M_{00000}$"),
    ("M_30000/M_20000", r"$d_{32}=M_{30000}/M_{20000}$"),
    ("M_01000/M_00000", r"$\overline{\chi}=M_{01000}/M_{00000}$"),
    ("M_30100/M_00000", r"$P_x=M_{30100}/M_{00000}$"),
    ("M_30010/M_00000", r"$P_y=M_{30010}/M_{00000}$"),
    ("M_30001/M_00000", r"$P_z=M_{30001}/M_{00000}$"),
    ("Theta_x", r"$\Theta_x$"),
    ("Theta_y", r"$\Theta_y$"),
    ("Theta_z", r"$\Theta_z$"),
]


def main() -> None:
    data = np.genfromtxt(Path("output") / "moment_comparison.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    figure, axes = plt.subplots(3, 3, figsize=(11, 8), constrained_layout=True)
    for axis, (quantity, label) in zip(axes.flat, QUANTITIES):
        selected = data[data["quantity"] == quantity]
        axis.plot(selected["time"], selected["qmomkit"], label="QMOMKit")
        axis.plot(selected["time"], selected["monte_carlo"], "o", fillstyle="none", label="Monte Carlo")
        axis.set_title(label)
        axis.set_ylabel("quantity")
    figure.supxlabel("time")
    axes.flat[0].legend()
    figure.savefig("non_spherical_particles_moments.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
