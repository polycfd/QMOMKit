#!/usr/bin/env python3
"""Plot the CQMOM-versus-Monte-Carlo comparison for reactive spheres."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


QUANTITIES = [
    ("M_00000", r"$M_{00000}$"),
    ("M_30000/M_20000", r"$d_{32}$"),
    ("mean_diameter_m", r"$\overline{d}$"),
    ("diameter_variance_m2", r"$\mathrm{Var}(d)$"),
    ("mean_u_m_per_s", r"$\overline{u}$"),
    ("velocity_u_variance_m2_per_s2", r"$\mathrm{Var}(u)$"),
    ("mean_v_m_per_s", r"$\overline{v}$"),
    ("velocity_v_variance_m2_per_s2", r"$\mathrm{Var}(v)$"),
    ("mean_w_m_per_s", r"$\overline{w}$"),
    ("velocity_w_variance_m2_per_s2", r"$\mathrm{Var}(w)$"),
    ("mean_temperature_kK", r"$\overline{T}$"),
    ("temperature_variance_kK2", r"$\mathrm{Var}(T)$"),
    ("diameter_temperature_covariance_m_kK", r"$\mathrm{Cov}(d,T)$"),
    ("M_31000/M_00000", r"$M_{31000}/M_{00000}$"),
]


def main() -> None:
    data = np.genfromtxt(Path("output") / "moment_comparison.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    figure, axes = plt.subplots(4, 4, figsize=(12, 9), constrained_layout=True)
    for axis, (quantity, label) in zip(axes.flat, QUANTITIES):
        selected = data[data["quantity"] == quantity]
        axis.plot(selected["time"], selected["qmomkit"], label="QMOMKit")
        axis.plot(selected["time"], selected["monte_carlo"], "o", fillstyle="none", label="Monte Carlo")
        axis.set_title(label)
        axis.set_ylabel("quantity")
    for axis in axes.flat[len(QUANTITIES) :]:
        axis.set_visible(False)
    figure.supxlabel("time (s)")
    axes.flat[0].legend()
    figure.savefig("reactive_spherical_particles_moments.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
