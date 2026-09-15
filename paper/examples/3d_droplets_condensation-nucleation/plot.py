#!/usr/bin/env python3
"""Plot the paper moment comparison for the condensing-droplet case.

QMOMKit and Monte Carlo moments are read from the native histories the run
writes; the compared quantities are derived here.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# Phase-space order is (r, w, T).
MOMENTS = ["M[0,0,0]", "M[3,0,0]", "M[2,0,0]", "M[0,1,0]", "M[0,0,1]"]

QUANTITIES = [
    (r"$M_{000}$", lambda m: m["M[0,0,0]"]),
    (r"$M_{300}/M_{200}$", lambda m: m["M[3,0,0]"] / m["M[2,0,0]"]),
    (r"$M_{010}/M_{000}$", lambda m: m["M[0,1,0]"] / m["M[0,0,0]"]),
    (r"$M_{001}/M_{000}$", lambda m: m["M[0,0,1]"] / m["M[0,0,0]"]),
]


def read_native_history(path: Path, labels: list[str]) -> dict[str, np.ndarray]:
    """Read selected raw moments from a QMOMKit-native moment history.

    The header is a commented, column-aligned label row, so the comment marker
    is not at the start of the line.
    """

    with path.open(encoding="utf-8") as handle:
        header = handle.readline().strip().lstrip("#").split()

    if not header or header[0] != "time":
        raise ValueError(f"{path} does not contain a moment-history header")

    values = np.loadtxt(path, comments="#", ndmin=2)
    return {label: values[:, header.index(label)] for label in ["time", *labels]}


def main() -> None:
    output = Path("output")
    qmomkit = read_native_history(output / "moment_history.txt", MOMENTS)
    monte_carlo = read_native_history(output / "moment_history_mc.txt", MOMENTS)
    figure, axes = plt.subplots(2, 2, figsize=(8, 6), constrained_layout=True)
    for axis, (label, evaluate) in zip(axes.flat, QUANTITIES):
        axis.plot(qmomkit["time"], evaluate(qmomkit), label="QMOMKit")
        axis.plot(monte_carlo["time"], evaluate(monte_carlo), "o", fillstyle="none", label="Monte Carlo")
        axis.set_title(label)
        axis.set_ylabel("quantity")
    figure.supxlabel("time")
    axes.flat[0].legend()
    figure.savefig("condensing_droplets_moments.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
